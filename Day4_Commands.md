# Day 4 Commands

This practical introduces the strainphlan tool for viewing the strains present in each species, and then further explains the interpretation and analysis of results from the HUMAnN tool from the previous day.

---

## Workflow Overview

Run the strainphlan → Analyze the strains of one species → Interpretation of results from the HUMAnN tool and its analysis → Case Study → Human Microbiome Databases

---
## 🔹 Session 1: Run Strainphlan on the SAM files generated from the Metaphlan and Show Strains of One species


### Extracts clade-specific consensus marker genes from MetaPhlAn alignment files for each sample

```bash
mkdir -p consensus_marker
sample2markers.py -i sams/*.sam.bz2 -o consensus_marker -n 20
```

### Identifies clades/species that have sufficient marker representation across at least 10 samples

```bash
# print everything from the below command in print_clades_only.tsv (species names and also other things)
strainphlan -s consensus_marker/*.pkl -o consensus_marker --print_clades_only --marker_in_n_samples 10 > consensus_marker/print_clades_only.tsv

# extract only species names 
grep -o "s__[^:]*" consensus_marker/print_clades_only.tsv > consensus_marker/clades.txt
```

### Extracts the MetaPhlAn reference marker genes corresponding to one clade/species

```bash
mkdir -p db_marker output

# extract only species names 
extract_markers.py -c clade_name -o db_marker
```

### Performs strain-level phylogenetic reconstruction for the selected clade/species (using sample consensus markers and reference markers)

```bash
strainphlan -s consensus_marker/*.pkl -m db_marker/clade_name.fna -o output -c clade_name --marker_in_n_samples 10 --n 20
```


### Run strainphlan along with Metaphlan in one batch script (Only for your reference)

#### store this script in one bash file (.sh) and then run help option to know how it runs
```bash
#!/bin/bash

# Function to display help
show_help() {
    echo "Usage: $0 <directory> <paired|single> [steps]"
    echo ""
    echo "This script processes a directory with fastq files (paired-end or single-end) through a series of steps:"
    echo "1. Metaphlan (MetaPhlAn 3) for each file or pair of forward and reverse files"
    echo "2. Generate consensus marker."
    echo "3. Print clades for further process"
    echo "4. Extract markers for each clade and run StrainPhlAn (v3)"
    echo ""
    echo "Examples:"
    echo "  Run all steps: bash script.sh /path/to/directory paired 1 2 3 4"
    echo "  Run specific steps: bash script.sh /path/to/directory single 1 3"
    echo ""
    echo "Options:"
    echo "  -h, --help    Show this help message and exit"
}

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

if [ $# -lt 3 ]; then
    echo "Usage: $0 <directory> <paired|single> [steps]"
    exit 1
fi

input_dir="$1"
file_type="$2"
shift 2
steps=($@)

cd "$input_dir" || { echo "Error: Unable to access directory $input_dir"; exit 1; }

step1_paired() {
    echo "Running Step 1 (paired-end): MetaPhlAn 3"
    mkdir -p bowtie2 sams
    for file in *_1.fastq.gz; do
        sample_name=$(basename "$file" _1.fastq.gz)
        pigz -p 20 -dc "${sample_name}_1.fastq.gz" "${sample_name}_2.fastq.gz" > "${sample_name}.merged.fastq"
        metaphlan "${sample_name}.merged.fastq" --input_type fastq -s "sams/${sample_name}.sam.bz2" --bowtie2out "bowtie2/${sample_name}.bowtie2.bz2" --nproc 20 -o "${sample_name}_profiled.txt"
        rm "${sample_name}.merged.fastq"
    done
    rm -rf bowtie2
}

step1_single() {
    echo "Running Step 1 (single-end): MetaPhlAn 3"
    mkdir -p bowtie2 sams

    for file in *.fastq.gz; do
        # Skip files ending in _1.fastq.gz or _2.fastq.gz
        if [[ "$file" == *_1.fastq.gz || "$file" == *_2.fastq.gz ]]; then
            continue
        fi

        sample_name=$(basename "$file" .fastq.gz)
       
        metaphlan "$file" --input_type fastq -s "sams/${sample_name}.sam.bz2" --bowtie2out "bowtie2/${sample_name}.bowtie2.bz2" --nproc 20 -o "${sample_name}_profiled.txt"
    done
    rm -rf bowtie2
}


step2() {
    echo "Running Step 2: sample2markers"
    mkdir -p consensus_marker
    sample2markers.py -i sams/*.sam.bz2 -o consensus_marker -n 20
}

step3() {
    echo "Running Step 3: strainphlan --print_clades_only"
    strainphlan -s consensus_marker/*.pkl -o consensus_marker --print_clades_only --marker_in_n_samples 10 > consensus_marker/print_clades_only.tsv

    grep -o "s__[^:]*" consensus_marker/print_clades_only.tsv > consensus_marker/clades.txt
    echo "Saved clades to consensus_marker/clades.txt"
}

step4() {
    echo "Running Step 4: extract_markers and strainphlan"
    clades_file="consensus_marker/clades.txt"
    if [ ! -f "$clades_file" ]; then
        echo "Error: $clades_file not found. Please run step 3 first."
        exit 1
    fi

    mkdir -p db_marker output
    while read -r clade; do
        echo "Processing clade: $clade"
        extract_markers.py -c "$clade" -o db_marker
        strainphlan -s consensus_marker/*.pkl -m "db_marker/${clade}.fna" -o output -c "$clade" --marker_in_n_samples 10 --n 20
    done < "$clades_file"
}

for step in "${steps[@]}"; do
    case $step in
        1)
            [ "$file_type" == "paired" ] && step1_paired || step1_single
            ;;
        2) step2 ;;
        3) step3 ;;
        4) step4 ;;
        *) echo "Invalid step: $step" ;;
    esac
done
```


---
## 🔹 Session 2: Functional 

## 🔹 Session 2: HUMAnN Output Processing 

### 1. Join Individual Sample Tables
```bash
mkdir -p combined_tables

# Join all gene families files into one matrix
humann_join_tables -i . -o combined_tables/all_genefamilies.tsv --file_name genefamilies

# Join all pathway abundance files into one matrix
humann_join_tables -i . -o combined_tables/all_pathabundance.tsv --file_name pathabundance

# Join all pathway coverage files into one matrix
humann_join_tables -i . -o combined_tables/all_pathcoverage.tsv --file_name pathcoverage

mkdir -p normalized_tables

# Renormalize the gene families table
humann_renorm_table -i combined_tables/all_genefamilies.tsv \
                    -o normalized_tables/all_genefamilies_cpm.tsv \
                    --units cpm --update-snames

# Renormalize the pathway abundance table
humann_renorm_table -i combined_tables/all_pathabundance.tsv \
                    -o normalized_tables/all_pathabundance_cpm.tsv \
                    --units cpm --update-snames

# Create a folder for your project locally
mkdir -p ~/Desktop/Metagenomics_Analysis

# Securely copy the normalized tables to your local machine
scp "sagar@192.168.3.254:~/human_outputs/normalized_tables/*.tsv" ~/Desktop/Metagenomics_Analysis/



---


## 🔹 Session 3: Case Study

We will use Whole Genome Sequencing (WGS) data to conduct a hands-on session based on a real case study.

![PEN-UC Study Design](images/PEN_UC_StudyDesign.png)
---

