# Day 4 Commands

This practical introduces the strainphlan tool to view the strains present for each species, and then further explains the interpretation and analysis of results from the HUMAnN tool from the previous day.

---

## Workflow Overview

Run the strainphlan → Analyze the strains of one species → Interpretation of results from the HUMAnN tool and its analysis → Case Study → Human Microbiome Databases

---
## 🔹 STEP 1: Run Strainphlan on the SAM files generated from the Metaphlan and Show Strains of One species


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

---
