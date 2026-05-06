# Day 2 Commands

This practical covers downstream microbiome analysis using processed taxonomic profiles, including abundance normalization, alpha diversity, beta diversity, and differential abundance statistics.

---

# Workflow Overview

Processed Taxonomic Profile → Normalization → Diversity Analysis → Statistical Testing → Visualization

---

## 🔹 STEP 1: Processing Taxonomic Profile from SPINGO and Data Normalization

### Purpose
Merge all SPINGO output files into a single abundance matrix and normalize the species abundance profiles for downstream microbiome analysis.

---

### Merge SPINGO Outputs into Abundance Matrix

### create a text file containing all the names of output files
```bash
ls *.spingo.out.txt > all_spingo_outputs.txt
```

### Now give all_spingo_outputs.txt to a perl script which merge all the samples and create abundance profile (matrix)
```bash
perl /path_of_the_file/create_species_matrix.pl all_spingo_outputs.txt > Abundance_Profile.txt
```

### Now see the brief view of the (matrix) created
```bash
head Abundance_Profile.txt
```

---

## transfer Abundance_Profile.txt to laptop/local system

### On your local system open command prompt (cmd)

### Then open the new tab using Windows PowerShell and go in Downloads directory (cd Downloads/)
```bash
cd Downloads/
```

### transfer abundance profile from server to local system
```bash
scp -r user1@192.168.22.173:/Path_to_the_abundance_profile/Abundance_Profile.txt .
```

---

# 🔹 STEP 2: Data Normalization in RStudio

## Purpose
Import abundance matrix and metadata into RStudio and normalize the abundance profile.

---

## Open the RStudio on your laptop/local system

### Now Import the Abundance_Profile.txt using below command in RStudio
```r
species_matrix <- read.delim(
  "C:/Users/ompra/Downloads/Abundance_Profile.txt",
  row.names = 1
)
```

---

## Import the metadata
```r
library(readxl)

metadata_df <- data.frame(
  read_excel("C:/Users/ompra/Downloads/metadata_df.xlsx")
)
```

### Add rownames to the metadata dataframe
```r
rownames(metadata_df) <- metadata_df$sampleID
```

---

## Now we need to have same rownames in the species_matrix and metadata_df (patients matched)

### Replace the '_' and everything after that with nothing ''
```r
rownames(species_matrix) <- sub("_.*","",rownames(species_matrix))
```

### Replace the starting '0' with nothing ''. '^' is to mention at the beginning of the name
```r
rownames(species_matrix) <- sub("^0","",rownames(species_matrix))
```

### Now Add 'MT' at the beginning of every name.
```r
rownames(species_matrix) <- gsub("^","MT",rownames(species_matrix))
```

---

## now confirm all the sample ids are same in species profile and metadata dataframe
```r
setdiff(rownames(species_matrix), rownames(metadata_df))

setdiff(rownames(metadata_df), rownames(species_matrix))
```

---

## Normalize the species profile (RowSum Normalization) and order the rownames in both dataframe
```r
species_matrix_norm <- species_matrix / rowSums(species_matrix)

species_matrix_norm <- species_matrix_norm[rownames(metadata_df),]
```

---

# 🔹 STEP 3: Alpha Diversity and Visualization Techniques (SHANNON + PIELOU)

## Purpose
Compute microbial diversity and evenness between Control and UC samples.

---

## install necessary packages
```r
install.packages("vegan")    # required for computing alpha diversity

install.packages("vioplot")  # required for violin plot visualization
```

---

## load necessary packages
```r
library(vegan)

library(vioplot)
```

---

## INPUT
- metadata_df:
  - sampleID
  - study_condition (Control / UC)

- species_matrix_norm:
  - rows = sampleID
  - columns = species

---

## make sure the rownames of species profile and metadata are same and in same order
```r
all(rownames(species_matrix_norm) %in% rownames(metadata_df))
```

---

## computing SHANNON INDEX
```r
shannon_index <- diversity(
  species_matrix_norm,
  index = "shannon"
)
```

---

## computing PIELOU INDEX
```r
calculate_pielou <- function(mat) {

  H <- diversity(mat, index = "shannon")

  S <- rowSums(mat > 0)

  ifelse(S > 1, H / log(S), NA)
}

pielou_index <- calculate_pielou(species_matrix_norm)
```

---

## create a dataframe for storing alpha diversity results
```r
alpha_diversity_df <- data.frame(
  study_condition = metadata_df$study_condition,
  shannon = shannon_index,
  pielou = pielou_index
)
```

---

# 🔹 STEP 4: Visualization of Alpha Diversity

## SHANNON INDEX

### 1) boxplot
```r
boxplot(
  alpha_diversity_df$shannon[
    alpha_diversity_df$study_condition == "Control"
  ],

  alpha_diversity_df$shannon[
    alpha_diversity_df$study_condition == "UC"
  ],

  names = c("Control", "UC"),

  col = c("turquoise3", "pink2"),

  ylab = "Shannon Index",

  outline = FALSE
)
```

---

### 2) violin plot
```r
vioplot(
  alpha_diversity_df$shannon[
    alpha_diversity_df$study_condition == "Control"
  ],

  alpha_diversity_df$shannon[
    alpha_diversity_df$study_condition == "UC"
  ],

  names = c("Control", "UC"),

  col = c("turquoise3", "pink2"),

  ylab = "Shannon Index"
)
```

---

## wilcox test for showing significance between the Control and UC group
```r
wilcox.test(
  alpha_diversity_df$shannon[
    alpha_diversity_df$study_condition == "Control"
  ],

  alpha_diversity_df$shannon[
    alpha_diversity_df$study_condition == "UC"
  ]
)
```

---

## PIELOU INDEX

### 1) boxplot
```r
boxplot(
  alpha_diversity_df$pielou[
    alpha_diversity_df$study_condition == "Control"
  ],

  alpha_diversity_df$pielou[
    alpha_diversity_df$study_condition == "UC"
  ],

  names = c("Control", "UC"),

  col = c("turquoise3", "pink2"),

  ylab = "Pielou Index",

  outline = FALSE
)
```

---

### 2) violin plot
```r
vioplot(
  alpha_diversity_df$pielou[
    alpha_diversity_df$study_condition == "Control"
  ],

  alpha_diversity_df$pielou[
    alpha_diversity_df$study_condition == "UC"
  ],

  names = c("Control", "UC"),

  col = c("turquoise3", "pink2"),

  ylab = "Pielou Index"
)
```

---

## wilcox test for showing significance between the Control and UC group
```r
wilcox.test(
  alpha_diversity_df$pielou[
    alpha_diversity_df$study_condition == "Control"
  ],

  alpha_diversity_df$pielou[
    alpha_diversity_df$study_condition == "UC"
  ]
)
```

---


# 🔹 STEP 5: Beta Diversity and Visualization Techniques

## Purpose
Evaluate differences in microbial community composition between groups.

---


---

# 🔹 STEP 6: Differential Abundance Statistics (Species associated with IBD Disease patients and Control Patients)

## separate the control and disease metadata, and then separate the species profile using their metadata
```r
control_metadata <- metadata_df[
  metadata_df$study_condition == "Control",
]

disease_metadata <- metadata_df[
  metadata_df$study_condition == "UC",
]
```

---

## subset species profiles
```r
control_species_profile <- species_matrix_norm[
  rownames(control_metadata),
]

disease_species_profile <- species_matrix_norm[
  rownames(disease_metadata),
]
```

---

## Now use wilcox.test for one species (Coprococcus_catus is Health Associated Core Species (from previous references))
```r
wilcox.test(
  disease_species_profile$Coprococcus_catus,

  control_species_profile$Coprococcus_catus
)
```

---

## mean abundance of Coprococcus_catus in disease samples
```r
mean(disease_species_profile$Coprococcus_catus)
```

---

## mean abundance of Coprococcus_catus in control samples
```r
mean(control_species_profile$Coprococcus_catus)
```

---

## See the difference of mean abundance of Coprococcus_catus in disease and control
## (disease minus control >> If -ve then Control associated, If +ve then Disease associated)
```r
mean(disease_species_profile$Coprococcus_catus) - mean(control_species_profile$Coprococcus_catus)
```

---

# 🔹 STEP 7: Use wilcox.test for all the species using wilcox batch function

## install required package
```r
install.packages("dplyr")
```

---

## load required package
```r
library(dplyr)
```

---

## Wilcoxon batch function
```r
wilcox_batch = function(x,y)
{
  p_array <- NULL
  type_array <- NULL
  mean1_array <- NULL
  mean2_array <- NULL

  x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,]
  y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,]

  z <- intersect(rownames(x),rownames(y))

  for(i in 1:length(z))
  {
    p_array[i] <- wilcox.test(
      as.numeric(x[z[i],]),
      as.numeric(y[z[i],])
    )$p.value

    type_array[i] <- ifelse(
      mean(as.numeric(x[z[i],]),na.rm=TRUE) >
      mean(as.numeric(y[z[i],]),na.rm=TRUE),

      1,

      ifelse(
        mean(as.numeric(x[z[i],]),na.rm=TRUE) <
        mean(as.numeric(y[z[i],]),na.rm=TRUE),

        -1,

        0
      )
    )

    mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE)

    mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE)
  }

  out <- as.data.frame(
    cbind(
      p_array,
      type_array,
      p.adjust(p_array,method="fdr"),
      mean1_array,
      mean2_array
    )
  )

  rownames(out) <- z

  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)))

  out <- data.frame(t(out))

  colnames(out)[3] <- "p_adjust"

  out <- out %>%
    mutate(
      directionality = ifelse(
        type_array == 0,
        0,

        ifelse(
          p_adjust <= 0.05,
          3 * sign(type_array),

          ifelse(
            p_array <= 0.05,
            2 * sign(type_array),

            type_array
          )
        )
      )
    )

  return(out)
}
```

---

## Run the function for our data
```r
wilcox_result <- wilcox_batch(
  t(disease_species_profile),
  t(control_species_profile)
)
```

---

# 🔹 STEP 8: Volcano Plot Visualization

## install required packages
```r
install.packages("ggplot2")

install.packages("ggrepel")
```

---

## load required packages
```r
library(ggplot2)

library(ggrepel)
```

---

## Add the log2 fold change
## (value 1 = mean abundance in disease is double the mean abundance in control)
```r
wilcox_result$log2_fold_change <- log2(
  wilcox_result$mean1_array /
  wilcox_result$mean2_array
)
```

---

## Add the significance column based on p.adjust i.e q-value
```r
wilcox_result$significance <- ifelse(
  wilcox_result$p_adjust <= 0.05,
  "significant",
  "non-significant"
)
```

---

## For proper visualization we can use negative log of q-value
## (Higher the value >> more is the significance)
```r
wilcox_result$logQValue <- -log10(
  wilcox_result$p_adjust
)
```

---

## add the color code based on the log2 fold change
## (Only species with double mean abundance change will be colored)
```r
wilcox_result$color_code <- ifelse(
  wilcox_result$log2_fold_change > 1,
  "Enriched",

  ifelse(
    wilcox_result$log2_fold_change < -1,
    "Depleted",

    "Moderate Change"
  )
)
```

---

## Filter the dataframe where only significant species are present (q-value <= 0.05)
```r
wilcox_result_filt <- wilcox_result[
  which(wilcox_result$significance == "significant"),
]
```

---

## Using q-value <= 0.05, we have lot of species,
## so Visualization will look congested
## therefore filter species having q-value <= 0.001
```r
wilcox_result_filt2 <- wilcox_result_filt[
  which(wilcox_result_filt$p_adjust <= 0.001),
]
```

---

## Now plot the Volcano plot
```r
pdf("Wilcox_results.pdf", height = 8, width = 15)

ggplot(
  wilcox_result_filt2,
  aes(x = log2_fold_change, y = logQValue)
) +

  geom_point(
    aes(color = color_code),
    size = 3
  ) +

  geom_text_repel(
    data = wilcox_result_filt2,
    aes(label = rownames(wilcox_result_filt2)),
    size = 4,
    max.overlaps = 45
  ) +

  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey",
    linewidth = 0.7
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +

  scale_color_manual(
    values = c(
      "Enriched" = "royalblue3",
      "Depleted" = "palevioletred2",
      "Moderate Change" = "gray3"
    )
  ) +

  theme_minimal() +

  labs(
    x = "Log2 Fold Change (IBD/Control)",
    y = "-logQValue"
  ) +

  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),

    axis.text = element_text(size = 12),

    axis.title = element_text(size = 12)
  )

dev.off()
```

---

# 💡 Key Notes

- Normalize abundance profiles before diversity analysis
- Use Wilcoxon test for non-parametric microbiome comparisons
- Use FDR-adjusted p-values (q-values) for multiple testing correction
- Volcano plots help identify biologically important taxa
