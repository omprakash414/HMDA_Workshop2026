
# Day 2 Commands

This practical covers downstream microbiome analysis using processed taxonomic profiles, including abundance normalization, alpha diversity, beta diversity, and differential abundance statistics.

---

## Workflow Overview

Processed Taxonomic Profile → Normalization → Diversity Analysis → Statistical Testing → Visualization

---

## 🔹 STEP 1: Processing Taxonomic Profile from SPINGO and Data Normalization

### Purpose  
Merge all SPINGO output files into a single abundance matrix and normalize the species abundance profiles for downstream microbiome analysis.

---

### Merge SPINGO Outputs into Abundance Matrix

#### create a text file containing all the names of output files

```bash
ls *.spingo.out.txt > all_spingo_outputs.txt
```

#### Now give all_spingo_outputs.txt to a perl script which merge all the samples and create abundance profile (matrix)

```bash
perl /path_of_the_file/create_species_matrix.pl all_spingo_outputs.txt > Abundance_Profile.txt
```

#### Now see the brief view of the matrix created

```bash
head Abundance_Profile.txt
```

---

### Transfer Abundance_Profile.txt to laptop/local system

#### On your local system open command prompt (cmd)

#### Then open the new tab using Windows PowerShell and go in Downloads directory

```bash
cd Downloads/
```

#### transfer abundance profile from server to local system

```bash
scp -r user1@192.168.22.173:/Path_to_the_abundance_profile/Abundance_Profile.txt .
```

---

## 🔹 STEP 2: Data Normalization in RStudio

### Purpose

Import abundance matrix and metadata into RStudio and normalize the abundance profile.

---

### Open the RStudio on your laptop/local system

#### Import the Abundance_Profile.txt using below command in RStudio

```r
species_matrix <- read.delim("C:/Users/ompra/Downloads/Abundance_Profile.txt", row.names = 1)
```

---

### Import the metadata

```r
library(readxl)
metadata_df <- data.frame(read_excel("C:/Users/ompra/Downloads/metadata_df.xlsx"))
```

#### Add rownames to the metadata dataframe

```r
rownames(metadata_df) <- metadata_df$sample_id
```

---

### We need to have same rownames in species_matrix and metadata_df (patients matched)

#### Replace the '_' and everything after that with nothing ''

```r
rownames(species_matrix) <- sub("_.*", "", rownames(species_matrix))
```

#### Replace the starting '0' with nothing ''

#### '^' means beginning of the name

```r
rownames(species_matrix) <- sub("^0", "", rownames(species_matrix))
```

#### Add 'MT' at the beginning of every sample name

```r
rownames(species_matrix) <- gsub("^", "MT", rownames(species_matrix))
```

---

### Confirm all the sample ids are same in species profile and metadata dataframe

```r
setdiff(rownames(species_matrix), rownames(metadata_df))

setdiff(rownames(metadata_df), rownames(species_matrix))
```

---

### Normalize the species profile (RowSum Normalization)

### and order the rownames in both dataframe

```r
species_matrix_norm <- species_matrix / rowSums(species_matrix)
species_matrix_norm <- species_matrix_norm[rownames(metadata_df),]
```

---

## 🔹 STEP 3: Alpha Diversity and Visualization Techniques (SHANNON + PIELOU)

### Purpose

Compute microbial diversity and evenness between Control and UC samples.

---

### Install necessary packages

```r
install.packages("vegan")
install.packages("vioplot")
```

---

### Load necessary packages

```r
library(vegan)
library(vioplot)
```

---

### INPUT

* metadata_df
  * sample_id
  * study_condition (Control / UC)

* species_matrix_norm
  * rows = sample_id
  * columns = species

---

### Make sure the rownames of species profile and metadata are same and in same order

```r
all(rownames(species_matrix_norm) %in% rownames(metadata_df))
```

---

### Computing SHANNON INDEX

```r
shannon_index <- diversity(species_matrix_norm, index = "shannon")
```

---

### Computing PIELOU INDEX

```r
calculate_pielou <- function(mat) {
  H <- diversity(mat, index = "shannon")
  S <- rowSums(mat > 0)
  ifelse(S > 1, H / log(S), NA)
}

pielou_index <- calculate_pielou(species_matrix_norm)
```

---

### Create a dataframe for storing alpha diversity results

```r
alpha_diversity_df <- data.frame(study_condition = metadata_df$study_condition, shannon = shannon_index, pielou = pielou_index)
```

---

## 🔹 STEP 4: Visualization of Alpha Diversity

### Purpose

Visualize and compare alpha diversity between Control and UC samples.

---

### SHANNON INDEX

#### 1) boxplot

```r
boxplot(shannon ~ study_condition, data = alpha_diversity_df, col = c("turquoise3", "pink2"))
```

---

#### 2) violin plot

```r
vioplot(shannon ~ study_condition, data = alpha_diversity_df, col = c("turquoise3", "pink2"))
```

---

### Wilcox test for showing significance between the Control and UC group

```r
wilcox.test(
  alpha_diversity_df$shannon[alpha_diversity_df$study_condition == "Control"],
  alpha_diversity_df$shannon[alpha_diversity_df$study_condition == "UC"])
```

---

### PIELOU INDEX

#### 1) boxplot

```r
boxplot(pielou ~ study_condition, data = alpha_diversity_df, col = c("turquoise3", "pink2"))
```

---

#### 2) violin plot

```r
vioplot(pielou ~ study_condition, data = alpha_diversity_df, col = c("turquoise3", "pink2"))
```

---

### Wilcox test for showing significance between the Control and UC group

```r
wilcox.test(
  alpha_diversity_df$pielou[alpha_diversity_df$study_condition == "Control"],
  alpha_diversity_df$pielou[alpha_diversity_df$study_condition == "UC"])
```

---

## 🔹 STEP 5: Beta Diversity and Visualization Techniques

### Purpose

Evaluate differences in microbial community composition between groups using distance-based methods and PCoA visualization.

---

### Load required packages

```r
library(vegan)
library(ggplot2)
```

---

### Convert normalized abundance profile to matrix

```r
mat <- as.matrix(species_matrix_norm)
meta <- metadata_df
```

---

### Confirm sample order is same in metadata and abundance matrix

```r
all(rownames(mat) == meta$sample_id)
```

---

## Bray-Curtis Distance Based Beta Diversity

### Compute Bray-Curtis distance

```r
dist_bray <- vegdist(mat, method = "bray")
```

---

### Perform PCoA

```r
pcoa_bray <- cmdscale(dist_bray, k = 2, eig = TRUE)
```

---

### Create dataframe for plotting

```r
pcoa_bray_df <- as.data.frame(pcoa_bray$points)
colnames(pcoa_bray_df) <- c("PCoA1", "PCoA2")
pcoa_bray_df$condition <- meta$study_condition
```

---

### Plot Bray-Curtis PCoA

```r
ggplot(
  pcoa_bray_df,
  aes(PCoA1, PCoA2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCoA - Bray Curtis")
```

---

### PERMANOVA Statistics

```r
adonis_bray <- adonis2(dist_bray ~ study_condition, data = meta)
adonis_bray
```

---

## Euclidean Distance Based Beta Diversity

### Compute Euclidean distance

```r
dist_euc <- dist(mat, method = "euclidean")
```

---

### Perform PCoA

```r
pcoa_euc <- cmdscale(dist_euc, k = 2, eig = TRUE)
```

---

### Create dataframe for plotting

```r
pcoa_euc_df <- as.data.frame(pcoa_euc$points)
colnames(pcoa_euc_df) <- c("PCoA1", "PCoA2")
pcoa_euc_df$condition <- meta$study_condition
```

---

### Plot Euclidean PCoA

```r
ggplot(
  pcoa_euc_df,
  aes(PCoA1, PCoA2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCoA - Euclidean")
```

---

### PERMANOVA Statistics

```r
adonis_euc <- adonis2(dist_euc ~ study_condition, data = meta)
adonis_euc
```

---

## Kendall Distance Based Beta Diversity

### Load required package

```r
library(pcaPP)
```

---

### Compute Kendall based distance

```r
dist_kendall <- as.dist(1 - cor.fk(t(mat)) / 2)
```

---

### Perform PCoA

```r
pcoa_kendall <- cmdscale(dist_kendall, k = 2, eig = TRUE)
```

---

### Create dataframe for plotting

```r
pcoa_kendall_df <- as.data.frame(pcoa_kendall$points)
colnames(pcoa_kendall_df) <- c("PCoA1", "PCoA2")
pcoa_kendall_df$condition <- meta$study_condition
```

---

### Plot Kendall Distance PCoA

```r
ggplot(
  pcoa_kendall_df,
  aes(PCoA1, PCoA2, color = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCoA - Kendall Distance")
```

---

### PERMANOVA Statistics

```r
adonis_kendall <- adonis2(dist_kendall ~ study_condition, data = meta)
adonis_kendall
```

---

## Advanced Beta Diversity Visualization

### Load required package

```r
library(ade4)
```

---

### Plot confidence ellipses with group separation

```r
pcoa_coords <- pcoa_kendall_df[,c("PCoA1", "PCoA2")]
grp <- as.factor(pcoa_kendall_df$condition)

s.class(
  pcoa_coords,
  fac = grp,
  col = c("#00A087","#E64B35"),
  cellipse = 1.5,
  cstar = 1,
  axesell = FALSE
)
```

---

### PCoA with statistical ellipses

```r
ggplot(
  pcoa_kendall_df,
  aes(PCoA1, PCoA2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) +
  theme_minimal()
```

---

## Beta Diversity Through Group Centroids

### Load required package

```r
library(dplyr)
```

---

### Compute centroids

```r
centroids <- pcoa_kendall_df %>%
  group_by(condition) %>%
  summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2))
```

---

### Merge centroid coordinates with original dataframe

```r
df2 <- merge(
  pcoa_kendall_df,
  centroids,
  by = "condition",
  suffixes = c("", "_cent"))
```

---

### Plot centroid-based beta diversity visualization

```r
ggplot(
  df2,
  aes(PCoA1, PCoA2, color = condition)) +
  geom_point(size = 3) +
  geom_segment(aes(xend = PCoA1_cent, yend = PCoA2_cent), alpha = 0.9) +
  theme_minimal()
```

---

## 🔹 STEP 6: Differential Abundance Statistics (Species associated with IBD Disease patients and Control Patients)

### Purpose

Identify species associated with disease and control groups using Wilcoxon statistical testing.

---

### Separate the control and disease metadata

```r
control_metadata <- metadata_df[metadata_df$study_condition == "Control",]
disease_metadata <- metadata_df[metadata_df$study_condition == "UC",]
```

---

### Subset species profiles

```r
control_species_profile <- species_matrix_norm[rownames(control_metadata),]
disease_species_profile <- species_matrix_norm[rownames(disease_metadata),]
```

---

### Use wilcox.test for one species

### Coprococcus_catus is Health Associated Core Species

```r
wilcox.test(disease_species_profile$Coprococcus_catus, control_species_profile$Coprococcus_catus)
```

---

### Mean abundance of Coprococcus_catus in disease samples

```r
mean(disease_species_profile$Coprococcus_catus)
```

---

### Mean abundance of Coprococcus_catus in control samples

```r
mean(control_species_profile$Coprococcus_catus)
```

---

### Difference of mean abundance

#### Disease minus Control

#### If -ve then Control associated

#### If +ve then Disease associated

```r
mean(disease_species_profile$Coprococcus_catus) - mean(control_species_profile$Coprococcus_catus)
```

---

## 🔹 STEP 7: Use wilcox.test for all the species using wilcox batch function

### Purpose

Perform batch differential abundance testing for all microbial species.

---

### Install required package

```r
install.packages("dplyr")
```

---

### Load required package

```r
library(dplyr)
```

---

### Wilcoxon batch function

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
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value
    type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0))
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

  out <- out %>% mutate(directionality = ifelse(type_array == 0, 0,
                               ifelse(p_adjust <= 0.05, 3 * sign(type_array),
                               ifelse(p_array <= 0.05, 2 * sign(type_array),
                               type_array))))
  return(out)
}
```

---

### Run the function for our data

```r
wilcox_result <- wilcox_batch(t(disease_species_profile), t(control_species_profile))
```

---

## 🔹 STEP 8: Volcano Plot Visualization

### Purpose

Visualize significantly enriched and depleted microbial species using volcano plots.

---

### Install required packages

```r
install.packages("ggplot2")
install.packages("ggrepel")
```

---

### Load required packages

```r
library(ggplot2)
library(ggrepel)
```

---

### Add the log2 fold change

### value 1 means disease abundance is double than control

```r
wilcox_result$log2_fold_change <- log2(wilcox_result$mean1_array / wilcox_result$mean2_array)
```

---

### Add significance column based on q-value

```r
wilcox_result$significance <- ifelse(wilcox_result$p_adjust <= 0.05, "significant", "non-significant")
```

---

### Use negative log10 q-value for visualization

### Higher value means more significant

```r
wilcox_result$logQValue <- -log10(wilcox_result$p_adjust)
```

---

### Add color code based on log2 fold change

```r
wilcox_result$color_code <- ifelse(wilcox_result$log2_fold_change > 1, "Enriched",
                                    ifelse(wilcox_result$log2_fold_change < -1, "Depleted", "Moderate Change"))
```

---

### Filter dataframe where only significant species are present

```r
wilcox_result_filt <- wilcox_result[which(wilcox_result$significance == "significant"),]
```

---

### Further filter species having q-value <= 0.001

```r
wilcox_result_filt2 <- wilcox_result_filt[which(wilcox_result_filt$p_adjust <= 0.001),]
```

---

### Plot Volcano Plot

```r
pdf("Wilcox_results.pdf", height = 8, width = 15)

ggplot(wilcox_result_filt2,
  aes(x = log2_fold_change, y = logQValue)) +
  geom_point(aes(color = color_code), size = 3) +
  geom_text_repel(data = wilcox_result_filt2, aes(label = rownames(wilcox_result_filt2)), size = 4, max.overlaps = 45) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey", linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(
    values = c("Enriched" = "royalblue3", "Depleted" = "palevioletred2", "Moderate Change" = "gray3")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (IBD/Control)", y = "-logQValue") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12))

dev.off()
```

---

## 💡 Key Notes

* Normalize abundance profiles before diversity analysis
* Use Wilcoxon test for non-parametric microbiome comparisons
* Use FDR-adjusted p-values (q-values) for multiple testing correction
* Bray-Curtis and Kendall distance are commonly used in microbiome beta diversity
* PCoA helps visualize microbial community separation
* Volcano plots help identify biologically important taxa

````
