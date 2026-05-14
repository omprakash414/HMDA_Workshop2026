load("C:\\Users\\abc\\Downloads\\PEN_UC_data.RData")

library(vegan)
library(dplyr)
library(ggplot2)
library(pcaPP)
library(ade4)
library(stats)
library(pheatmap)
library(ggrepel)

rownames(Species_Profile_PEN_UC_norm)[order(rownames(Species_Profile_PEN_UC_norm))]

# remove S35 row, as it has no follow-up
Species_Profile_PEN_UC_norm <- Species_Profile_PEN_UC_norm[rownames(Species_Profile_PEN_UC_norm) != "S35", ]

Species_Profile_PEN_UC_norm_w4 <- Species_Profile_PEN_UC_norm[grep("_I$", rownames(Species_Profile_PEN_UC_norm)), ]
Species_Profile_PEN_UC_norm_baseline <- Species_Profile_PEN_UC_norm[!(rownames(Species_Profile_PEN_UC_norm)%in%rownames(Species_Profile_PEN_UC_norm_w4)),]

rownames(Species_Profile_PEN_UC_norm_w4) <- gsub("_I","", rownames(Species_Profile_PEN_UC_norm_w4))

Metadata_PEN_UC_baseline <- subset(Metadata_PEN_UC, select = c(study_code,age,gender,clinical_remission_w4,disease_score_baseline,albumin_baseline,CRP_baseline,FCP_baseline,Hb_baseline,diff_disease_score))
Metadata_PEN_UC_w4 <- subset(Metadata_PEN_UC, select = c(study_code,age,gender,clinical_remission_w4,disease_score_w4,albumin_w4,CRP_w4,FCP_w4,Hb_w4,diff_disease_score))

Species_Profile_PEN_UC_norm_baseline <- Species_Profile_PEN_UC_norm_baseline[rownames(Metadata_PEN_UC_baseline),]
Species_Profile_PEN_UC_norm_w4 <- Species_Profile_PEN_UC_norm_w4[rownames(Metadata_PEN_UC_w4),]


### Computing Alpha diversity

## Baseline
shannon_index_baseline <- diversity(Species_Profile_PEN_UC_norm_baseline, index = "shannon")

shannon_df_baseline <- data.frame(
  Sample = rownames(Species_Profile_PEN_UC_norm_baseline),
  Shannon_Index = shannon_index_baseline,
  clinical_remission_w4 = Metadata_PEN_UC_baseline$clinical_remission_w4)

shannon_df_baseline$clinical_remission_w4 <- ifelse(shannon_df_baseline$clinical_remission_w4 == 1,"Remission","No Remission")

# Week 4
shannon_index_w4 <- diversity(Species_Profile_PEN_UC_norm_w4, index = "shannon")
shannon_df_w4 <- data.frame(
  Sample = rownames(Species_Profile_PEN_UC_norm_w4),
  Shannon_Index = shannon_index_w4,
  clinical_remission_w4 = Metadata_PEN_UC_w4$clinical_remission_w4)

shannon_df_w4$clinical_remission_w4 <- ifelse(shannon_df_w4$clinical_remission_w4 == 1,"Remission","No Remission")

boxplot(shannon_df_baseline$Shannon_Index,shannon_df_w4$Shannon_Index, names = c("Baseline","Week4"), col = c("royalblue","gold"), outline = F)
wilcox.test(shannon_df_baseline$Shannon_Index,shannon_df_w4$Shannon_Index, paired = T)


###### Beta Diversity ######

### Between Baseline and Week 4

mat <- as.matrix(Species_Profile_PEN_UC_norm)
dist_kendall <- as.dist(1 - cor.fk(t(mat)) / 2)

pcoa_kendall <- dudi.pco(dist_kendall, scannf = FALSE, nf = 20)
pcoa_kendall_df <- as.data.frame(pcoa_kendall$li[,1:2])
colnames(pcoa_kendall_df) <- c("PCoA1", "PCoA2")
pcoa_kendall_df$condition <- ifelse(grepl("_I$", rownames(Species_Profile_PEN_UC_norm)),"week4","baseline")

grp <- as.factor(pcoa_kendall_df$condition)

s.class(
  pcoa_kendall_df,
  fac = grp,
  col = c("royalblue1","orange"),
  cellipse = 1.5,
  cstar = 1,
  axesell = FALSE,
  cpoint  = 2)

adonis_kendall <- adonis2(dist_kendall ~ condition, data = pcoa_kendall_df)
adonis_kendall

##### Baseline Responders vs Non-responders 

mat2 <- as.matrix(Species_Profile_PEN_UC_norm_baseline)
dist_kendall_bs <- as.dist(1 - cor.fk(t(mat2)) / 2)

pcoa_kendall_bs <- dudi.pco(dist_kendall_bs, scannf = FALSE, nf = 20)
pcoa_kendall_df_bs <- as.data.frame(pcoa_kendall_bs$li[,1:2])
colnames(pcoa_kendall_df_bs) <- c("PCoA1", "PCoA2")
pcoa_kendall_df_bs$condition <- Metadata_PEN_UC_baseline$clinical_remission_w4
pcoa_kendall_df_bs$condition <- ifelse(pcoa_kendall_df_bs$condition == 1,"responder","non-responder")

grp <- as.factor(pcoa_kendall_df_bs$condition)

s.class(
  pcoa_kendall_df_bs,
  fac = grp,
  col = c("royalblue1","orange"),
  cellipse = 1.5,
  cstar = 1,
  axesell = FALSE,
  cpoint  = 2)

adonis_kendall_bs <- adonis2(dist_kendall_bs ~ condition, data = pcoa_kendall_df_bs)
adonis_kendall_bs


##### Week4 Responders vs Non-responders 

mat3 <- as.matrix(Species_Profile_PEN_UC_norm_w4)
dist_kendall_w4 <- as.dist(1 - cor.fk(t(mat3)) / 2)

pcoa_kendall_w4 <- dudi.pco(dist_kendall_w4, scannf = FALSE, nf = 20)
pcoa_kendall_df_w4 <- as.data.frame(pcoa_kendall_w4$li[,1:2])
colnames(pcoa_kendall_df_w4) <- c("PCoA1", "PCoA2")
pcoa_kendall_df_w4$condition <- Metadata_PEN_UC_w4$clinical_remission_w4
pcoa_kendall_df_w4$condition <- ifelse(pcoa_kendall_df_w4$condition == 1,"responder","non-responder")

grp <- as.factor(pcoa_kendall_df_w4$condition)

s.class(
  pcoa_kendall_df_w4,
  fac = grp,
  col = c("royalblue1","orange"),
  cellipse = 1.5,
  cstar = 1,
  axesell = FALSE,
  cpoint  = 2)

adonis_kendall_w4 <- adonis2(dist_kendall_w4 ~ condition, data = pcoa_kendall_df_w4)
adonis_kendall_w4


##################################################################
## Mean Abundance wise top Species ##

species_mean_abundance <- function(abundance_matrix_norm){
  abundance_mean <- as.data.frame(apply(abundance_matrix_norm,2,mean))
  abundance_mean$species <- rownames(abundance_mean)
  rownames(abundance_mean) <- NULL
  colnames(abundance_mean)[1] <- "mean_abundance"
  abundance_mean <- abundance_mean[,c(2,1)]
  abundance_mean <- abundance_mean[order(abundance_mean$mean_abundance, decreasing = T),]
  return(abundance_mean)
}

baseline_species_mean_abundance <- species_mean_abundance(Species_Profile_PEN_UC_norm_baseline)
w4_species_mean_abundance <- species_mean_abundance(Species_Profile_PEN_UC_norm_w4)

#pdf("G:\\.shortcut-targets-by-id\\1cCOYITKBwsY8RsVn82_B6lk4CikVsrwU\\PEN_UC_Project\\Plots\\abundance_prevalence_plots\\baseline_abundance_top20.pdf",height = 4,width = 5)
pheatmap(Species_Profile_PEN_UC_norm_baseline[,baseline_species_mean_abundance$species[1:20]],
         cluster_rows = F,cluster_cols = F, fontsize_row = 8, fontsize_col = 10, labels_col = gsub("_"," ", baseline_species_mean_abundance$species[1:20]))
#dev.off()

#pdf("G:\\.shortcut-targets-by-id\\1cCOYITKBwsY8RsVn82_B6lk4CikVsrwU\\PEN_UC_Project\\Plots\\abundance_prevalence_plots\\w4_abundance_top20.pdf",height = 4,width = 5)
pheatmap(Species_Profile_PEN_UC_norm_w4[,w4_species_mean_abundance$species[1:20]],
         cluster_rows = F,cluster_cols = F, fontsize_row = 8, fontsize_col = 10, labels_col = gsub("_"," ", w4_species_mean_abundance$species[1:20]))
#dev.off()

#species <- c("Segatella_copri","Enterococcus_avium","Intestinibacter_bartlettii","Enterococcus_gallinarum")

taxa_prof_bs <- Species_Profile_PEN_UC_norm_baseline
taxa_prof_bs$response <- Metadata_PEN_UC_baseline$clinical_remission_w4
taxa_prof_bs$response <- ifelse(taxa_prof_bs$response == 1,"responder","non_responder")

boxplot(taxa_prof_bs[taxa_prof_bs$response == 'responder','Intestinibacter_bartlettii'],taxa_prof_bs[taxa_prof_bs$response == 'non_responder','Intestinibacter_bartlettii'])
boxplot(taxa_prof_bs[taxa_prof_bs$response == 'responder','Segatella_copri'],taxa_prof_bs[taxa_prof_bs$response == 'non_responder','Segatella_copri'])

wilcox.test(taxa_prof_bs[taxa_prof_bs$response == 'responder','Intestinibacter_bartlettii'],taxa_prof_bs[taxa_prof_bs$response == 'non_responder','Intestinibacter_bartlettii'])

##################################################################

boxplot(Metadata_PEN_UC_baseline$disease_score_baseline, Metadata_PEN_UC_w4$disease_score_w4, names = c("Baseline","Follow-up"), col = c("royalblue1","coral"), ylab = "Disease Score")
wilcox.test(Metadata_PEN_UC_baseline$disease_score_baseline, Metadata_PEN_UC_w4$disease_score_w4, paired = T)

wilcox_batch_paired_mean = function(x,y)
{
  p_array <- NULL;
  type_array <- NULL;
  mean1_array <- NULL;
  mean2_array <- NULL;
  x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
  y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
  z <- intersect(rownames(x),rownames(y));
  for(i in 1:length(z))
  {
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]),paired=TRUE)$p.value;
    type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
    mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
    i <- i + 1;
  }
  out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array, method = 'BH'),mean1_array,mean2_array));
  rownames(out) <- z;
  colnames(out)[3] <- "p_adj"
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  return(t(out));
}

pairwise_wilcox_mean_bs_w4 <- as.data.frame(wilcox_batch_paired_mean(t(Species_Profile_PEN_UC_norm_baseline),t(Species_Profile_PEN_UC_norm_w4)))

### volcano plot ###

pairwise_wilcox_mean_bs_w4$transformed_p <- -log10(pairwise_wilcox_mean_bs_w4$p_array)
pairwise_wilcox_mean_bs_w4$log2_fold_change <- log2(pairwise_wilcox_mean_bs_w4$mean2_array / pairwise_wilcox_mean_bs_w4$mean1_array)

pairwise_wilcox_mean_bs_w4$color_code <- ifelse(pairwise_wilcox_mean_bs_w4$log2_fold_change > 1, "color1",ifelse(pairwise_wilcox_mean_bs_w4$log2_fold_change < -1, "color2","color3"))
pairwise_wilcox_mean_bs_w4_trimmed <- pairwise_wilcox_mean_bs_w4[pairwise_wilcox_mean_bs_w4$transformed_p >1,]

View(pairwise_wilcox_mean_bs_w4_trimmed)

table(sign(pairwise_wilcox_mean_bs_w4_trimmed$type_array))


ggplot(pairwise_wilcox_mean_bs_w4_trimmed, aes(x = log2_fold_change, y = transformed_p)) +
  geom_point(aes(color = color_code), size = 3) +  # Increase size of the dots
  geom_text_repel(data = pairwise_wilcox_mean_bs_w4_trimmed,
                  aes(label = rownames(pairwise_wilcox_mean_bs_w4_trimmed)), size = 4, max.overlaps = 45) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey", size = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  scale_color_manual(values = c("color1" = "royalblue3", "color2" = "deeppink1", "color3" = "gray3")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (Week4/Baseline)", y = "-Log10 (p-value)") +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))


##################################################################

cor.test(Species_Profile_PEN_UC_norm$Eubacterium_rectale, c(Metadata_PEN_UC_baseline$disease_score_baseline,Metadata_PEN_UC_w4$disease_score_w4), method = 'spearman')


disease_score <- c(Metadata_PEN_UC_baseline$disease_score_baseline,Metadata_PEN_UC_w4$disease_score_w4)


corr_results <- data.frame(
  rho = numeric(ncol(Species_Profile_PEN_UC_norm)),
  pval = numeric(ncol(Species_Profile_PEN_UC_norm)),
  row.names = colnames(Species_Profile_PEN_UC_norm))


for(i in seq_len(ncol(Species_Profile_PEN_UC_norm))) {
  
  species_values <- Species_Profile_PEN_UC_norm[, i]
  
  # Remove NA pairs
  valid <- complete.cases(species_values, disease_score)
  
  # Skip if too few observations
  if(sum(valid) < 3) next
  
  # Skip constant vectors
  if(length(unique(species_values[valid])) < 2) next
  
  # Spearman correlation
  test <- cor.test(
    species_values[valid],
    disease_score[valid],
    method = "spearman"
  )
  
  # Store results
  corr_results$rho[i] <- test$estimate
  corr_results$pval[i] <- test$p.value
}

corr_results <- na.omit(corr_results)

View(corr_results)

taxa_prof <- Species_Profile_PEN_UC_norm
taxa_prof$disease_score <- c(Metadata_PEN_UC_baseline$disease_score_baseline,Metadata_PEN_UC_w4$disease_score_w4)




