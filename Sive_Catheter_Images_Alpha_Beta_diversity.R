#Install and call libraries
#install.packages("vegan")
#install.packages("phyloseq")
#install.packages("ggplot2")
#install.packages("dplyr")

library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(microbiome)

#import relative abundance table
abund_data <- read.csv("C:/Users/omx419/Desktop/Siva_Catheter/Siva_Catheter_Figures/Siva_Catheter_Figures/05192023_Redo_Analysis_Figures/05192023_Catheter_data_for_R.csv")

#import metadata table
metadata <- read.csv("C:/Users/omx419/Desktop/Siva_Catheter/Siva_Catheter_Figures/Siva_Catheter_Figures/05022023_Catheter_List_and_Category.csv")

#set row names
rownames(abund_data) <- abund_data$X
abund_data$X <- NULL


#Alpha Diversity
# Calculate alpha diversity for each treatment
alpha_diversity <- list()
for (Category in unique(metadata$Category)) {
  subset_abund_data <- abund_data[, metadata$Category == Category]
  alpha_diversity[[Category]] <- diversity(subset_abund_data, index="shannon")
}

# Combine alpha diversity values into a data frame for plotting
alpha_df <- data.frame(Category=names(alpha_diversity), 
                       alpha=unlist(alpha_diversity))

# Create boxplot of alpha diversity by treatment
boxplot(alpha ~ Category, data=alpha_df, 
        xlab="Category", ylab="Alpha diversity", main="Alpha diversity by Category")
        
        
#Beta Diversity
# Calculate Bray-Curtis dissimilarity matrix for all samples
bray_curtis <- vegdist(abund_data, method="bray")

# Perform PCoA on the Bray-Curtis dissimilarity matrix
pcoa <- cmdscale(bray_curtis, k=2, eig=TRUE)

# Combine PCoA coordinates and metadata into a data frame for plotting
pcoa_df <- data.frame(sample=rownames(abund_data), 
                      PC1=pcoa$points[,1], 
                      PC2=pcoa$points[,2], 
                      Category=metadata$Category)
rownames(pcoa_df) <- pcoa$sample
pcoa_df$sample <- NULL                      

write.csv(pcoa_df, file = "C:/Users/omx419/Desktop/Siva_Catheter/Siva_Catheter_Figures/Siva_Catheter_Figures/05192023_Siva_beta_diversity_pcoa_table.csv")

# Create PCoA plot with different symbols for each treatment
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Category)) +
  geom_point() +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2") +
  theme_bw()


#Evenness and Richness

rownames(metadata) <- metadata$Catheter
metadata$CAtheter <- NULL

comm_matrix <- otu_table(abund_data, taxa_are_rows=TRUE)
category_comm_matrix <- comm_matrix[, metadata$Category %in% c("Asymptomatic", "Symptomatic")]

richness <- rowSums(abund_data > 0)
abundance_data_T = t(abund_data)
evenness <- evenness(abundance_data_T, index = "all", zeroes = TRUE, detection = 0)


alpha_diversity <- data.frame(Richness=richness, Evenness=evenness$evar, Category=metadata$Category)

write.csv(alpha_diversity, file = "C:/Users/omx419/Desktop/Siva_Catheter/Siva_Catheter_Figures/Siva_Catheter_Figures/05192023_Siva_alpha_diversity.csv")

par(mfrow=c(1,2))
boxplot(Richness ~ Category, data=alpha_diversity, main="Richness")
boxplot(Evenness ~ Category, data=alpha_diversity, main="Evenness")