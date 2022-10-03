#Import data
file <- read.csv("OTU_limit.tsv", sep = "\t")

df <- as.data.frame(file)
View(df)

#Modify the name of file column

library(stringr)

df$file <- str_extract(df$file, "[:digit:]{7}")

#We are going to have 4 dataframes, one per each sample

df1 <- df[df$file == 8490409, c("taxon_name", "reads")]
df2 <- df[df$file == 8490411, c("taxon_name", "reads")]
df3 <- df[df$file == 8490412, c("taxon_name", "reads")]
df4 <- df[df$file == 8490413, c("taxon_name", "reads")]

library(ggplot2)
library(cowplot)

p1 <- ggplot(df1, aes(reorder(taxon_name, reads), reads)) +
  geom_col(aes(fill = reads)) +
  scale_fill_gradient(low = "blue",  high = "Purple") +
  ggtitle("Sample 8490409") +
  xlab("Taxon") + ylab("Reads") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip()

p2 <- ggplot(df2, aes(reorder(taxon_name, reads), reads)) +
  geom_col(aes(fill = reads)) +
  scale_fill_gradient(low = "blue",  high = "Purple") +
  ggtitle("Sample 8490411") +
  xlab("Taxon") + ylab("Reads") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip()

p3 <- ggplot(df3, aes(reorder(taxon_name, reads), reads)) +
  geom_col(aes(fill = reads)) +
  scale_fill_gradient(low = "blue",  high = "Purple") +
  ggtitle("Sample 8490412") +
  xlab("Taxon") + ylab("Reads") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip()

p4 <- ggplot(df4, aes(reorder(taxon_name, reads), reads)) +
  geom_col(aes(fill = reads)) +
  scale_fill_gradient(low = "blue",  high = "Purple") +
  ggtitle("Sample 8490413") +
  xlab("Taxon") + ylab("Reads") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip()

p1 <- as_grob(p1)
p2 <- as_grob(p2)
p3 <- as_grob(p3)
p4 <- as_grob(p4)

cowplot::plot_grid(p1, p2, p3, p4, nrow = 2)


###########
#DIVERSITY#
###########

file <- read.csv("OTU.tsv", sep = "\t")

df <- as.data.frame(file)
View(df)

#Modify the name of file column

library(stringr)

df$file <- str_extract(df$file, "[:digit:]{7}")

#We are going to have 4 dataframes, one per each sample

df1 <- df[df$file == 8490409, c("file", "taxon_id", "reads")]
df2 <- df[df$file == 8490411, c("file", "taxon_id", "reads")]
df3 <- df[df$file == 8490412, c("file", "taxon_id", "reads")]
df4 <- df[df$file == 8490413, c("file", "taxon_id", "reads")]

#taxon_id as columns, reads as values

library(tidyverse)

dfv1 <- df1 %>% group_by(file) %>% pivot_wider(names_from = taxon_id, values_from = reads)
dfv2 <- df2 %>% group_by(file) %>% pivot_wider(names_from = taxon_id, values_from = reads)
dfv3 <- df3 %>% group_by(file) %>% pivot_wider(names_from = taxon_id, values_from = reads)
dfv4 <- df4 %>% group_by(file) %>% pivot_wider(names_from = taxon_id, values_from = reads)


dfv <- dplyr::bind_rows(dfv1, dfv2, dfv3, dfv4) #Join the tables in one
dfv[is.na(dfv)] <- 0 #NA to 0
names <- dfv$file #to keep the row names for later
dfv <- dfv[,-1] #remove file column
dfv <- as.matrix(dfv) #convert to matrix
row.names(dfv) <- names #file as rowname

#Richness diversity

apply(dfv[]>0,1,sum) #total in each sample


library(vegan)

#Shannon alpha-diversity

exp(diversity(dfv, index = "shannon"))

#Bray-Curtis beta-diversity

vegdist(dfv, method = "bray")
