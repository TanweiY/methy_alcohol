library(readxl)
scores <- read_excel("alchohol_methy/scores.xlsx", 
                     sheet = "alcohol_socre_details")

summary(as.factor(scores$Study))
#Elliott2014      EpiTob   McCartney         SSt   Zhang2016 
#187           5         233         121           4 

library(VennDiagram)

# Assuming your dataframe is named df
# Convert the dataframe into a list of CpG vectors by Methylation_Score
cpg_list <- split(scores$CpG, scores$Study)

# find the one cpgs in all scores
Reduce(intersect, cpg_list) # only one "cg06690548" was repeated reported

## find the cpgs that were in  3 groups
all_cpgs <- unlist(cpg_list)
cpg_counts <- table(all_cpgs)
names(cpg_counts[cpg_counts == 3])

# Hexadecimal color specification 
library(RColorBrewer)
brewer.pal(n = 6, name = "Blues")

if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library("ggplot2")
ggVennDiagram(cpg_list, label_alpha = 0, set_size = 4, label = c("count"),
              order.intersect.by= 'size', order.set.by = 'size',
              edge_size = 1)+
              scale_fill_gradient(low = "#C6DBEF", high = "#3182BD") 

# 6*4 






