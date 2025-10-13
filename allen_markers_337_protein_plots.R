library(dplyr)
library(pheatmap)
library(gridExtra)

d <- read.csv("Supplementary S5_ALlen gene whehter in proteomics.csv",
              header = T, stringsAsFactors = F)
head(d)
d <- filter(d, Detection_in_Proteomics == TRUE)

my_gene_list <- d$Gene_Symbol
length(my_gene_list)

# there are only 249 unique genes 
length(unique(my_gene_list))

my_gene_list <- unique(my_gene_list)


table(d$Layer)

# generate protein spatial maps for these 337 proteins
dim(d1)
colnames(d1)

# use at least 15%, 2687 proteins
337/20

# make sure do not use dataset without the empty voxles.
# read in the 2D dataset
# d1 <- read.csv("MOP_with_clustering.csv",
#                row.names = 1, header = T, stringsAsFactors = F)
# 
# colnames(d1)
# 
# d1$PG.Genes <- gsub(";.*", "", d1$PG.Genes)
# 
# table(d1$PG.Genes %in% my_gene_list)
# 
# d1 <- d1[d1$PG.Genes %in% my_gene_list, ]
# dim(d1)
# colnames(d1)




d1 <- read.csv("MOp_norm_for_Supplemental_Data.csv",
               row.names = 1, header = T, stringsAsFactors = F)

colnames(d1)

d1$PG.Genes <- gsub(";.*", "", d1$PG.Genes)

table(d1$PG.Genes %in% my_gene_list)

d1 <- d1[d1$PG.Genes %in% my_gene_list, ]
dim(d1)
colnames(d1)

row.names(d1) <- d1$PG.Genes
d1 <- select(d1, R1:R493)
d1 <- d1[match( my_gene_list, row.names(d1)), ]
row.names(d1)

summary(apply(d1, 1, function(x) sum(!is.na(x))))

# remove rows that has less than 20% 
d1 <- d1[apply(d1, 1, function(x) sum(!is.na(x))) > 493*0.2,]

# finished the data prep

dim(d1)[1]/20
spatial_protein_plot_log2_scale("Prdx3")


for(k in 1:1) {
  png(paste0("MoP_allen_marker_at_leastp20_set_", k, ".png" ),
      width = 4500, height = 3000, res = 300)
  
  
  start <- k*20 - 19
  end <- k*20
  
  
  heatmap_list <- vector("list", 20)
  
  for(i in start:end) {
    plot_tmp <- spatial_protein_plot_log2_scale(row.names(d1)[i])$gtable
    
    # turn start number to 1, 2, 3, to 20
    heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
    if( heat_map_index == 0) {heat_map_index <- 20}
    
    heatmap_list[[heat_map_index]] <- plot_tmp
    
  }
  
  grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)
  
  dev.off()
}


k = 1
################# the remaining 7 ###################


png(paste0("MoP_All_protein_protein_set_", 135, ".png" ),
    width = 4500, height = 3000, res = 300)


start <- 2681
end <- 2687


heatmap_list <- vector("list", 7)

for(i in start:end) {
  plot_tmp <- spatial_protein_plot_log2_scale(row.names(d1)[i])$gtable
  
  # turn start number to 1, 2, 3, to 20
  heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
  if( heat_map_index == 0) {heat_map_index <- 20}
  
  heatmap_list[[heat_map_index]] <- plot_tmp
  
}

grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)

dev.off()



2687 - 134*20

heatmap_list <- vector("list", 7)

for(i in 2681:2687) {
  plot_tmp <- spatial_protein_plot_log2_scale(row.names(d1)[i])$gtable
  
  # turn start number to 1, 2, 3, to 20
  heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
  if( heat_map_index == 0) {heat_map_index <- 20}
  
  heatmap_list[[heat_map_index]] <- plot_tmp
  
}

grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)

dev.off()



# chekc the function with a gene
spatial_protein_plot_log2_scale("Prdx3")
my_gene <- "Prdx3"

spatial_protein_plot_log2_scale <- function(my_gene) {
  
  
  dt <- d1[my_gene,]
  
  
  # dt <- 2^dt
  # max <- max(dt, na.rm = T)
  # dt <- 100* dt / max
  
  exp_data <- as.numeric(dt)
  
  exp_data <- exp_data - median(exp_data, na.rm = T)
  
  # put the data into a matrix
  pm <- matrix(exp_data, nrow = nrow, ncol = ncol)
  
  
  data <- pm
  
  
  # Generate heatmap using pheatmap with handling of NA values
  p <-   pheatmap(data, scale = "none", 
                  # breaks = breaksList,
                  # color = my_palette_with_na, 
                  cluster_rows = F,
                  cluster_cols = F,
                  #display_numbers = mat2,
                  legend = T,
                  main = gsub("_MOUSE", "", my_gene))
  return(p)
  
}



plot(1:10, 1:10)
