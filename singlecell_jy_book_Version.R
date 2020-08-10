library(ggplot2)
library(cowplot)
library(dplyr)
library(monocle)
install.packages("Rmisc",dependencies=TRUE)
install.packages("http://s3-us-west-2.amazonaws.com/10x.files/code/cellrangerRkit-2.0.0.tar.gz",repos=NULL)
library(cellrangerRkit)
library(HSMMSingleCell)

my_dir="/home/node01/Desktop/Chapter_V. RNAseq 전사체 데이터 분석"
gbm <-load_cellranger_matrix(my_dir) # load data from 10x cellranger pipeline 
class(gbm)
dim(exprs(gbm))

my_feat <-fData(gbm)
names(my_feat) <-c('id','gene_short_name')
my_cds <-newCellDataSet(exprs(gbm),
                        phenoData=new("AnnotatedDataFrame",data=pData(gbm)),
                        featureData =new("AnnotatedDataFrame",data=my_feat),
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
my_cds
my_cds<-estimateSizeFactors(my_cds)
my_cds<-estimateDispersions(my_cds)

my_cds <- detectGenes(my_cds,min_expr = 0.1)
head(fData(my_cds))

summary(fData(my_cds)$num_cells_expressed)
sum((exprs(my_cds["ENSG00000239945",])))
sum((exprs(my_cds["ENSG00000238009",])))
summary(pData(my_cds)$num_genes_expressed)
x<-pData(my_cds)$num_genes_expressed
x_1 <-(x-mean(x))/sd(x)
summary(x_1)
df <-data.frame(x=x_1)
ggplot(df,aes(x))+geom_histogram(bins=50)+geom_vline(xintercept=c(-2,2),linetype="dotted",color="red")

#umi
pData(my_cds)$UMI <-Matrix::colSums(exprs(my_cds)) #count data
ggplot(pData(my_cds),aes(num_genes_expressed,UMI))+geom_point()

disp_table <-dispersionTable(my_cds)
head(disp_table)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <-subset(disp_table,mean_expression>=0.1)
my_cds <-setOrderingFilter(my_cds,unsup_clustering_genes$gene_id)
plot_ordering_genes(my_cds)

#cluster
my_cds<-reduceDimension(my_cds,max_components = 2,num_dim=10,reduction_method = 'tSNE',verbose = TRUE)
#whether verbose is true or not isnt important
?reduceDimension
my_cds <-clusterCells(my_cds,num_clusters=15)
head(pData(my_cds))
my_cluster_dim_10 <- pData(my_cds)$Cluster
plot_cell_clusters(my_cds)

#
my_vector <-rep('no',nrow(pData(my_cds)))
my_vector
my_vector[pData(my_cds)$Cluster==1] <-rep('yes',sum(pData(my_cds)$Cluster==1))
my_vector[F]
pData(my_cds)$test <-my_vector
head(pData(my_cds))
length(unsup_clustering_genes$gene_id)
de_cluster_one <- differentialGeneTest(my_cds[unsup_clustering_genes$gene_id,],fullModelFormulaStr = '~test',cores=4)
dim(de_cluster_one)
de_cluster_one %>% arrange(qval) %>% head()
plot_genes_jitter(my_cds['ENSG00000163131',],grouping="Cluster")
?plot_genes_jitter

#
pData(my_cds)$my_colour <-pData(my_cds)$Cluster ==4 |pData(my_cds)$Cluster ==5 |pData(my_cds)$Cluster ==13 
plot_cell_clusters(my_cds,color_by="my_colour")
pData(my_cds)$my_colour
#
expressed_genes <- row.names(subset(fData(my_cds),num_cells_expressed>=10))
expressed_genes

my_cds_subset <-my_cds[expressed_genes,pData(my_cds)$my_colour]
my_cds_subset <- detectGenes(my_cds_subset,min_expr = 0.1)
fData(my_cds_subset)$use_for_ordering <-fData(my_cds_subset)$num_cells_expressed >0.05
ncol(my_cds_subset)
table(fData(my_cds_subset)$use_for_ordering)
my_cds_subset <- reduceDimension(my_cds_subset,max_commponents=2,norm_method = 'log',num_dim=10,reuction_method="tSNE",verbose = T)

my_cds_subset <-clusterCells(my_cds_subset, verbose=F)






