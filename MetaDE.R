rm(list=ls()) #������ܿ� �ִ� environment�� �ִ� dataset�� �� ���������� ���� �� 
setwd('C:/Users/82104/Documents/prostate')
#--------------------importing data into R------------------------------------------------------#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

#metade�� ����� impute�� combinate�� ��� �־���Ѵ� metade�� cran ���� �ٿ� �޾Ƽ� ���� ��ġ�Ѵ�
library(MetaDE)

#meta qc�� ����� proto foreach iterators�� �ʿ��ϴ�
library(MetaQC)
study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace","Nanni","Dhanasekaran", "Tomlins")
#������ �� ����� �̸�/ 9���� ������� ���� �� �� ����� �������ڷ���� �������� �̸����� ǥ���� �� 
prostate.raw<-MetaDE.Read(study.names,skip=rep(1,9),via="txt",matched=T,log=F)
#MetaDE.Read�� �����͸� �����ͼ� meatde pakage���ʿ��� ���·� �ٲٴ� �� 
#ù��° ������ ���ϸ�! �׷��� ������ ���ϸ��� study.name���� ��� ǥ���Ѱ�!! �׷��� �� ���ϸ��� �´� �� �о���̴� ���� 
#via- txt/ csv �̰� ������ �����Ͱ� ,�� ���� �� ������ csv�� 
#skip���� 2�� survival data�� 1�̸� ������/ rep(1,9)->1�� 9�� �ݺ� ��
#matched �� genesymbol�̶� probeid�� ��ġ�ϳĴ� �Ű�
#log�� �����Ϳ� log2�� ����� �۾��� �ʿ��ϳĴ� �� 

#--------------------merge and filter data---------------------------------------------------------#
prostate.merged<-MetaDE.merge(prostate.raw)#microarray data set�� �ұ�Ģ�� ������ ��ġ����/ prostate.raw-list�������� �����ϱ� �̰� ��ġ�� �� ���� 
dim(prostate.merged[[1]][[1]])# [[]], 5���� �����͸� �ϳ��� �����ͷ� ��ġ�� [[1]]�̷��� �ϸ� ��ģ 5�� �߿� ù��°�ž�! dim�� ũ�⸦ ����� ��
prostate.filtered<-MetaDE.filter(prostate.merged,c(0.3,0.3))#ù������ list�� �ް� �״����� �����ڰ� ��������� ���͸� �ž� �ϴ��� ��Ÿ��
dim(prostate.filtered[[1]][[1]])#�׷��� ������ �پ�� 

#--------------------MetaQC------------------------------------------------------------------------#
library(MetaQC)
Data.QC<-list() #�� ����Ʈ�� ���� 

for(i in 1:9){
  colnames(prostate.filtered[[i]][[1]])<-prostate.filtered[[i]][[2]] #ù���� column�� normal�ε� �װ� ���ְ� ������ 
  Data.QC[[i]]<-impute.knn(prostate.filtered[[i]][[1]])$data
  #impute.knn�Լ�-���� ����� ���� ����� ����Ͽ� ���� ǥ���� �����͸� �ͼӽ�Ű�� �Լ�.

  print(dim(Data.QC[[1]]))
}
names(Data.QC)<-names(prostate.filtered)
start<-Sys.time()#��¥�� ����ð��� ��Ÿ�� 
ProstateQC<-MetaQC(Data.QC, "c2.all.v3.0.symbols.gmt", filterGenes=F, verbose=TRUE,isParallel=FALSE, resp.type="Twoclass")
#metaqc ���״�� qc�� �������ִ� ��
#ù��° ���-�� �����ʹ� ���� ǥ�� ������, ��� �Ǵ� Ŭ���� ��, ������ ID�� ��Ÿ���� x, y �� geneid(��Ʈ���� x�� �� �̸����� �����ڸ� �ٲ� �� ����) ��Ҹ� ������ �ϴ� ������� ǥ�õȴ�. 
#gene symbol lists�� ������ �ִ� ������ ��ġ��
#filtergenes-gene filtering�� ������ ������
#verbose-�α� ������ �μ����� ����
#isparallel - �����ϱ����� �ھ ������ ������
#ncores-������ true�� �׷��� core�� � ������
#resp.type- type of response variable 

runQC(ProstateQC, B=1e4, fileForCQCp="c2.all.v3.0.symbols.gmt")

#B-permutation testȽ��

QC_time<-Sys.time()-start
png(filename = "Prostate_QC0421.png", width = 3500, height = 3500,res=600)
plot(ProstateQC)
dev.off()
jpeg(filename = "Prostate_QC0421.jpeg", width = 3500, height = 3500,res=600)
plot(ProstateQC)
dev.off()

#���⼭ png�� �ϸ� �ڿ� ���� plot()�� png ���·� �����ϰڴٴ°Ű� dev.off�ϸ� ������� ���� plot�� ��� �ϳ��� png ���·� �����ϴ°ž�
#jpeg�� png�� ���� ����
#png plot plot dev.off �ϸ� plot�ΰ��� ���ÿ� png�� ��������°���
#�׷��� �߰��� plot()�� �ص� ���� �Ⱥ����ִµ� �� ������ �������� �׳� plot�� �����ϸ� �׸� ������ 

#-----------------------------------------------------------------------------------------------#
# (1) To remove the three studies ("Nanni","Dhanasekaran","Tomlins") with bad quality
# (2) To remerge the remaining six studies
# (3) To re-filter the data
#-----------------------------------------------------------------------------------------------#
study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace")

data.QC.raw<-list()
for(i in 1:length(study.names)){
  data.QC.raw[[i]]<-prostate.raw[[study.names[[i]]]]
}
names(data.QC.raw)#null ������
names(data.QC.raw)<-study.names
names(data.QC.raw)#�̸��� �� 

data.QC.merged<-MetaDE.merge(data.QC.raw)
dim(data.QC.merged[[1]][[1]])
data.QC.filtered<-MetaDE.filter(data.QC.merged,c(0.2,0.2))
#low expression gene���̶� test control���� variation�� ���� ��� �Ÿ�
dim(data.QC.filtered[[1]][[1]])
#---------------------- MetaDE-----------------------------------------------------------------#
start<-Sys.time()
MetaDE.Res<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Fisher"),nperm=300,asymptotic=F,ind.tail = "abs")
#��ġ���� ����� ������ metade.res�� 

MetaDE.fisherH<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Fisher"),nperm=300,asymptotic=F,ind.tail = "high")
MetaDE.fisherL<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Fisher"),nperm=300,asymptotic=F,ind.tail = "low")

MetaDE.stouffer<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Stouffer"),nperm=300,asymptotic=F)
MetaDE.fem<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("FEM"),nperm=300,asymptotic=F)
MetaDE.rem<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("REM"),nperm=300,asymptotic=F)




#modt-moderate t test/ two sample t test�� ����� �� �� �� �µ��� ������ �������� �����ϰ�?
#nperm�� permutation test �װŸ� 300�� �ϰڴٴ� �� 
#ind->independent��../ ind.method�� ����� ����� ���ϴ� ��
#meta.method�� meta analysis����� ���ϴ� �Ű� �������� �������� �ٶ�� �� �� 4������
#rth�� rop�� rop.oc ��� �� �� �ʿ��ѵ� 
#asymptotic - ��Ÿ �м����� p-���� ����ϱ� ���� �Ķ��Ʈ�� ����� �����ߴ��� ���θ� �����ϴ� ������ ��. �⺻���� FALSE�Դϴ�.


'���� ���⼭ MetaDE.Res1<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),meta.method=c("Fisher","sto~","fem","rem"),nperm=300,asymptotic=F)
 �̷��� ���ȴµ� �̷��ϱ� paired�� �����ش޶�� ������ ���� 
 �� �׷��ĸ� fem�̶� rem�� �� �����͵��� ����ǥ������ �ƴ����� �����ؾ��Ѵ� �׷��� �̰� �׷��� �ƴ����� paired=false�̷������� �˷������
 �� ������ �̰� 4���� �ɰ� ������ � meta������� ���� ������ �� ���ư��µ� �ƴ� �͵��� �־ ���� ����
 �ٵ� �̷��� �ص� ������ �� 
 �ȿ� ��Ű������ low ���� �����;��ϴµ� column���� ��������? �׷� �� ���� �� �ִ��� 
 �׷��� ��Ű�� �ٽ� ��ġ�� ��ġ�� 
 MetaDE_revised ��¼�� 
'

#min ������ 
'fisher�� ���� ������ ����ϱ� ������ �������� ���� �Ŀ� ����/2���� �� ���� ���� ã�Ƽ� �װ� �ι����ش�'

MetaDE.fisher <-MetaDE.fisherH

statH <-MetaDE.fisherH$meta.analysis$stat
statL <-MetaDE.fisherL$meta.analysis$stat
for(i in 1:length(statH)){
    MetaDE.fisher$meta.analysis$stat[i]<-2*(min(statH[i],statL[i]))
  
  }

pvalH <-MetaDE.fisherH$meta.analysis$pval
pvalL <-MetaDE.fisherL$meta.analysis$pval
for(i in 1:length(pvalH)){
  MetaDE.fisher$meta.analysis$pval[i]<-2*(min(pvalH[i],pvalL[i]))
  
}

FDR_H <-MetaDE.fisherH$meta.analysis$FDR
FDR_L <-MetaDE.fisherL$meta.analysis$FDR
for(i in 1:length(FDR_H)){
  MetaDE.fisher$meta.analysis$FDR[i]<-2*(min(FDR_H[i],FDR_L[i]))
  
}

AW_H <-MetaDE.fisherH$meta.analysis$AW.weight
AW_L <-MetaDE.fisherL$meta.analysis$AW.weight
for(i in 1:length(AW_H)){
  MetaDE.fisher$meta.analysis$AW.weight[i]<-2*(min(AW_H[i],AW_L[i]))
  
}

#�̸���ġ��
# category<-list()
category <-c(names(MetaDE.rem$meta.analysis),names(MetaDE.fisher$meta.analysis),names(MetaDE.fem$meta.analysis),names(MetaDE.stouffer$meta.analysis))
category<-unique(category)
category

#list�� �����Ϸ��� [[]]/ vector�� c�̰Ű� �Ѵ� for������ ���ٰ�����
#list�� index��ȣ�ε� ���� �����ϰ� �� �������ε� ���� ���� ��� [[��ȣ/�̸�]]�̷��� �ؾ���! 

MetaDE.Res$meta.analysis[['pval']]
#���ο� metaanalysis
for (i in category){
  MetaDE.Res$meta.analysis[[i]]<-cbind(MetaDE.fisher$meta.analysis[[i]],MetaDE.rem$meta.analysis[[i]],MetaDE.fem$meta.analysis[[i]],MetaDE.stouffer$meta.analysis[[i]])
  # colnames(MetaDE.Res$meta.analysis[[i]])<-c("Fisher","REM","FEM","Stouffer")[which(c(is.null(MetaDE.fisher$meta.analysis[category[i]]),is.null(MetaDE.Res.REM$meta.analysis[[i]]),is.null(MetaDE.Res.FEM$meta.analysis[[i]]))==F),is.null(MetaDE.stouffer$meta.analysis[category[i]])]
  colnames(MetaDE.Res$meta.analysis[[i]])<-c("Fisher","REM","FEM","Stouffer")[which(c(is.null(MetaDE.fisher$meta.analysis[[i]]),is.null(MetaDE.rem$meta.analysis[[i]]),is.null(MetaDE.fem$meta.analysis[[i]]),is.null(MetaDE.stouffer$meta.analysis[[i]]))==F)]
}

'which(c(1,2,3,4,5,1)==1) �̷��� �ϸ� �� ���Ϳ��� ������ �����ϴ� ��ġ�� ��ȯ�� �׷� 1�̶� 6�̰���
���� for�� ���� ���� c(f,r,f,s)[which(c(����))]�̷��� �����ݾ�
�̶� �� ����� meta.analysis�� pbval�� �ִ��� ������ �������� is null�Լ��� ������ f�ϱ� ä�����ִ� ��ġ�� ��ȯ��
1,2,3,4 / 1,3/ ��Ÿ ��� �׷� �� index�� �տ� �ִ� c(f,r,f,s)���� ��ġ�� �°� ��ȯ�ϴ°ž�
fisher�� ���� ��� ���� colname�� �ְ� rem/fem�� ���ŵ� �׷��� �׸��� �׷��� ������ �̸��� ���� ������ ������ �ž�
�׷��� �̰� �������ַ��� �� �̸��� �������ذ���! 

'

#--------------------------------------------------------------------------------------------

b<-Sys.time()-start
print(b)
mylty<-rep(c(1,2),c(6,4))
mycol<-c(rep("black",6),c("red","green","blue","orange"))
mylwd<-rep(c(1,2),c(6,4))
mypch<-1:10
png(filename = "MetaDE_Prostate0421.png", width = 3500, height = 3500,res=600)
draw.DEnumber(MetaDE.Res,0.05,mlty=mylty,mcol=mycol,mlwd=mylwd,mpch=mypch,FDR=T)
dev.off()
#-----------------------------------------------------------------------------------------------#
# (1) To remove the three studies ("Nanni","Dhanasekaran","Tomlins") with bad quality
# (2) To remerge the remaining six studies
# (3) To re-filter the data
#-----------------------------------------------------------------------------------------------#
study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace","Nanni","Dhanasekaran",
               "Tomlins")
data.QC.raw<-list()
for(i in 1:length(study.names)){
  data.QC.raw[[i]]<-prostate.raw[[study.names[[i]]]]
}
names(data.QC.raw)<-study.names
data.QC.merged<-MetaDE.merge(data.QC.raw)
dim(data.QC.merged[[1]][[1]])
data.QC.filtered<-MetaDE.filter(data.QC.merged,c(0.2,0.2))
dim(data.QC.filtered[[1]][[1]])
#---------------------- MetaDE-----------------------------------------------------------------#
start<-Sys.time()
#MetaDE.Res<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",9),meta.method=c("Fisher","REM","FEM","stouffer"),nperm=300,asymptotic=F)
#���⵵ ���� ó�� 4������ ������ �ٽ� �����ϰ� ��ġ�� �׷� ������ ���ľ��� 
b<-Sys.time()-startprint(b)
mylty<-rep(c(1,2),c(9,1))
mycol<-c(rep("black",6),c("red","pink","blue","yellow"))
mylwd<-rep(c(1,2),c(9,1))
mypch<-1:10
png(filename = "MetaDE_Prostate0421.png", width = 3500, height = 3500,res=600)
draw.DEnumber(MetaDE.Res,0.05,mlty=mylty,mcol=mycol,mlwd=mylwd,mpch=mypch,FDR=T)
dev.off()