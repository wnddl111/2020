rm(list=ls()) #우측상단에 있는 environment에 있는 dataset을 다 지워버리고 싶을 때 
setwd('C:/Users/82104/Documents/prostate')
#--------------------importing data into R------------------------------------------------------#
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

#metade를 깔려면 impute랑 combinate도 깔려 있어야한다 metade는 cran 에서 다운 받아서 직접 설치한다
library(MetaDE)

#meta qc를 깔려면 proto foreach iterators가 필요하다
library(MetaQC)
study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace","Nanni","Dhanasekaran", "Tomlins")
#논문을 쓴 사람들 이름/ 9명의 사람들이 논문 쓸 때 사용한 데이터자료들을 논문저자 이름으로 표시한 것 
prostate.raw<-MetaDE.Read(study.names,skip=rep(1,9),via="txt",matched=T,log=F)
#MetaDE.Read는 데이터를 가져와서 meatde pakage에필요한 형태로 바꾸는 것 
#첫번째 변수가 파일명! 그래서 위에서 파일명을 study.name으로 묶어서 표시한것!! 그래서 저 파일명에 맞는 걸 읽어들이는 거지 
#via- txt/ csv 이거 구분은 데이터가 ,로 구분 돼 있으면 csv임 
#skip에서 2면 survival data고 1이면 나머지/ rep(1,9)->1을 9번 반복 함
#matched 는 genesymbol이랑 probeid가 일치하냐는 거고
#log는 데이터에 log2를 씌우는 작업이 필요하냐는 것 

#--------------------merge and filter data---------------------------------------------------------#
prostate.merged<-MetaDE.merge(prostate.raw)#microarray data set을 불규칙한 순서로 합치세요/ prostate.raw-list형식으로 있으니까 이걸 합치는 것 같음 
dim(prostate.merged[[1]][[1]])# [[]], 5개의 데이터를 하나의 데이터로 합치면 [[1]]이렇게 하면 합친 5개 중에 첫번째거야! dim은 크기를 물어보는 것
prostate.filtered<-MetaDE.filter(prostate.merged,c(0.3,0.3))#첫번쨰는 list를 받고 그다음은 유전자가 어느정도로 필터링 돼야 하는지 나타냄
dim(prostate.filtered[[1]][[1]])#그래서 개수가 줄어듬 

#--------------------MetaQC------------------------------------------------------------------------#
library(MetaQC)
Data.QC<-list() #빈 리스트를 만듬 

for(i in 1:9){
  colnames(prostate.filtered[[i]][[1]])<-prostate.filtered[[i]][[2]] #첫번쨰 column이 normal인데 그걸 없애고 싶은듯 
  Data.QC[[i]]<-impute.knn(prostate.filtered[[i]][[1]])$data
  #impute.knn함수-가장 가까운 인접 평균을 사용하여 결측 표현식 데이터를 귀속시키는 함수.

  print(dim(Data.QC[[1]]))
}
names(Data.QC)<-names(prostate.filtered)
start<-Sys.time()#날짜랑 현재시간을 나타냄 
ProstateQC<-MetaQC(Data.QC, "c2.all.v3.0.symbols.gmt", filterGenes=F, verbose=TRUE,isParallel=FALSE, resp.type="Twoclass")
#metaqc 말그대로 qc를 진행해주는 것
#첫번째 요소-각 데이터는 각각 표현 데이터, 결과 또는 클래스 라벨, 유전자 ID를 나타내는 x, y 및 geneid(매트릭스 x의 행 이름으로 유전자를 바꿀 수 있음) 요소를 가져야 하는 목록으로 표시된다. 
#gene symbol lists를 가지고 있는 파일의 위치램
#filtergenes-gene filtering을 쓸건지 말건지
#verbose-로그 정보를 인쇄할지 말지
#isparallel - 빨리하기위해 코어를 여러게 쓸건지
#ncores-위에거 true면 그래서 core를 몇개 쓸건지
#resp.type- type of response variable 

runQC(ProstateQC, B=1e4, fileForCQCp="c2.all.v3.0.symbols.gmt")

#B-permutation test횟수

QC_time<-Sys.time()-start
png(filename = "Prostate_QC0421.png", width = 3500, height = 3500,res=600)
plot(ProstateQC)
dev.off()
jpeg(filename = "Prostate_QC0421.jpeg", width = 3500, height = 3500,res=600)
plot(ProstateQC)
dev.off()

#여기서 png를 하면 뒤에 나올 plot()을 png 형태로 저장하겠다는거고 dev.off하면 여기까지 나온 plot을 모두 하나의 png 형태로 저장하는거야
#jpeg도 png랑 같은 원리
#png plot plot dev.off 하면 plot두개가 동시에 png로 만들어지는거지
#그래서 중간에 plot()을 해도 눈에 안보여주는데 저 순서로 하지말고 그냥 plot만 실행하면 그림 보여줌 

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
names(data.QC.raw)#null 상태임
names(data.QC.raw)<-study.names
names(data.QC.raw)#이름이 들어감 

data.QC.merged<-MetaDE.merge(data.QC.raw)
dim(data.QC.merged[[1]][[1]])
data.QC.filtered<-MetaDE.filter(data.QC.merged,c(0.2,0.2))
#low expression gene들이랑 test control간에 variation이 적은 놈들 거름
dim(data.QC.filtered[[1]][[1]])
#---------------------- MetaDE-----------------------------------------------------------------#
start<-Sys.time()
MetaDE.Res<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Fisher"),nperm=300,asymptotic=F,ind.tail = "abs")
#합치려고 만드는 임의의 metade.res임 

MetaDE.fisherH<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Fisher"),nperm=300,asymptotic=F,ind.tail = "high")
MetaDE.fisherL<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Fisher"),nperm=300,asymptotic=F,ind.tail = "low")

MetaDE.stouffer<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("Stouffer"),nperm=300,asymptotic=F)
MetaDE.fem<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("FEM"),nperm=300,asymptotic=F)
MetaDE.rem<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),paired=rep(FALSE,6),meta.method=c("REM"),nperm=300,asymptotic=F)




#modt-moderate t test/ two sample t test랑 비슷한 거 더 잘 맞도록 수정한 현대적인 개념일걸?
#nperm은 permutation test 그거를 300번 하겠다는 것 
#ind->independentㅎ../ ind.method는 통계적 방법을 말하는 듯
#meta.method는 meta analysis방법을 말하는 거고 교수님이 돌려보길 바라는 게 저 4가지임
#rth는 rop나 rop.oc 방법 쓸 때 필요한듯 
#asymptotic - 메타 분석에서 p-값을 계산하기 위해 파라메트릭 방법을 선택했는지 여부를 지정하는 논리적 값. 기본값은 FALSE입니다.


'원래 여기서 MetaDE.Res1<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",6),meta.method=c("Fisher","sto~","fem","rem"),nperm=300,asymptotic=F)
 이렇게 돌렸는데 이러니까 paired를 설정해달라는 오류가 떴어 
 왜 그러냐면 fem이랑 rem은 이 데이터들이 대응표본인지 아닌지를 구분해야한대 그래서 이게 그런지 아닌지를 paired=false이런식으로 알려줘야함
 그 다음에 이걸 4개로 쪼갠 이유는 어떤 meta방법들은 같이 돌려도 잘 돌아갔는데 아닌 것들이 있어서 따로 돌림
 근데 이렇게 해도 오류가 나 
 안에 패키지에서 low 수를 가져와야하는데 column수를 가져오는? 그런 잘 못된 게 있더래 
 그래서 패키지 다시 고치고 설치함 
 MetaDE_revised 어쩌고 
'

#min 고르기 
'fisher는 단측 검정을 사용하기 때문에 양쪽으로 해준 후에 알파/2값이 더 작은 것을 찾아서 그걸 두배해준다'

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

#이름합치기
# category<-list()
category <-c(names(MetaDE.rem$meta.analysis),names(MetaDE.fisher$meta.analysis),names(MetaDE.fem$meta.analysis),names(MetaDE.stouffer$meta.analysis))
category<-unique(category)
category

#list에 접근하려면 [[]]/ vector는 c이거고 둘다 for문으로 접근가능함
#list는 index번호로도 접근 가능하고 그 네임으로도 접근 가능 대신 [[번호/이름]]이렇게 해야함! 

MetaDE.Res$meta.analysis[['pval']]
#새로운 metaanalysis
for (i in category){
  MetaDE.Res$meta.analysis[[i]]<-cbind(MetaDE.fisher$meta.analysis[[i]],MetaDE.rem$meta.analysis[[i]],MetaDE.fem$meta.analysis[[i]],MetaDE.stouffer$meta.analysis[[i]])
  # colnames(MetaDE.Res$meta.analysis[[i]])<-c("Fisher","REM","FEM","Stouffer")[which(c(is.null(MetaDE.fisher$meta.analysis[category[i]]),is.null(MetaDE.Res.REM$meta.analysis[[i]]),is.null(MetaDE.Res.FEM$meta.analysis[[i]]))==F),is.null(MetaDE.stouffer$meta.analysis[category[i]])]
  colnames(MetaDE.Res$meta.analysis[[i]])<-c("Fisher","REM","FEM","Stouffer")[which(c(is.null(MetaDE.fisher$meta.analysis[[i]]),is.null(MetaDE.rem$meta.analysis[[i]]),is.null(MetaDE.fem$meta.analysis[[i]]),is.null(MetaDE.stouffer$meta.analysis[[i]]))==F)]
}

'which(c(1,2,3,4,5,1)==1) 이렇게 하면 저 벡터에서 조건을 만족하는 위치를 반환해 그럼 1이랑 6이겠지
지금 for문 안을 보면 c(f,r,f,s)[which(c(조건))]이렇게 되있잖아
이때 각 방법의 meta.analysis에 pbval가 있는지 없는지 본다음에 is null함수에 조건이 f니까 채워져있는 위치를 반환해
1,2,3,4 / 1,3/ 기타 등등 그럼 이 index를 앞에 있는 c(f,r,f,s)여기 위치에 맞게 반환하는거야
fisher는 돌린 결과 값에 colname이 있고 rem/fem은 없거든 그래서 그림을 그려도 누구는 이름이 없고 누구는 없었던 거야
그래서 이걸 보정해주려고 열 이름을 지정해준거지! 

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
#여기도 위에 처럼 4개각각 돌리고 다시 저장하고 합치고 그런 과정을 거쳐야해 
b<-Sys.time()-startprint(b)
mylty<-rep(c(1,2),c(9,1))
mycol<-c(rep("black",6),c("red","pink","blue","yellow"))
mylwd<-rep(c(1,2),c(9,1))
mypch<-1:10
png(filename = "MetaDE_Prostate0421.png", width = 3500, height = 3500,res=600)
draw.DEnumber(MetaDE.Res,0.05,mlty=mylty,mcol=mycol,mlwd=mylwd,mpch=mypch,FDR=T)
dev.off()
