english <-c(90,80,90,70)
math <-c(50,60,100,20)
student <-c("a","b","c","d")
class <-c(1,1,2,2)

midterm <-cbind(english,math,class) #column방향으로 합치는 것 
midterm
class(midterm)#데이터타입 확인

rownames(midterm)<-student #rownames를 지정안해주면 숫자로 나옴

midterm

midterm[2,2]
midterm[,-1]#1열만 빼고 부름 
midterm[,1]


midterm <-cbind(english,math,class,student)
midterm
#student를 cbind로 합치면 갑자기 모든 데이터에 ""이거 생김
#cbind는 같은 데이터 형식으로 묶어서 그렇게 됨
#data frame을 쓰면 다른 데이터 형식 끼리 묶을 수 있음 

midterm.df=data.frame(english,math,class,student)
midterm.df

#selection
midterm.df$english
midterm.df$math
midterm.df$student #자동으로 factor로 저장함

midterm.df$student <-as.character(student)
midterm.df$student #그래서 이걸 as. character를 써서 문자로 바꿔줘야함 

#df형식이면 새로운 열 추가 하는 게 쉬움 / df이름$열 이름 지정 <- 새로운 거 
midterm.df$newcolumn <-c(1,1,1,1)
midterm.df

#
dim(midterm.df) #행열 크기를 나타냄 몇 행 몇 열인지
summary(midterm.df)
View(midterm.df)

#vector,factor [] 
english[c(1,3)] #1,3 이렇게 하면 에러남 각행의 정보 가지고 오고 싶으면 c로 묶어야함 
#english[1,3]#error??
#list는 [[]]이렇게 해서 가져옴/ list는 길이가 다른 벡터들의 조합임
a <- c(1,2,3)
b <-c("a")
c<-c(T,F)

list1 <-list(a,b,c)
list2 <-list(numeric=a,charater=b,logi=c)
list2

list2$numeric
list2[[1]] #a 가 선택됨 

#df
midterm.df$student[midterm.df$math>=60]

#문자열 합치기 
str <-c("i","like","you")
paste(str, collapse = ",")

#
data(iris) #iris r에서 제공하는 데이터! 이걸 가져오는 작업 
class(iris)
dim(iris)
summary(iris)


#
mean(iris$Sepal.Width[iris$Species=="setosa"])
max(iris$Sepal.Width[iris$Species=="setosa"])
min(iris$Sepal.Width[iris$Species=="setosa"])

#if

x<--3
if (x<0) {
  print("x is a negative number")
}

x<-5
if (x<0){
  print("x is a negative number")
} else {
  print( "x is either a positive number or zero")
}


x<-2

if (x <0) {
  print("negative")
}else if (x ==0){
  print("x is zero")
}else{
  print("positive")
}


medium <-"instagram"
num_vies <-14

if (medium =="instagram"){
  print("showing instagram information")
}else{
  print ("unknown medium")
}

if (num_vies>15){
  print("you're popular")
}else{
  print("try to be more viribele")
}

#for
cities <-c("new yotk","paris","london","tokyo","rio de","cape town")

for (city in cities){
  print(city)
}

for(i in 1:length(cities)){
  print(paste(cities[i], "is on position", i, "in the cities vector"))
}

#
instagram <-c(16,9,13,5,2,17,14)

for(li in instagram){
  if (li >10){
    print("you're popular")
  }
}
# for 문 대신해서 쓸 수 있는 apply family

nyc <- list(pop = 8405837, boroughs = c("Manhattan", "Bronx", "Brooklyn", "Queens", "Staten Island"), capital = FALSE)

lapply(nyc, class)
#lapply(list, 수행하고 싶은 function)

#
num_chars <-c()
for(i in 1:length(cities)){
  num_chars[i] <-nchar(cities[i])
}
num_chars
#nchar- 문자개수 세는 거 

lapply(cities, nchar) #list로 결과 값 반환
unlist(lapply(cities, nchar)) #list 풀어주는 것 

#
pioneers <- c("GAUSS:1777", "BAYES:1702", "PASCAL:1623", "PEARSON:1857")
split <-strsplit(pioneers, split=":")
split_low <-lapply(split, tolower)
split_low

#함수만들기
'함수명 <- function(x ){
  x[1] 수행하고 싶은 거 
}'

select_first <- function(x){
  x[1]
}

names <-lapply(split_low, select_first)
names 

select_second <- function(x){
  x[2]
}

years <- lapply(split_low, select_second)
years


lapply(split_low, function(x) x[1]) #한줄로 함수도 만들고 for문도 돌리고 할 수 있는 것 

#argument 2개 이상일떄
lapply(list(1,2,3), log, base=2)
#log함수는 인자가 두갠데 그냥 쓰면 됨 , 찍고 쭉쭉

#sapply lapply랑 비슷 결과 값을 벡터로 반환
#unlist할 필요 없음 





























