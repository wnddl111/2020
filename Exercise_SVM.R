#SVM e1071 exercise
#Linear case

install.packages('e1071')
library('e1071')

data1=seq(1,10,by=2)
data1
class1=factor(c('a','a','a','b','b'))
class1

model1=svm(data1,class1,type='C',kernel='linear')
?svm
model1

test1=c(1.5,8,6.2,10,2.3)

predict(model1,test1)


#Nonlinear case

data2 = seq(1,10)
class2 = factor(c('b','b','b','a','a','a','a','b','b','b'))
#Build the model with a linear kernel.
model2 = svm(data2,class2,type='C',kernel='linear')
predict(model2,data2)
table(predict(model2,data2), class2) # this function shows frequency of data

model2 = svm(data2,class2,type='C',kernel='radial')
predict(model2,data2)

model2 = svm(data2,class2,type='C',kernel='polynomial')
predict(model2,data2)


#Example
#Many Classes >=3

data(iris)
attach(iris)

x = subset(iris, select = -Species)#except species
y = Species

model = svm(x, y)

pred = predict(model, x)
pred
table(pred, y)

#Try other kernels and parameters
model = svm(x, y, kernel='radial', gamma=5)
model = svm(x, y, kernel='radial', gamma=7) #or gamma=9
p = predict(model, x)
table(p, y)

#%in% in or not in vector / return t/f / x %in% y -> in x, whether there is y or not 
#Multi-dimensional scaling
plot(cmdscale(dist(iris[,-5])),
     col = as.integer(iris[,5]),
     pch = c("o","+")[1:150 %in% model$index + 1])
# (-) means remove 5 column 
#pch: plot character, see '?points'
#col: color
#cmdscale: Classical (Metric) Multidimensional Scaling


## try regression mode on two dimensions

x = seq(0.1, 5, by = 0.05)
x
y = log(x) + rnorm(x, sd = 0.2)
y
plot(x, y)
# estimate model and predict input values
m   = svm(x, y)
new = predict(m, x)

# visualize
#plot(x, y)
points(x, log(x), col = 2)
points(x, new, col = 4)

Exercise !
##Breast cancer classification
#ER positive vs. negative patients
#perform the five-fold cross validation
#Select top 5 up and down regulated genes
#Based on the 10 dimensional data, train the SVM model
#predict the class of the test samples using the trained model
       
		     
#Optimizing SVM parameters
#tune(svm, train.x=databctrain, train.y=classesbctrain, validation.x=databctest, validation.y=classesbctest, ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)), control = tune.control(sampling = "fix"))

#Try to find the optimal parameters manually
#Try both the linear and radial basis kernel functions


## Naive Bayes Classification
data(iris)
m = naiveBayes(Species ~ ., data = iris)
## alternatively:
m = naiveBayes(iris[,-5], iris[,5])
m
table(predict(m, iris), iris[,5])
