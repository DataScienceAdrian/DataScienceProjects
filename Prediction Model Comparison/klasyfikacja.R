library(rpart)
library(caret)
library(e1071)
library(class)
library(randomForest)
data2 <- read.csv("heart.csv")


#Podział na train i test
rows <- sample.int(nrow(data2), size =
                     round(nrow(data2)/3), replace = F)


data.train <- data2[-rows,]
data.test = data2[rows,]

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}


data.train$target <- as.factor(data.train$target)
##############
#Klasyfikacja KNN
##############

prc_test_pred <- knn(train = data.train, test = data.test, cl = data.train$target, k=15, prob=TRUE)

knn.model <- table(data.test$target, prc_test_pred)

accuracy(knn.model)


############
#Klasyfikacja Random Forest
############
model1 <- randomForest(target ~ ., data = data.train, ntree = 500, mtry = 6, importance = TRUE)
model1



predTrain <- predict(model1, data.train, type = "class")
table(predTrain, data.train$target) 
importance(model1) 





model_dt = train(target ~ ., data = data.train, method = "rpart")
model_dt_1 = predict(model_dt, data = data.train)
rf.model <- table(model_dt_1, data.train$target)
accuracy(rf.model)


model_dt_vs = predict(model_dt, newdata = data.test)
rf.model2 <- table(model_dt_vs, data.test$target)
accuracy(rf.model2)