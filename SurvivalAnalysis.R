library(tidyverse)
library(survival)
library(ggfortify)
library(ranger)
library(boot)

setwd("D:/Programowanie/R Projectss/GEOD")

b <- read.table("GSM268274_sample_table.txt",header = TRUE,sep = "\t" )
list.files <- list.files()
data <- data.frame()
for(i in 1:length(list.files)){
  d <- read.table(list.files[[i]],header = TRUE,sep = "\t")
  data1 <- d[,2]
  data <- rbind(data,data1)
}
rownames(data) <- c(list.files)
colnames(data) <- b[,1]

row.names <- rownames(data)
new.name <- list()


for(i in 1:length(row.names)){
  list <- list()
  name <- as.character(row.names[i])
  list <- as.list(strsplit((name), "_")[[1]])
  name <- list[1]
  new.name <- append(new.name , name)
}


row.names(data) <- c(new.name)
save(data , file = "D:/Programowanie/R Projectss/FinalData/data.processed.RDAta")


clinical.data <- read.delim("D:/Programowanie/R Projectss/FinalData/a.txt", header = TRUE, sep = "\t", quote = '"', dec = ".", fill = TRUE, comment.char = "")
status <- as.character(clinical.data[,22])
row.names.clinical <- as.character(clinical.data[,33])

status <- as.numeric(clinical.data[,22])
age.at.RRP <- as.numeric(clinical.data[,5])
stage <- as.numeric(clinical.data[,54])
ploidy <- as.integer(clinical.data[,10])
gleason.score <- as.numeric(clinical.data[,13])
control.group.var <- as.numeric(clinical.data[,6])
years.of.observation.after.RRP.var <- as.integer(clinical.data[,18])

table <- as.data.frame(cbind(row.names.clinical,age.at.RRP,years.of.observation.after.RRP.var,control.group.var, stage, gleason.score, ploidy, status))
make.intersect.patients <- intersect(as.character(row.names(data)) , as.character(table[,1]))

unique.data.var <- table[!duplicated(table[,1]),]
table <- unique.data.var[1:260,]
row.names(table) <- table[,1]
table$row.names.clinical <- NULL
data <- data[row.names(table),]

data$age.at.RRP <- table[,1]
data$years.of.observation.after.RRP.var <- table[,2]
data$control.group.var <- table[,3]
data$stage <- table[,4]
data$gleason.score <- table[,5]
data$ploidy <- table[,6]
data$status <- table[,7]


data.to.analyse <- data


data.to.analyse[,1:533] <- lapply(data.to.analyse[,1:533], as.character)
data.to.analyse[,1:533] <- lapply(data.to.analyse[,1:533], as.integer)

age.divided.in.half <- ifelse((data.to.analyse$age.at.RRP <= 63), "under.63", "over.63")
length(age.divided.in.half)

names(data.to.analyse)[528]<-"days.of.observation.after.RRP"
data.to.analyse[,528] <- data.to.analyse[,528] * 365

  colname.var <- as.character(colnames(data.to.analyse))
  new.names.var <- str_replace_all(colname.var, c("-|_"), ".")
  
  
  new.names.var1 <- str_replace(new.names.var, c("1557685"), "GIA")
  new.names.var2 <- str_replace(new.names.var1, c("1560225"), "GIB")
  new.names.var3 <- str_replace(new.names.var2, c("1561073"), "GIC")
  new.names.var4 <- str_replace(new.names.var3, c("213310"), "GID")
  new.names.var5 <- str_replace(new.names.var4, c("214174"), "GIE")
  new.names.var6 <- str_replace(new.names.var5, c("213310"), "GII")
  new.names.var7 <- str_replace(new.names.var6, c("214174"), "GIJ")
  new.names.var8 <- str_replace(new.names.var7, c("214384"), "GIK")
  new.names.var9 <- str_replace(new.names.var8, c("216584"), "GIL")
  new.names.var10 <- str_replace(new.names.var9, c("225311"), "GIM")
  new.names.var11 <- str_replace(new.names.var10, c("228178"), "GIN")
  new.names.var12 <- str_replace(new.names.var11, c("228687"), "GIO")
  new.names.var13 <- str_replace(new.names.var12, c("228798"), "GIP")
  new.names.var14 <- str_replace(new.names.var13, c("229455"), "GIR")
  new.names.var15 <- str_replace(new.names.var14, c("231181"), "GIS")
  new.names.var16 <- str_replace(new.names.var15, c("231651"), "GIT")
  new.names.var17 <- str_replace(new.names.var16, c("233424"), "GIU")
  new.names.var18 <- str_replace(new.names.var17, c("235196"), "GIX")
  new.names.var19 <- str_replace(new.names.var18, c("236338"), "GIY")
  new.names.var20 <- str_replace(new.names.var19, c("236472"), "GIZ")
  new.names.var21 <- str_replace(new.names.var20, c("240093"), "GIZ")
  new.names.var22 <- str_replace(new.names.var21, c("240236"), "GIZA")
  new.names.var23 <- str_replace(new.names.var22, c("240543"), "GIZB")
  new.names.var24 <- str_replace(new.names.var23, c("241944"), "GIZC")
  new.names.var25 <- str_replace(new.names.var24, c("242763"), "GIZD")
  new.names.var26 <- str_replace(new.names.var25, c("243636"), "GIZE")
  new.names.var27 <- str_replace(new.names.var26, c("244579"), "GIZF")

  
  colnames(data.to.analyse) <- new.names.var27



complete.cases(data.to.analyse)
data.to.analyse <- data.to.analyse[complete.cases(data.to.analyse), ]

age.divided.in.half <- ifelse((data.to.analyse$age.at.RRP <= 63), "under.63", "over.63")


keplan.meier <- with(data.to.analyse, Surv(days.of.observation.after.RRP, status))
head(keplan.meier,80)

keplan.meier.fit <- survfit(Surv(days.of.observation.after.RRP, status) ~ 1 , data=data.to.analyse)
summary(keplan.meier.fit)
autoplot(keplan.meier.fit)


keplan.meier.fit.cg <- survfit(Surv(days.of.observation.after.RRP, status) ~ control.group.var, data=data.to.analyse) 
summary(keplan.meier.fit.cg)
autoplot(keplan.meier.fit.cg)



keplan.meier.fit.st <- survfit(Surv(days.of.observation.after.RRP, status) ~ stage, data=data.to.analyse)
summary(keplan.meier.fit.st)
autoplot(keplan.meier.fit.st)



keplan.meier.fit.age <- survfit(Surv(days.of.observation.after.RRP, status) ~ age.divided.in.half, data=data.to.analyse)
summary(keplan.meier.fit.age)
autoplot(keplan.meier.fit.age)



attributes.comparison.fit <-aareg(Surv(days.of.observation.after.RRP, status) ~ control.group.var + stage + ploidy  + age.at.RRP , 
                                  data = data.to.analyse)
attributes.comparison.fit
autoplot(attributes.comparison.fit)

cox.model <- coxph(Surv(days.of.observation.after.RRP, status) ~  control.group.var + stage + age.at.RRP, data = data.to.analyse)

summary(cox.model)

cox.model.fit <- survfit(cox.model)
autoplot(cox.model.fit)


random.forest.fit <- ranger(Surv(days.of.observation.after.RRP, status) ~ . ,
                data = data.to.analyse,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)

death.times <- random.forest.fit$unique.death.times 
surv.prob <- data.frame(random.forest.fit$survival)
avg.prob <- sapply(surv.prob,mean)


plot(random.forest.fit$unique.death.times,random.forest.fit$survival[1,], 
     type = "l", 
     ylim = c(0,1),
     col = "red",
     xlab = "Dni",
     ylab = "szanse przeżycia",
     main = "Krzywa Przeżycia Pacjentów")



var.importance <- data.frame(sort(round(random.forest.fit$variable.importance, 4), decreasing = TRUE))
names(var.importance) <- "importance"
head(var.importance)

cat("Prediction Error = 1 - Harrell's c-index = ", random.forest.fit$prediction.error)


kmi <- rep("KM",length(keplan.meier.fit$time))
km.df <- data.frame(keplan.meier.fit$time,keplan.meier.fit$surv,kmi)
names(km.df) <- c("Time","Surv","Model")

coxi <- rep("Cox",length(cox.model.fit$time))
cox.df <- data.frame(cox.model.fit$time,cox.model.fit$surv,coxi)
names(cox.df) <- c("Time","Surv","Model")

randomfi <- rep("Random Forest",length(random.forest.fit$unique.death.times))
randomfdf <- data.frame(random.forest.fit$unique.death.times,avg.prob,randomfi)
names(randomfdf) <- c("Time","Surv","Model")

plot.df <- rbind(km.df,cox.df,randomfdf)

p <- ggplot(plot.df, aes(x = Time, y = Surv, color = Model))
p + geom_line()


summary(keplan.meier.fit)
summary(cox.model.fit)
summary(random.forest.fit)


log.rank.age.at.RRP <- survdiff(Surv(days.of.observation.after.RRP, status) ~ age.at.RRP ,data=data.to.analyse)
log.rank.age.at.RRP


log.rank.control.group.var <- survdiff(Surv(days.of.observation.after.RRP, status) ~ control.group.var, data = data.to.analyse)
log.rank.control.group.var

log.rank.stage<- survdiff(Surv(days.of.observation.after.RRP, status) ~ stage, data = data.to.analyse)
log.rank.stage

boootstrap.function <- function(data, indices, cor.type){
  dt<-data[indices,]
  c(
    cor(dt[,529],dt[,531],  method=cor.type),
    median(dt[,529]),
    median(dt[,531])
    
  )
}


bootstrapping <- boot(data.to.analyse, boootstrap.function, R=1000, cor.type='s')


plot(bootstrapping, index=1)


table.of.indices<-boot.array(bootstrapping, indices=T)

table.of.appearances<-boot.array(bootstrapping)

bootstrap.summary<-apply(table.of.indices, 1, boootstrap.function, data=data.to.analyse, cor.type='s')


all(t(bootstrap.summary)==bootstrapping$t)
