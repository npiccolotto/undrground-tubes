library(pwr)
library(MASS)
library(ggplot2)
library(car)
library(psych)
library(plyr)
library(moments)

pwr.anova.test(k=7, f=0.25, sig.level = 0.05, power=0.8)
pwr.t.test(d=0.75, sig.level = 0.05, power=0.8)


## read survey data to dataset
read_dataset <- read.csv2("/home/markus/PycharmProjects/ensemble-set-rendering/eval/survey_long.csv", colClasses = c("numeric","numeric","factor","factor","factor","factor"), sep=",", dec = ".")

# split time and accuracy in separate datasets
dataset_t <- subset(read_dataset, read_dataset$metric == "t")
dataset_a <- subset(read_dataset, read_dataset$metric == "a")

# join time and accuracy
dat <- join(dataset_t, dataset_a, by=c('participant', 'task', 'alpha', 'pipeline'))
# remove metric column
dat <- dat[,c(-6,-8)]

# rename metric column
colnames(dat)[2] <- "time"
colnames(dat)[6] <- "f1.error"
# set dataset
dataset <- dat

dataset_task1 <- subset(dataset, dataset$task == "1")
dataset_task2 <- subset(dataset, dataset$task == "2")
dataset_task3 <- subset(dataset, dataset$task == "3")
dataset_task4 <- subset(dataset, dataset$task == "4")

plot_dist <- function(df, column="time", binwidth=3000) {
  # Histogram overlaid with kernel density curve
  ggplot(df, aes_string(x=column)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=binwidth,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
}

##### Task 1
describeBy(dataset$time, group=dataset$task)
describeBy(dataset$f1.error, group=dataset$task)

## α = 0 is more accurate and faster than the alternatives for T1.
boxplot(dataset_task1$time ~ dataset_task1$alpha)
qqnorm(dataset_task1$time)
qqline(dataset_task1$time)
shapiro.test(dataset_task1$time)

plot_dist(dataset_task1)

boxplot(dataset_task1$f1.error ~ dataset_task1$alpha)
describeBy(dataset_task1$f1.error, group=dataset_task1$alpha)
qqnorm(dataset_task1$f1.error)
qqline(dataset_task1$f1.error)
shapiro.test(dataset_task1$f1.error)

dataset_task1$time_adj <- sqrt(dataset_task1[,'time'])
plot_dist(dataset_task1, "time_adj", binwidth = 10)

dataset_task1$f1_adj <- 1/(max(dataset_task1[,'f1.error']+1)-dataset_task1[,'f1.error'])
plot_dist(dataset_task1, "f1_adj", binwidth = 0.01)

par(mfrow=c(2,1))
qqnorm(dataset_task1$time)
qqline(dataset_task1$time)
qqnorm(dataset_task1$time_adj)
qqline(dataset_task1$time_adj)
shapiro.test(dataset_task1$time_adj)

leveneTest(lm(time_adj ~ alpha, data = dataset_task1))
describeBy(dataset_task1$time_adj, group=dataset_task1$alpha)

boxplot(dataset_task1$time_adj ~ dataset_task1$alpha)
anova(lm(time_adj ~ alpha*pipeline, data = dataset_task1))

0.05*0.05

t.test(dataset_task1[dataset_task1$alpha=='0','time_adj'],
       dataset_task1[dataset_task1$alpha=='1','time_adj'],
        var.equal = T, alternative = "less")
t.test(dataset_task1[dataset_task1$alpha=='0','time_adj'],
       dataset_task1[dataset_task1$alpha=='5','time_adj'],
        var.equal = T, alternative = "less")

##wilcox.test(value_t~alpha, data = dataset_task1[dataset_task1$alpha,])

######################################################
## α = 1 is faster than the alternatives for T2.
par(mfrow=c(1,1))
boxplot(dataset_task2$time ~ dataset_task2$alpha)
qqnorm(dataset_task2$time)
qqline(dataset_task2$time)
plot_dist(dataset_task2)
describeBy(dataset_task2$time, group=dataset_task2$alpha)
shapiro.test(dataset_task2$time)

b <- boxcox(lm(time~1, data=dataset_task2))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task2$time_adj <- log10(dataset_task2[,'time'])
plot_dist(dataset_task2, "time_adj", binwidth = 0.1)
shapiro.test(dataset_task2$time_adj)

dataset_task2$time_adj <- 1/sqrt(dataset_task2$time/10000)
plot_dist(dataset_task2, "time_adj", binwidth = 0.1)
shapiro.test(dataset_task2$time_adj)

qqnorm(dataset_task2$time_adj)
qqline(dataset_task2$time_adj)

leveneTest(lm(time_adj ~ alpha, data = dataset_task2))
describeBy(dataset_task2$time_adj, group=dataset_task2$alpha)

boxplot(dataset_task2$time_adj ~ dataset_task2$alpha)
anova(lm(time_adj ~ alpha*pipeline, data = dataset_task2))

t.test(dataset_task2[dataset_task2$alpha=='0','time_adj'],
       dataset_task2[dataset_task2$alpha=='1','time_adj'],
        var.equal = T)
t.test(dataset_task2[dataset_task2$alpha=='1','time_adj'],
       dataset_task2[dataset_task2$alpha=='5','time_adj'],
        var.equal = T)

#####################################################################
## H3: α = 1 is more accurate and faster than the alternatives for T3.
boxplot(dataset_task3$time ~ dataset_task3$alpha)
qqnorm(dataset_task3$time)
qqline(dataset_task3$time)
plot_dist(dataset_task3)
describeBy(dataset_task3$time, group=dataset_task3$alpha)
shapiro.test(dataset_task3$time)

b <- boxcox(lm(time~1, data=dataset_task3))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task3$time_adj <- 1/sqrt(dataset_task3$time/10000)
plot_dist(dataset_task3, "time_adj", binwidth = 0.1)
shapiro.test(dataset_task3$time_adj)

qqnorm(dataset_task3$time_adj)
qqline(dataset_task3$time_adj)

leveneTest(lm(time_adj ~ alpha, data = dataset_task3))
describeBy(dataset_task3$time_adj, group=dataset_task3$alpha)

boxplot(dataset_task3$time_adj ~ dataset_task3$alpha)
anova(lm(time_adj ~ alpha*pipeline, data = dataset_task3))

anova(lm(time_adj ~ alpha, data = dataset_task3[dataset_task3$pipeline=="o",]))
anova(lm(time_adj ~ alpha, data = dataset_task3[dataset_task3$pipeline=="h",]))

par(mfrow=c(2,1))
boxplot(dataset_task3$time_adj ~ dataset_task3$alpha, data=dataset_task3[dataset_task3$pipeline=="o",])
boxplot(dataset_task3$time_adj ~ dataset_task3$alpha, data=dataset_task3[dataset_task3$pipeline=="h",])

t.test(dataset_task3[dataset_task3$alpha=='0','time_adj'],
       dataset_task3[dataset_task3$alpha=='1','time_adj'],
        var.equal = T)

t.test(dataset_task3[dataset_task3$alpha=='1','time_adj'],
       dataset_task3[dataset_task3$alpha=='5','time_adj'],
        var.equal = T)

t.test(dataset_task3[dataset_task3$alpha=='0','time_adj'],
       dataset_task3[dataset_task3$alpha=='5','time_adj'],
        var.equal = T)


boxplot(dataset_task3$f1.error ~ dataset_task3$alpha)

##########################################################
## H4: α = 0 is more accurate and faster than the alternatives for T4
boxplot(dataset_task4$time ~ dataset_task4$alpha)
qqnorm(dataset_task4$time)
qqline(dataset_task4$time)
plot_dist(dataset_task4)
describeBy(dataset_task4$time, group=dataset_task4$alpha)
shapiro.test(dataset_task4$time)

b <- boxcox(lm(time~1, data=dataset_task4))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task4$time_adj <- sqrt(dataset_task4$time)
plot_dist(dataset_task4, "time_adj", binwidth = 10)
shapiro.test(dataset_task4$time_adj)

qqnorm(dataset_task4$time_adj)
qqline(dataset_task4$time_adj)

leveneTest(lm(time_adj ~ alpha, data = dataset_task4))
describeBy(dataset_task4$time_adj, group=dataset_task4$alpha)

boxplot(dataset_task4$time_adj ~ dataset_task4$alpha)
anova(lm(time_adj ~ alpha*pipeline, data = dataset_task4))

0.05*0.05

t.test(dataset_task4[dataset_task4$alpha=='0','time_adj'],
       dataset_task4[dataset_task4$alpha=='1','time_adj'],
        var.equal = T)

t.test(dataset_task4[dataset_task4$alpha=='0','time_adj'],
       dataset_task4[dataset_task4$alpha=='5','time_adj'],
        var.equal = T)

##########################################################
## H5: α = 0.5 is neither the slowest nor least accurate for any task


##########################################################
## H6: Optimal UT drawings lead to faster answers for all tasks.


############################################

qq_data <- subset(dataset_task2$value, dataset_task2$alpha == "1")
qq_data <- dataset_task2$value
qqnorm(qq_data)
qqline(qq_data)
shapiro.test(qq_data)

### box cox transform
b <- boxcox(lm(value ~ 1, data = dataset_task2))
lambda <- b$x[which.max(b$y)]
lambda
lambda <- -0.5

qq_data_bc <- 1 / sqrt(qq_data)
qqnorm(qq_data_bc)
qqline(qq_data_bc)
shapiro.test(qq_data_bc)

dataset_task2$boxcoxvalue <- qq_data_bc

anova(lm(value ~ alpha*pipeline, data = dataset_task2))

anova(lm(value ~ alpha, data = dataset_task2[dataset_task2$pipeline == "o",]))
anova(lm(value ~ alpha, data = dataset_task2[dataset_task2$pipeline == "h",]))

## α = 1 is more accurate and faster than the alternatives for T3.
boxplot(dataset_task3$value ~ dataset_task3$alpha)

qq_data <- subset(dataset_task3$value, dataset_task3$alpha == "1")

hist(qq_data, freq=F)
fit<-fitdistr(qq_data,"log-normal")$estimate
lines(dlnorm(0:max(qq_data),fit[1],fit[2]), lwd=3)

qq_data <- dataset_task3$value
qq_data <- log(qq_data)
qqnorm(qq_data)
qqline(qq_data)
shapiro.test(qq_data)

anova(lm(value ~ alpha*pipeline, data = dataset_task3))


anova(lm(value ~ alpha, data = dataset_task3[dataset_task3$pipeline == "o",]))
anova(lm(value ~ alpha, data = dataset_task3[dataset_task3$pipeline == "h",]))

### check all data
qq_data <- dataset_t$value
qqnorm(qq_data)
qqline(qq_data)

shapiro.test(qq_data)



anova(lm(value ~ pipeline*task, data = dataset_t))

describeBy(dataset_task4$value, group=dataset_task4$pipeline)







##########################################
#############################################
### box cox transform
b <- boxcox(lm(value ~ 1, data = dat))
lambda <- b$x[which.max(b$y)]
lambda
lambda <- 0.5

qq_data_bc <- 1 / sqrt(qq_data)

dat$logvalue <- sqrt(dat$value)

head(dat)
plot_dist(dat, c=logvalue)

ggplot(dataset, aes(x=task, y=time, fill=task)) + geom_boxplot() +
    guides(fill=FALSE)


anova(lm(value ~ alpha*pipeline, data = dataset_t[dataset_t$task == "4",]))

anova(lm(value ~ task*alpha*pipeline, data = dataset[dataset$metric == "t",]))

anova(lm(value ~ pipeline*task*alpha, data = dataset[dataset$metric == "t",]))

boxplot((value/1000) ~ task, data =  dataset[dataset$metric == "t",])

boxplot((value/1000) ~ pipeline, data =  dataset[dataset$metric == "t",])
boxplot((value/1000) ~ alpha + task, data =  dataset[dataset$metric == "t",])

boxplot(dataset$value ~ dataset$pipeline)

boxplot(dataset$value ~ dataset$task)

boxplot(dataset$value ~ dataset$alpha)

qq_data <- dataset


qq_data <- dataset[(dataset$metric == "t") & (dataset$task == "4") &
                     (dataset$pipeline == "o") & (dataset$alpha == "5"),]
qqnorm(qq_data$value)
qqline(qq_data,)

##############################################
## OUTLIER

### remove outlier
summary(dataset_task1$value)

out <- which(dataset_task1$value >= quantile(dataset_task1$value, 0.99))
out
plot_dist(dataset_task1[-out,])
qqnorm(dataset_task1[-out,'value'])
qqline(dataset_task1[-out,'value'])

skewness(dataset_task1[-out,'value'])
kurtosis(dataset_task1[-out,'value'])
