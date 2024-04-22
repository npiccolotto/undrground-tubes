library(pwr)
library(MASS)
library(ggplot2)
library(car)
library(psych)

pwr.anova.test(k=7, f=0.25, sig.level = 0.05, power=0.8)
pwr.t.test(d=0.75, sig.level = 0.05, power=0.8)



dataset <- read.csv2("/home/markus/PycharmProjects/ensemble-set-rendering/eval/survey_long.csv", colClasses = c("numeric","numeric","factor","factor","factor","factor"), sep=",", dec = ".")

dataset[dataset$metric == "t",]

dataset_t <- subset(dataset, dataset$metric == "t")
#dataset_t[dataset_t$task == "1",]

plot_dist <- function(df, c=value, binwidth=3000) {
  # Histogram overlaid with kernel density curve
  ggplot(df, aes(x=value)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=binwidth,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
}

plot_dist_trans <- function(df, c=value, binwidth=3000) {
  # Histogram overlaid with kernel density curve
  ggplot(df, aes(x=value_t)) +
    geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                   binwidth=binwidth,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
}


##### Task 1
dat <- dataset_t[dataset_t$task == "4",]

library(moments)
skewness(dat$value)

### box cox transform
b <- boxcox(lm(value ~ 1, data = dat))
lambda <- b$x[which.max(b$y)]
lambda
lambda <- 0.5

qq_data_bc <- 1 / sqrt(qq_data)

dat$logvalue <- sqrt(dat$value)

head(dat)
plot_dist(dat, c=logvalue)

ggplot(dataset_t, aes(x=task, y=value, fill=task)) + geom_boxplot() +
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


dataset_task1 <- subset(dataset_t, dataset_t$task == "1")
dataset_task2 <- subset(dataset_t, dataset_t$task == "2")
dataset_task3 <- subset(dataset_t, dataset_t$task == "3")
dataset_task4 <- subset(dataset_t, dataset_t$task == "4")

qq_data <- dataset


qq_data <- dataset[(dataset$metric == "t") & (dataset$task == "4") &
                     (dataset$pipeline == "o") & (dataset$alpha == "5"),]
qqnorm(qq_data$value)
qqline(qq_data,)

## α = 0 is more accurate and faster than the alternatives for T1.
boxplot(dataset_task1$value ~ dataset_task1$alpha)
qqnorm(dataset_task1$value)
qqline(dataset_task1$value)
shapiro.test(dataset_task1$value)

plot_dist(dataset_task1)

### remove outlier
summary(dataset_task1$value)

out <- which(dataset_task1$value >= quantile(dataset_task1$value, 0.99))

plot_dist(dataset_task1[-out,])
qqnorm(dataset_task1[-out,'value'])
qqline(dataset_task1[-out,'value'])

skewness(dataset_task1[-out,'value'])
kurtosis(dataset_task1[-out,'value'])

dataset_task1$value_t <- sqrt(dataset_task1[,'value'])
plot_dist_trans(dataset_task1, binwidth = 10)

# dataset_task1$value_t <- log10(dataset_task1[,'value'])
# plot_dist_trans(dataset_task1, binwidth = 0.05)

skewness(dataset_task1[-out,'value_t'])
kurtosis(dataset_task1[-out,'value_t'])

qqnorm(dataset_task1[-out,'value_t'])
qqline(dataset_task1[-out,'value_t'])
shapiro.test(dataset_task1[-out,'value_t'])
shapiro.test(dataset_task1$value_t)

leveneTest(lm(value_t ~ alpha, data = dataset_task1))
describeBy(dataset_task1$value_t, group=dataset_task1$alpha)

boxplot(dataset_task1$value_t ~ dataset_task1$alpha)
anova(lm(value_t ~ alpha*pipeline, data = dataset_task1))

0.05*0.05

t.test(dataset_task1[dataset_task1$alpha=='0','value_t'],
       dataset_task1[dataset_task1$alpha=='1','value_t'],
        var.equal = T, alternative = "less")
t.test(dataset_task1[dataset_task1$alpha=='0','value_t'],
       dataset_task1[dataset_task1$alpha=='5','value_t'],
        var.equal = T, alternative = "less")

##wilcox.test(value_t~alpha, data = dataset_task1[dataset_task1$alpha,])

## α = 1 is faster than the alternatives for T2.
boxplot(dataset_task2$value ~ dataset_task2$alpha)
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
