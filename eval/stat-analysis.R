library(pwr)
library(pwr2)
library(MASS)
library(ggplot2)
library(car)
library(psych)
library(plyr)
library(moments)
library(rstatix)
library(emmeans)
library(ggpubr)

setwd('/home/markus/PycharmProjects/ensemble-set-rendering/eval/')


prepareDataset <- function(filename) {
  ## read survey data to dataset
  read_dataset <-
    read.csv2(
      filename,
      colClasses = c("numeric", "numeric", "factor", "factor", "factor", "factor"),
      sep = ",",
      dec = "."
    )

  # split time and accuracy in separate datasets
  dataset_t <- subset(read_dataset, read_dataset$metric == "t")
  dataset_a <- subset(read_dataset, read_dataset$metric == "a")

  # join time and accuracy
  dat <-
    join(dataset_t,
         dataset_a,
         by = c('participant', 'task', 'alpha', 'pipeline'))
  # remove metric column
  dat <- dat[, c(-6, -8)]

  # rename metric column
  colnames(dat)[2] <- "time"
  colnames(dat)[6] <- "f1.error"
  # set dataset
  return (dat)
}

plot_dist <- function(df,
                      column = "time",
                      binwidth = 3000) {
  # Histogram overlaid with kernel density curve
  ggplot(df, aes(x = .data[[column]])) +
    geom_histogram(
      aes(y = after_stat(density)),
      # Histogram with density instead of count on y-axis
      binwidth = binwidth,
      colour = "black",
      fill = "white"
    ) +
    geom_density(alpha = .2, fill = "#FF6666")  # Overlay with transparent density plot

  # ggplot(df, aes_string(x=column)) +
  #   geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
  #                 binwidth=binwidth,
  #                colour="black", fill="white") +
  # geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot
}


### Pre Study Analysis
pilot <-
    prepareDataset("pilot_long.csv")

pilot_task1 <- subset(pilot, pilot$task == "1")
pilot_task2 <- subset(pilot, pilot$task == "2")
pilot_task3 <- subset(pilot, pilot$task == "3")
pilot_task4 <- subset(pilot, pilot$task == "4")

describeBy(pilot_task1$time, group = pilot_task1$alpha)

pilot_task1 %>%
  cohens_d(time ~ alpha) %>%
  as.data.frame()

describeBy(pilot_task2$time, group = pilot_task2$alpha)

pilot_task2 %>%
  cohens_d(time ~ alpha) %>%
  as.data.frame()

describeBy(pilot_task3$time, group = pilot_task3$alpha)

pilot_task3 %>%
  cohens_d(time ~ alpha) %>%
  as.data.frame()

describeBy(pilot_task4$time, group = pilot_task4$alpha)

pilot_task4 %>%
  cohens_d(time ~ alpha) %>%
  as.data.frame()


## Power Analysis
ss.2way(
  a=3,
  b=2,
  alpha=0.05,
  f.A = 0.4,
  f.B = 0.2,
  beta=0.2,
  B = 100
)

pwr.2way(
  a=3,
  b=2,
  alpha=.05,
  size.A=49,
  size.B=49,
  f.A=.4,
  f.B=.2
)


## Main Study Analysis
dataset <- prepareDataset("survey_long.csv")

dataset_task1 <- subset(dataset, dataset$task == "1")
dataset_task2 <- subset(dataset, dataset$task == "2")
dataset_task3 <- subset(dataset, dataset$task == "3")
dataset_task4 <- subset(dataset, dataset$task == "4")

##### Task 1
describeBy(dataset$time, group = dataset$task)
describeBy(dataset$f1.error, group = dataset$task)

############################################################
## H1: α = 0 is more accurate and faster than the alternatives for T1.
boxplot(dataset_task1$time ~ dataset_task1$alpha)
qqnorm(dataset_task1$time)
qqline(dataset_task1$time)
shapiro.test(dataset_task1$time)

plot_dist(dataset_task1)

boxplot(dataset_task1$f1.error ~ dataset_task1$alpha)
describeBy(dataset_task1$f1.error, group = dataset_task1$alpha)
qqnorm(dataset_task1$f1.error)
qqline(dataset_task1$f1.error)
shapiro.test(dataset_task1$f1.error)

b <- boxcox(lm(time ~ 1, data = dataset_task1))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task1$time_adj <- sqrt(dataset_task1[, 'time'])
plot_dist(dataset_task1, "time_adj", binwidth = 10)

dataset_task1$f1_adj <-
  1 / (max(dataset_task1[, 'f1.error'] + 1) - dataset_task1[, 'f1.error'])
plot_dist(dataset_task1, "f1_adj", binwidth = 0.01)
boxplot(dataset_task1$f1_adj ~ dataset_task1$alpha)

par(mfrow = c(2, 1))
qqnorm(dataset_task1$time)
qqline(dataset_task1$time)
qqnorm(dataset_task1$time_adj)
qqline(dataset_task1$time_adj)
shapiro.test(dataset_task1$time_adj)

leveneTest(lm(time_adj ~ alpha * pipeline, data = dataset_task1))
describeBy(dataset_task1$time_adj, group = dataset_task1$alpha)

### check for normal distribution and homogenity of variances
par(mfrow = c(1, 1))
nv <- lm(time_adj ~ alpha * pipeline, data = dataset_task1)
ggqqplot(residuals(nv))
plot(nv, 1)

boxplot(dataset_task1$time_adj ~ dataset_task1$alpha)
anova(lm(time_adj ~ alpha * pipeline, data = dataset_task1))
anova_test(time_adj ~ alpha * pipeline,
           data = dataset_task1,
           effect.size = "pes")

### pairwise t-test
dataset_task1 %>%
  pairwise_t_test(time_adj ~ alpha,
                  pool.sd = F,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task1 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

##wilcox.test(value_t~alpha, data = dataset_task1[dataset_task1$alpha,])
describeBy(dataset_task1$f1.error, group = dataset_task1$alpha)

# test statistic for f1 score difference
# followed by https://bjoernwalther.com/kruskal-wallis-test-in-r-rechnen/
k_test <- kruskal.test(dataset_task1$f1.error ~ dataset_task1$alpha)
print(k_test)


# post hoc analysis of f1 score
#pairwise.wilcox.test(dataset_task1$f1.error,dataset_task1$alpha, p.adjust="bonferroni")
# bzw. für den Dunn's Test:
d <-
  dunn_test(f1.error ~ alpha, data = dataset_task1, p.adjust.method = "bonferroni")
print(d)
? anova
? dunn_test

kruskal.effect <- function (chi, k, n) {
  eta_sqd <- (chi - k + 1) / (n - k)
  return (sqrt(eta_sqd / (1 - eta_sqd)))
}

### calculate effect of the Kruskal Walis Test
f_value <-
  kruskal.effect(unname(k_test$statistic),
                 nlevels(dataset_task1$alpha),
                 nrow(dataset_task1))
print(f_value)
## F_value
#    Ab 0,1 ist es ein schwacher Effekt,
#    ab 0,25 ein mittlerer und
#    ab 0,4 ein starker Effekt.

### effect of dunn's test
d$r <- d$statistic / sqrt(d$n1 + d$n2)
d

######################################################
## H2: α = 1 is faster than the alternatives for T2.
par(mfrow = c(1, 1))
boxplot(dataset_task2$time ~ dataset_task2$alpha)
qqnorm(dataset_task2$time)
qqline(dataset_task2$time)
plot_dist(dataset_task2)
describeBy(dataset_task2$time, group = dataset_task2$alpha)
shapiro.test(dataset_task2$time)

#https://www.datanovia.com/en/lessons/transform-data-to-normal-distribution-in-r/
b <- boxcox(lm(time ~ 1, data = dataset_task2))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task2$time_adj <- 1 / sqrt(dataset_task2$time / 10000)
plot_dist(dataset_task2, "time_adj", binwidth = 0.1)
shapiro.test(dataset_task2$time_adj)

qqnorm(dataset_task2$time_adj)
qqline(dataset_task2$time_adj)

leveneTest(lm(time_adj ~ alpha, data = dataset_task2))
describeBy(dataset_task2$time_adj, group = dataset_task2$alpha)

### check for normal distribution and homogenity of variances
par(mfrow = c(1, 1))
nv <- lm(time_adj ~ alpha * pipeline, data = dataset_task2)
ggqqplot(residuals(nv))
plot(nv, 1)

boxplot(dataset_task2$time_adj ~ dataset_task2$alpha)
anova(lm(time_adj ~ alpha * pipeline, data = dataset_task2))
anova_test(time_adj ~ alpha * pipeline,
           data = dataset_task2,
           effect.size = "pes")

### pairwise t-test
dataset_task2 %>%
  pairwise_t_test(time_adj ~ alpha,
                  pool.sd = F,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task2 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()


## separate t-tests
alpha_adj <- 0.05 * 0.05 * 0.05

t.test(dataset_task2[dataset_task2$alpha == '0', 'time_adj'],
       dataset_task2[dataset_task2$alpha == '1', 'time_adj'],
       var.equal = T)
t.test(dataset_task2[dataset_task2$alpha == '1', 'time_adj'],
       dataset_task2[dataset_task2$alpha == '5', 'time_adj'],
       var.equal = T)
t.test(dataset_task2[dataset_task2$alpha == '0', 'time_adj'],
       dataset_task2[dataset_task2$alpha == '5', 'time_adj'],
       var.equal = T)

##wilcox.test(value_t~alpha, data = dataset_task1[dataset_task1$alpha,])
describeBy(dataset_task2$f1.error, group = dataset_task2$alpha)

boxplot(dataset_task2$f1.error ~ dataset_task2$alpha)

# test statistic for f1 score difference
# followed by https://bjoernwalther.com/kruskal-wallis-test-in-r-rechnen/
k_test <- kruskal.test(dataset_task2$f1.error ~ dataset_task2$alpha)
print(k_test)
#kruskal.test(dataset_task1$f1_adj~dataset_task1$alpha)

# post hoc analysis of f1 score
#pairwise.wilcox.test(dataset_task1$f1.error,dataset_task1$alpha, p.adjust="bonferroni")
# bzw. für den Dunn's Test:
d <-
  dunn_test(f1.error ~ alpha, data = dataset_task2, p.adjust.method = "bonferroni")
print(d)

## F_value
#    Ab 0,1 ist es ein schwacher Effekt,
#    ab 0,25 ein mittlerer und
#    ab 0,4 ein starker Effekt.

### effect of dunn's test
d$r <- d$statistic / sqrt(d$n1 + d$n2)
d

#####################################################################
## H3: α = 1 is more accurate and faster than the alternatives for T3.
boxplot(dataset_task3$time ~ dataset_task3$alpha)
qqnorm(dataset_task3$time)
qqline(dataset_task3$time)
plot_dist(dataset_task3)
describeBy(dataset_task3$time, group = dataset_task3$alpha)
shapiro.test(dataset_task3$time)

b <- boxcox(lm(time ~ 1, data = dataset_task3))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task3$time_adj <- 1 / sqrt(dataset_task3$time / 10000)
plot_dist(dataset_task3, "time_adj", binwidth = 0.1)
shapiro.test(dataset_task3$time_adj)

qqnorm(dataset_task3$time_adj)
qqline(dataset_task3$time_adj)

leveneTest(lm(time_adj ~ alpha * pipeline, data = dataset_task3))
describeBy(dataset_task3$time_adj, group = dataset_task3$alpha)

### check for normal distribution and homogenity of variances
par(mfrow = c(1, 1))
nv <- lm(time_adj ~ alpha * pipeline, data = dataset_task3)
ggqqplot(residuals(nv))
plot(nv, 1)

boxplot(dataset_task3$time_adj ~ dataset_task3$alpha)
anova(lm(time_adj ~ alpha * pipeline, data = dataset_task3))
anova_test(time_adj ~ alpha * pipeline,
           data = dataset_task3,
           effect.size = "pes")

anova(lm(time_adj ~ alpha, data = dataset_task3[dataset_task3$pipeline ==
                                                  "o", ]))
anova(lm(time_adj ~ alpha, data = dataset_task3[dataset_task3$pipeline ==
                                                  "h", ]))

par(mfrow = c(2, 1))
boxplot(dataset_task3$time_adj ~ dataset_task3$alpha, data = dataset_task3[dataset_task3$pipeline ==
                                                                             "o", ])
boxplot(dataset_task3$time_adj ~ dataset_task3$alpha, data = dataset_task3[dataset_task3$pipeline ==
                                                                             "h", ])

t.test(dataset_task3[dataset_task3$alpha == '0', 'time_adj'],
       dataset_task3[dataset_task3$alpha == '1', 'time_adj'],
       var.equal = T)

t.test(dataset_task3[dataset_task3$alpha == '1', 'time_adj'],
       dataset_task3[dataset_task3$alpha == '5', 'time_adj'],
       var.equal = T)

t.test(dataset_task3[dataset_task3$alpha == '0', 'time_adj'],
       dataset_task3[dataset_task3$alpha == '5', 'time_adj'],
       var.equal = T)

### posthoc analysis for interaction effect
# from https://bjoernwalther.com/zweifaktorielle-anova-in-r-rechnen-und-interpretieren/
dataset_task3 %>%
  group_by(alpha) %>%
  anova_test(time_adj ~ pipeline, effect.size = "pes") %>%
  as.data.frame()

# bei signifikanter Interaktion (alpha*pipeline)
dataset_task3 %>%
  group_by(pipeline) %>%
  anova_test(time_adj ~ alpha, effect.size = "pes") %>%
  as.data.frame()


# Post-hoc-Analyse der ANOVAs aus A.

dataset_task3 %>%
  group_by(alpha) %>%
  emmeans_test(time_adj ~ pipeline, p.adjust.method = "bonferroni") %>%
  as.data.frame()


# Post-hoc-Analyse der ANOVAs aus B.

 dataset_task3 %>%
 group_by(pipeline) %>%
 emmeans_test(time_adj ~ alpha, p.adjust.method = "bonferroni") %>%
 as.data.frame()

#### prepare plots for post hoc

# Resultate speichern
ph <- dataset_task3 %>%
  group_by(alpha) %>%
  emmeans_test(time_adj ~ pipeline, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "alpha")

# Resultate in einen Boxplot übergeben
ggboxplot(dataset_task3,
          x = "alpha",
          y = "time_adj",
          color = "pipeline") +
  stat_pvalue_manual(ph) +
  labs(subtitle = get_test_label
       (
         anova_test(
           data = dataset_task3,
           time_adj ~ alpha * pipeline,
           effect.size = "pes"
         ),
         detailed = TRUE
       ))


#### Analyse the alpha differences for each pipeline separately
### check for normal distribution and homogenity of variances
par(mfrow = c(1, 1))
nv <-
  lm(time ~ alpha, data = dataset_task3[dataset_task3$pipeline == "o", ])
ggqqplot(residuals(nv))
plot(nv, 1)

par(mfrow = c(1, 1))
nv <-
  lm(time ~ alpha, data = dataset_task3[dataset_task3$pipeline == "h", ])
ggqqplot(residuals(nv))
plot(nv, 1)

plot_dist(dataset_task3[dataset_task3$pipeline == "h", ], "time_adj", binwidth = 0.05)
plot_dist(dataset_task3[dataset_task3$pipeline == "o", ], "time_adj", binwidth = 0.05)


anova_test(time_adj ~ alpha, data = dataset_task3[dataset_task3$pipeline == "h", ], effect.size =
             "pes")
anova_test(time_adj ~ alpha, data = dataset_task3[dataset_task3$pipeline == "o", ], effect.size =
             "pes")

### pairwise t-test
dataset_task3[dataset_task3$pipeline == "h", ] %>%
  pairwise_t_test(time_adj ~ alpha,
                  pool.sd = TRUE,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task3[dataset_task3$pipeline == "h", ] %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

dataset_task3[dataset_task3$pipeline == "o", ] %>%
  pairwise_t_test(time_adj ~ alpha,
                  pool.sd = TRUE,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task3[dataset_task3$pipeline == "o", ] %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()


##wilcox.test(value_t~alpha, data = dataset_task1[dataset_task1$alpha,])
describeBy(dataset_task3$f1.error, group = dataset_task3$alpha)

boxplot(dataset_task3$f1.error ~ dataset_task3$alpha)

# test statistic for f1 score difference
# followed by https://bjoernwalther.com/kruskal-wallis-test-in-r-rechnen/
k_test <- kruskal.test(dataset_task3$f1.error ~ dataset_task3$alpha)
print(k_test)
#kruskal.test(dataset_task1$f1_adj~dataset_task1$alpha)

# post hoc analysis of f1 score
#pairwise.wilcox.test(dataset_task1$f1.error,dataset_task1$alpha, p.adjust="bonferroni")
# bzw. für den Dunn's Test:
d <-
  dunn_test(f1.error ~ alpha, data = dataset_task3, p.adjust.method = "bonferroni")
print(d)

## F_value
#    Ab 0,1 ist es ein schwacher Effekt,
#    ab 0,25 ein mittlerer und
#    ab 0,4 ein starker Effekt.

### effect of dunn's test
d$r <- d$statistic / sqrt(d$n1 + d$n2)
d

##########################################################
## H4: α = 0 is more accurate and faster than the alternatives for T4
boxplot(dataset_task4$time ~ dataset_task4$alpha)
qqnorm(dataset_task4$time)
qqline(dataset_task4$time)
plot_dist(dataset_task4)
describeBy(dataset_task4$time, group = dataset_task4$alpha)
shapiro.test(dataset_task4$time)

b <- boxcox(lm(time ~ 1, data = dataset_task4))
lambda <- b$x[which.max(b$y)]
lambda

dataset_task4$time_adj <- sqrt(dataset_task4$time)
plot_dist(dataset_task4, "time_adj", binwidth = 10)
shapiro.test(dataset_task4$time_adj)

qqnorm(dataset_task4$time_adj)
qqline(dataset_task4$time_adj)

leveneTest(lm(time_adj ~ alpha * pipeline, data = dataset_task4))
describeBy(dataset_task4$time_adj, group = dataset_task4$alpha)

### check for normal distribution and homogenity of variances
par(mfrow = c(1, 1))
nv <- lm(time_adj ~ alpha * pipeline, data = dataset_task4)
ggqqplot(residuals(nv))
plot(nv, 1)

boxplot(dataset_task4$time_adj ~ dataset_task4$alpha)
anova(lm(time_adj ~ alpha * pipeline, data = dataset_task4))
anova_test(time_adj ~ alpha * pipeline,
           data = dataset_task4,
           effect.size = "pes")

### pairwise t-test
dataset_task4 %>%
  pairwise_t_test(time_adj ~ alpha,
                  pool.sd = F,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task4 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

0.05 * 0.05

t.test(dataset_task4[dataset_task4$alpha == '0', 'time_adj'],
       dataset_task4[dataset_task4$alpha == '1', 'time_adj'],
       var.equal = T)

t.test(dataset_task4[dataset_task4$alpha == '0', 'time_adj'],
       dataset_task4[dataset_task4$alpha == '5', 'time_adj'],
       var.equal = T)

##wilcox.test(value_t~alpha, data = dataset_task1[dataset_task1$alpha,])
describeBy(dataset_task4$f1.error, group = dataset_task4$alpha)

boxplot(dataset_task4$f1.error ~ dataset_task4$alpha)

# test statistic for f1 score difference
# followed by https://bjoernwalther.com/kruskal-wallis-test-in-r-rechnen/
k_test <- kruskal.test(dataset_task4$f1.error ~ dataset_task4$alpha)
print(k_test)
#kruskal.test(dataset_task1$f1_adj~dataset_task1$alpha)

# post hoc analysis of f1 score
#pairwise.wilcox.test(dataset_task1$f1.error,dataset_task1$alpha, p.adjust="bonferroni")
# bzw. für den Dunn's Test:
d <-
  dunn_test(f1.error ~ alpha, data = dataset_task4, p.adjust.method = "bonferroni")
print(d)

## F_value
#    Ab 0,1 ist es ein schwacher Effekt,
#    ab 0,25 ein mittlerer und
#    ab 0,4 ein starker Effekt.

### effect of dunn's test
d$r <- d$statistic / sqrt(d$n1 + d$n2)
d

##########################################################
## H5: α = 0.5 is neither the slowest nor least accurate for any task

### T1 - Time -> α = 0.5 is significant slower than α = 0 but equal to α = 1
### T1 - F1 error -> α = 0.5 has significant lower F1 score than α = 0 but equal to α = 1

### T2 - Time -> α = 0.5 is significant slower than α = 0 but equal to α = 1
### T2 - F1 error -> no difference in f1 error for any alpha combination

### T3 - Time -> α = 0.5 is between α = 1 and α = 0 in pipeline = h, slower than α = 1 and equal to α = 0 in pipeline = o
### T3 - F1 error -> α = 0.5 is more accurate than α = 0, but equal to α = 1

### T4 - time -> No difference
### T4 - F1 error -> No difference

##########################################################
## H6: Optimal UT drawings lead to faster answers for all tasks.

#### pairwise t-tests Task 1
dataset_task1 %>%
  group_by(alpha) %>%
  pairwise_t_test(time_adj ~ pipeline,
                  pool.sd = TRUE,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task1 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

ph <- dataset_task1 %>%
  group_by(alpha) %>%
  emmeans_test(time_adj ~ pipeline, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "alpha")

# Resultate in einen Boxplot übergeben
ggboxplot(dataset_task1,
          x = "alpha",
          y = "time_adj",
          color = "pipeline") +
  stat_pvalue_manual(ph) +
  labs(subtitle = get_test_label
       (
         anova_test(
           data = dataset_task1,
           time_adj ~ alpha * pipeline,
           effect.size = "pes"
         ),
         detailed = TRUE
       ))

#### pairwise t-tests Task 2
dataset_task2 %>%
  group_by(alpha) %>%
  pairwise_t_test(time_adj ~ pipeline,
                  pool.sd = TRUE,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task2 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

ph <- dataset_task2 %>%
  group_by(alpha) %>%
  emmeans_test(time_adj ~ pipeline, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "alpha")

# Resultate in einen Boxplot übergeben
ggboxplot(dataset_task2,
          x = "alpha",
          y = "time_adj",
          color = "pipeline") +
  stat_pvalue_manual(ph) +
  labs(subtitle = get_test_label
       (
         anova_test(
           data = dataset_task2,
           time_adj ~ alpha * pipeline,
           effect.size = "pes"
         ),
         detailed = TRUE
       ))

#### pairwise t-tests Task 3
dataset_task3 %>%
  group_by(alpha) %>%
  pairwise_t_test(time_adj ~ pipeline,
                  pool.sd = TRUE,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task3 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

ph <- dataset_task3 %>%
  group_by(alpha) %>%
  emmeans_test(time_adj ~ pipeline, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "alpha")

# Resultate in einen Boxplot übergeben
ggboxplot(dataset_task3,
          x = "alpha",
          y = "time_adj",
          color = "pipeline") +
  stat_pvalue_manual(ph) +
  labs(subtitle = get_test_label
       (
         anova_test(
           data = dataset_task3,
           time_adj ~ alpha * pipeline,
           effect.size = "pes"
         ),
         detailed = TRUE
       ))

#### pairwise t-tests Task 4
dataset_task4 %>%
  group_by(alpha) %>%
  pairwise_t_test(time_adj ~ pipeline,
                  pool.sd = TRUE,
                  p.adjust.method = "bonferroni") %>%
  as.data.frame()

dataset_task4 %>%
  cohens_d(time_adj ~ alpha) %>%
  as.data.frame()

ph <- dataset_task4 %>%
  group_by(alpha) %>%
  emmeans_test(time_adj ~ pipeline, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "alpha")

# Resultate in einen Boxplot übergeben
ggboxplot(dataset_task4,
          x = "alpha",
          y = "time_adj",
          color = "pipeline") +
  stat_pvalue_manual(ph) +
  labs(subtitle = get_test_label
       (
         anova_test(
           data = dataset_task4,
           time_adj ~ alpha * pipeline,
           effect.size = "pes"
         ),
         detailed = TRUE
       ))

############################################
############################################
############################################
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

anova(lm(value ~ alpha * pipeline, data = dataset_task2))

anova(lm(value ~ alpha, data = dataset_task2[dataset_task2$pipeline == "o", ]))
anova(lm(value ~ alpha, data = dataset_task2[dataset_task2$pipeline == "h", ]))

## α = 1 is more accurate and faster than the alternatives for T3.
boxplot(dataset_task3$value ~ dataset_task3$alpha)

qq_data <- subset(dataset_task3$value, dataset_task3$alpha == "1")

hist(qq_data, freq = F)
fit <- fitdistr(qq_data, "log-normal")$estimate
lines(dlnorm(0:max(qq_data), fit[1], fit[2]), lwd = 3)

qq_data <- dataset_task3$value
qq_data <- log(qq_data)
qqnorm(qq_data)
qqline(qq_data)
shapiro.test(qq_data)

anova(lm(value ~ alpha * pipeline, data = dataset_task3))


anova(lm(value ~ alpha, data = dataset_task3[dataset_task3$pipeline == "o", ]))
anova(lm(value ~ alpha, data = dataset_task3[dataset_task3$pipeline == "h", ]))

### check all data
qq_data <- dataset_t$value
qqnorm(qq_data)
qqline(qq_data)

shapiro.test(qq_data)



anova(lm(value ~ pipeline * task, data = dataset_t))

describeBy(dataset_task4$value, group = dataset_task4$pipeline)







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
plot_dist(dat, c = logvalue)

ggplot(dataset, aes(x = task, y = time, fill = task)) + geom_boxplot() +
  guides(fill = FALSE)


anova(lm(value ~ alpha * pipeline, data = dataset_t[dataset_t$task == "4", ]))

anova(lm(value ~ task * alpha * pipeline, data = dataset[dataset$metric == "t", ]))

anova(lm(value ~ pipeline * task * alpha, data = dataset[dataset$metric == "t", ]))

boxplot((value / 1000) ~ task, data =  dataset[dataset$metric == "t", ])

boxplot((value / 1000) ~ pipeline, data =  dataset[dataset$metric == "t", ])
boxplot((value / 1000) ~ alpha + task, data =  dataset[dataset$metric == "t", ])

boxplot(dataset$value ~ dataset$pipeline)

boxplot(dataset$value ~ dataset$task)

boxplot(dataset$value ~ dataset$alpha)

qq_data <- dataset


qq_data <- dataset[(dataset$metric == "t") & (dataset$task == "4") &
                     (dataset$pipeline == "o") &
                     (dataset$alpha == "5"), ]
qqnorm(qq_data$value)
qqline(qq_data, )

##############################################
## OUTLIER

### remove outlier
summary(dataset_task1$value)

out <-
  which(dataset_task1$value >= quantile(dataset_task1$value, 0.99))
out
plot_dist(dataset_task1[-out, ])
qqnorm(dataset_task1[-out, 'value'])
qqline(dataset_task1[-out, 'value'])

skewness(dataset_task1[-out, 'value'])
kurtosis(dataset_task1[-out, 'value'])
