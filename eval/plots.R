library(pwr)
library(MASS)
library(ggplot2)
library(car)
library(psych)
library(plyr)
library(moments)
library(rstatix)
library(emmeans)
library(ggpubr)
library(ggdist)


setwd('/Users/npiccolotto/Projects/cvast/bssvis/ensemble-set-rendering/eval/')

dataset <- read.csv('survey_long.csv')
dataset[dataset$metric == 't', ]$value <-
  dataset[dataset$metric == 't', ]$value / 1000
dataset <- dataset[dataset$metric == 't', ]
dataset[dataset$alpha == 5, ]$alpha <- 0.5

dataset_task1 <- subset(dataset, dataset$task == "1")
dataset_task2 <- subset(dataset, dataset$task == "2")
dataset_task3 <- subset(dataset, dataset$task == "3")
dataset_task4 <- subset(dataset, dataset$task == "4")

# Transform/adjust just as in stat-analysis.R
dataset[dataset$task == 1, ]$value <-
  sqrt(dataset[dataset$task == 1, ]$value)
dataset[dataset$task == 2, ]$value <-
  1 / sqrt(dataset[dataset$task == 2, ]$value)
dataset[dataset$task == 3, ]$value <-
  1 / sqrt(dataset[dataset$task == 3, ]$value)
dataset[dataset$task == 4, ]$value <-
  sqrt(dataset[dataset$task == 4, ]$value)
dataset[dataset$task == 1,]$task <- 'T1'
dataset[dataset$task == 2,]$task <- 'T2'
dataset[dataset$task == 3,]$task <- 'T3'
dataset[dataset$task == 4,]$task <- 'T4'

labeller_sqrt <- function(breaks) {
  return(purrr::map(breaks, function(b) {
    return(b*b)
  }))
}
labeller_sqrt_inv <- function(breaks) {
  return(purrr::map(breaks, function(b) {
    return(round((1/b)^2,2))
  }))
}

labeller_sqrt(c(0,1,2,3))
labeller_sqrt_inv(1:4)

## Plots using transformed time
ggplot(dataset[dataset$task == 'T1' | dataset$task == 'T4',]) +
  aes(
    y = factor(alpha),
    x = value,
    fill = factor(alpha),
    group = factor(task)
  ) +
  stat_eye(
    point_interval = 'mean_hdci',
    slab_alpha=.5,
    .width = c(.5, .75),
    position = 'dodgejust'
  ) +
  theme_bw()+
  theme(
    legend.position='none',
    axis.text = element_text(size=14),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    strip.text = element_text(size = 18)
  )+
  ylab('Weight')+
  xlab('Seconds')+
  scale_x_continuous(labels=labeller_sqrt)+
  scale_fill_discrete() +
  facet_grid(scales = 'free_x', rows = 'task')

ggplot(dataset[dataset$task == 'T2' | dataset$task == 'T3',]) +
  aes(
    y = factor(alpha),
    x = value,
    fill = factor(alpha),
    group = factor(task)
  ) +
  stat_eye(
    point_interval = 'mean_hdci',
    slab_alpha=.5,
    .width = c(.5, .75),
    position = 'dodgejust'
  ) +
  theme_bw()+
  theme(
    legend.position='none',
    axis.text = element_text(size=14),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    strip.text = element_text(size = 18)
  )+
  ylab('Weight')+
  xlab('Seconds')+
  scale_x_reverse(labels=labeller_sqrt_inv)+
  scale_fill_discrete() +
  facet_grid(scales = 'free_x', rows = 'task')

## Accuracy Plot
dataset_a <- read.csv('survey_long.csv')
dataset_a <- dataset_a[dataset_a$metric == 'a', ]
dataset_a[dataset_a$alpha == 5, ]$alpha <- 0.5
dataset_a[dataset_a$task == 1,]$task <- 'T1'
dataset_a[dataset_a$task == 2,]$task <- 'T2'
dataset_a[dataset_a$task == 3,]$task <- 'T3'
dataset_a[dataset_a$task == 4,]$task <- 'T4'

ggplot(dataset_a) +
  aes(
    x = factor(alpha),
    y = value,
    fill = factor(alpha),
    group = task
  ) +
  scale_y_continuous(breaks=seq(0,1,.25))+
  stat_histinterval(
    orientation='vertical',
    breaks = seq(0, 1, 0.1),
    point_interval = 'median_qi',
    slab_alpha=.5,
    .width = c(.5, .75),
    position = 'dodgejust',
    normalize = 'panels',
    side='both'
  ) +
  scale_fill_discrete() +
  theme_bw()+
  theme(
    legend.position='none',
    axis.text = element_text(size=14),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    strip.text = element_text(size = 18)
  )+
  ylab('F1 Score')+
  xlab('Weight')+
  facet_grid(cols=vars(task))

## Plot using measured time
dataset <- read.csv('survey_long.csv')
dataset[dataset$metric == 't', ]$value <-
  dataset[dataset$metric == 't', ]$value / 1000
dataset <- dataset[dataset$metric == 't', ]
dataset[dataset$alpha == 5, ]$alpha <- 0.5
dataset[dataset$task == 1,]$task <- 'T1'
dataset[dataset$task == 2,]$task <- 'T2'
dataset[dataset$task == 3,]$task <- 'T3'
dataset[dataset$task == 4,]$task <- 'T4'

ggplot(dataset) +
  aes(
    x = factor(alpha),
    y = value,
    fill = factor(alpha),
    group = factor(task)
  ) +
  stat_eye(
    point_interval = 'mean_hdci',
    slab_alpha=.5,
    normalize='panels',
    .width = c(.5, .75),
    position = 'dodgejust'
  ) +
  theme_bw()+
  theme(
    legend.position='none',
    axis.text = element_text(size=14),
    axis.title.x = element_text(size=18),
    axis.title.y = element_text(size=18),
    strip.text = element_text(size = 18)
  )+
  xlab('Weight')+
  ylab('Seconds')+
  scale_y_continuous(limits = c(0,90),breaks=c(0,10,20,30,45,60,90))+
  scale_fill_discrete() +
  facet_grid(scales = 'free_x', cols = vars(task))
