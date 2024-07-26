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


## Accuracy Plot
dataset_a <- read.csv('survey_long.csv')
dataset_a <- dataset_a[dataset_a$metric == 'a', ]
dataset_a[dataset_a$alpha == 5, ]$alpha <- 0.5
dataset_a[dataset_a$task == 1,]$task <- 'T1'
dataset_a[dataset_a$task == 2,]$task <- 'T2'
dataset_a[dataset_a$task == 3,]$task <- 'T3'
dataset_a[dataset_a$task == 4,]$task <- 'T4'

ggplot2::ggplot(dataset_a) +
  ggplot2::aes(
    x = factor(alpha),
    y = value,
    fill = factor(alpha),
    group = task
  ) +
  ggplot2::scale_y_continuous(breaks=seq(0,1,.2))+
  ggdist::stat_histinterval(
    orientation='vertical',
    breaks = seq(0, 1, 0.1),
    point_interval = 'median_qi',
    slab_alpha=.5,
    .width = c(.5, .75),
    position = 'dodgejust',
    normalize = 'panels',
    side='both'
  ) +
  ggplot2::scale_fill_discrete() +
  ggplot2::theme_bw()+
  ggplot2::theme(
    legend.position='none',
    axis.text = ggplot2::element_text(size=14),
    axis.title.x = ggplot2::element_text(size=18),
    axis.title.y = ggplot2::element_text(size=18),
    strip.text = ggplot2::element_text(size = 18)
  )+
  ylab('F1 Score')+
  xlab('Weight')+
  ggplot2::facet_grid(cols=vars(task))

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
  ggdist::stat_eye(
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
