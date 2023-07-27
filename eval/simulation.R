install.packages('pwr')
set.seed(1)

# some power analysis if we want to do t tests
pwr::pwr.t.test(
  d = 0.05 / 0.125,
  power = .8,
  sig.level = .05,
  type = 'two.sample',
  alternative = 'two.sided'
)

hist(truncnorm::rtruncnorm(
  n = 100,
  a = 0,
  b = 1,
  mean = 0.5,
  sd = 0.125
),
breaks = 11)$d * 0.125

# some simulated data and bayesian model if we want to go that route
# based on statistical rethinking by richard mcelreath: =
# https://catalogplus.tuwien.at/permalink/f/8j3js/UTW_alma2147901530003336

library('rethinking')

# SCIENTIFIC MODEL
# W -> EA <- U
# W -> SA <- U
# accuray (either) is influenced by weight and some unknown stuff
# later maybe introduce more variables like which design was shown - whatever we think may reasonably influence accuracy

# SIMULATION MODEL
# the accuracies are ranked choices compared to a reference element
# we model it via a latent variable that is linearly related to weight
# a = slope, b = intercept, sd_noise = SD of noise around the line
# the actual accuracy is then a discretized verion of that to the required amount of levels
sim_eacc <- function(weights, N, a, b, sd_noise) {
  U <- rnorm(length(weights), 0, sd_noise)
  ea <- (1-weights) * a + b + U
  ea <- pmax(pmin(ea, 1), 0)
  ea <- round(ea * (N - 1)) + 1
  return(ea)
}

sim_sacc <- function(weights, N, s, i, sd_noise) {
  U <- rnorm(length(weights), 0, sd_noise)
  sr <- (weights) * s + i + U
  sr <- pmax(pmin(sr, 1), 0)
  sr <- round(sr * (N)) + 1
  return(sr)
}

n_participants <- 25
n_questions <- 4
n_references <- 3
n_trials_per_participant <- n_questions*n_references # 4 questions x 3 reference elements

w <- sample(seq(0,1,0.1),replace = T, n_participants*n_trials_per_participant)
eacc <- sim_eacc(w, 19, 1, 0, .3)
sacc <- sim_sacc(w, 4, 1,0, .25)
df <- cbind(w, eacc, sacc, expand.grid(1:n_references, 1:n_questions, 1:n_participants))
colnames(df) <- c('w', 'ea', 'sa', 'ref', 'task', 'participant')

plot(df$sa ~ df$w, xlim=c(0,1), ylim=c(1,5))

plot(df$sa~df$w)
simplehist(df$sa)

ggplot2::ggplot(df) +
  ggplot2::geom_bin2d(ggplot2::aes(x=w, y=sa), bins=10) +
  ggplot2::scale_fill_gradient(low='#eeeeee', high='#000000') +
  ggplot2::theme_minimal()



# STATISTICAL MODEL
# EA = OrderedLogit(phi_i, alpha)
# SA = OrderedLogit(phi_j, alpha)
# note: not the same phi and alpha

# Model without predictor variables
statmodel.sa.empty <- map(
  alist(
    sa ~ dordlogit(phi, c(a1,a2,a3,a4)),
    phi <- 0,
    c(a1,a2,a3,a4) ~ dnorm(0,1)
  ),
  data=df,
  start=list(a1=-.25, a2=-.1, a3=.1, a4=.25)
)

# Model with predictor
statmodel.sa <- map(
  alist(
    sa ~ dordlogit(phi, c(a1,a2,a3,a4)),
    phi <- b*w,
    b ~ dnorm(0,1),
    c(a1,a2,a3,a4) ~ dnorm(0,1)
  ),
  data=df,
  start=list(a1=-.25, a2=-.1, a3=.1, a4=.25)
)

# Compare model with predictor to model without predictor
# If there's an effect of weight on accuracy, the predictor should be the better model
plot(compare(statmodel.sa.empty, statmodel.sa))

# Check sign of slope of weight if it matches with simulation
precis(statmodel.sa)
logistic(coef(statmodel.sa)[1:4])

# Draw from posterior one way or another and check if fit is reasonable

# Get samples of parameter posteriors
# Plot boundaries of accuracies
# No effect looks like a bunch of horizontal lines
# An effect has smooth slopey boundaries reminiscent of a CDF
post.samples <- extract.samples(statmodel.sa)

# TODO do the same plot based on data directly out of the simulation model
# see what it looks like, compare to posterior samples

head(post.samples)
W <- seq(0,1,.1)
plot( 1,1, type='n', xlab='weight', ylab='probability', xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,length(W)-1), yaxp=c(0,1,4) )
for (s in 1:100){
  p <- post.samples[s,]
  ak <- as.numeric(p[1:4])
  phi <- p$b*W
  pk <- pordlogit(1:4, a=ak, phi=phi)
  for (i in 1:4){
    lines(W, pk[,i], col=alpha(i,0.1))
  }
}


# Simulate observations given weights, plot histogram of accuracy
# get inspiration from here for actual plots https://will-ball.github.io/ordinalplots/
# No effect: Heatmap is same-ish gray everywhere (= same probability for each accuracy with all weights)
# Effect: The corners should be lighter than the center/diagonal as the probability of high/low accuracy in/decreases with weight
dfhist <- do.call(rbind, lapply(seq(0,1,.25), function(wi) {
  return(cbind(wi, as.vector(sim(statmodel.sa, data = list(w = wi)))))
}))
colnames(dfhist) <- c('w','sa')
dfhist <- as.data.frame(dfhist)
dfhist$w <- as.factor(dfhist$w)
dfhist$sa <- as.factor(dfhist$sa)
ggplot2::ggplot(dfhist) +
  ggplot2::geom_bin2d(aes(y=sa,x=w), bins = 5*5) +
  ggplot2::scale_fill_gradient(low='#eeeeee', high='#000000') +
  ggplot2::theme_minimal()
