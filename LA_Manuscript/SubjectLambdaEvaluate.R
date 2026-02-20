library(readr)
library(tibble)
library(dplyr)

claseFILETEST <- read_csv("Github/LossAversionManuscript/claseFILETEST.csv")
View(claseFILETEST)

choices <- tibble(
  participant = claseFILETEST$subjID,
  gain = claseFILETEST$riskygain,
  loss = claseFILETEST$riskyloss,
  choice = claseFILETEST$choice
)


choice_probability <- function(choiceset, parameters) {
  # A function to calculate the probability of taking a risky option
  # using a prospect theory model.
  # Assumes parameters are [rho, lambda, mu] as used in S-H 2009, 2013, 2015, etc.
  # Assumes choiceset has columns riskygain, riskyloss, and certainalternative.
  # Creates binary choices, with True/1 = risky, False/0 = safe.
  #
  # Peter Sokol-Hessner
  # July 2021
  
  # extract  parameters
  rho = as.double(parameters[1]); # risk attitudes
  lambda = as.double(parameters[2]); # loss aversion
  mu = as.double(parameters[3]); # choice consistency
  
  # calculate utility of the two options, making as few assumptions as possible about the +/- sign of the values
  utility_risky_1 = (choiceset$riskygain >= 0) * abs(choiceset$riskygain)^rho + 
    (choiceset$riskygain < 0) * -lambda * abs(choiceset$riskygain)^rho;
  
  utility_risky_2 = (choiceset$riskyloss >= 0) * abs(choiceset$riskyloss)^rho + 
    (choiceset$riskyloss < 0) * -lambda * abs(choiceset$riskyloss)^rho;
  
  utility_risky_option = 0.5 * utility_risky_1 + 0.5 * utility_risky_2
  
  utility_safe_option = (choiceset$certainalternative >= 0) * abs(choiceset$certainalternative)^rho + 
    (choiceset$certainalternative < 0) * -lambda * abs(choiceset$certainalternative)^rho;
  
  # normalize values using this term
  div <- max(choiceset)^rho; # decorrelates rho & mu
  
  # calculate the probability of selecting the risky option
  p = 1/(1+exp(-mu/div*(utility_risky_option - utility_safe_option)));
  
  return(p)
}

claseFILETEST <- dplyr::filter(
  claseFILETEST,
  !is.na(certainalternative),
  !is.na(riskygain),
  !is.na(riskyloss)
)

claseFILE2 = claseFILETEST[,c('certainalternative','riskygain',
                             'riskyloss')]

start1 <- c(lambda=1.2, rho=0.8, mu=1)
choiC = choice_probability(claseFILE2, start1)


# Negative log-likelihood function for Model 1 (lambda estimated)
nll_model1 <- function(params, gains, losses, choices) {
  lambda <- params[1]
  rho <- params[2]
  mu <- params[3]
  
  # 1) compute utilities (example)
  U_gain <- gains^lambda
  U_loss <- -mu * (abs(losses)^lambda)   # now 0^λ = 0, negative^λ never occurs
  V       <- U_gain + rho * U_loss
  
  p_raw <- 1 / (1 + exp(-V))
  eps   <- .Machine$double.eps
  p     <- pmin(pmax(p_raw, eps), 1 - eps)
  
  if (any(is.na(choices))) {
    stop("Missing values in choices!")
  }
  
  ll <- ifelse(choices == 1, log(p), log(1 - p))
  
  return(-sum(ll))
}


nll_model2 <- function(params, gains, losses, choices) {
  # fixed exponent
  lambda <- 1
  
  rho <- params[1]      # loss weight
  mu  <- params[2]      # decision "temperature"
  
  U_gain <- gains^lambda
  U_loss <- -rho * (abs(losses)^lambda)   # flip‐sign & abs‐before‐power
  
  V <- U_gain + U_loss
  
  p_raw <- 1 / (1 + exp(-mu * V))
  
  eps <- .Machine$double.eps
  p_accept<- pmin(pmax(p_raw, eps), 1 - eps)
  
  if (any(is.na(choices))) {
    stop("Missing values in choices!")
  }
  
  nll <- -sum(ifelse(choices == 1,
                           log(p_accept),
                           log(1 - p_accept)))
  
  return(nll)
}



# Fit both models per participant
results <- choices %>%
  group_by(participant) %>%
  group_modify(~{
    df <- .x %>% 
      filter(!is.na(choice), !is.na(gain), !is.na(loss))
    
    # Model 1 starting values
    start1 <- c(lambda=1.2, rho=0.8, mu=1)
    # Model 2 starting values
    start2 <- c(rho=0.8, mu=1)
    
    # Model 1 fit
    fit1 <- optim(
      par = start1,
      fn = nll_model1,
      gains = df$gain,
      losses = df$loss,
      choices = df$choice,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01, 0.01),
      upper = c(5, 5, 10)
    )
    
    # Model 2 fit
    fit2 <- optim(
      par = start2,
      fn = nll_model2,
      gains = df$gain,
      losses = df$loss,
      choices = df$choice,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01),
      upper = c(5, 10)
    )
    
    # Likelihood ratio test
    D = -2 * (fit2$value - fit1$value)
    p = pchisq(D, df=1, lower.tail=FALSE)
    
    tibble(
      ll_model1 = -fit1$value,
      ll_model2 = -fit2$value,
      lambda_est = fit1$par[1],
      rho1 = fit1$par[2],
      mu1 = fit1$par[3],
      rho2 = fit2$par[1],
      mu2 = fit2$par[2],
      LR_D = D,
      LR_p = p
    )
  })

# View results
print(results)


# Negative log-likelihood function for Model 2 (lambda fixed to 1)
nll_model2 <- function(params, gains, losses, choices) {
  lambda <- 1
  rho <- params[1]
  mu <- params[2]
  
  util_gain  <- compute_utility(gains, lambda, rho)
  util_loss  <- compute_utility(losses, lambda, rho)
  
  ev = util_gain + util_loss
  
  p_accept = 1 / (1 + exp(-mu * ev))
  
  p_accept = pmin(pmax(p_accept, 1e-10), 1 - 1e-10)
  
  -sum(choices * log(p_accept) + (1 - choices) * log(1 - p_accept))
}
