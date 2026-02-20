choice_probability_lambda1 <- function(choiceset, parameters) {
  # Same as choice_probability, but lambda = 1
  
  rho = as.double(parameters[1])
  lambda = 1
  mu = as.double(parameters[2])
  
  utility_risky_1 = (choiceset$riskygain >= 0) * abs(choiceset$riskygain)^rho + 
    (choiceset$riskygain < 0) * -lambda * abs(choiceset$riskygain)^rho;
  
  utility_risky_2 = (choiceset$riskyloss >= 0) * abs(choiceset$riskyloss)^rho + 
    (choiceset$riskyloss < 0) * -lambda * abs(choiceset$riskyloss)^rho;
  
  utility_risky_option = 0.5 * utility_risky_1 + 0.5 * utility_risky_2
  
  utility_safe_option = (choiceset$certainalternative >= 0) * abs(choiceset$certainalternative)^rho + 
    (choiceset$certainalternative < 0) * -lambda * abs(choiceset$certainalternative)^rho;
  
  div = max(choiceset)^rho;
  
  p = 1 / (1 + exp(-mu / div * (utility_risky_option - utility_safe_option)));
  
  return(p)
}