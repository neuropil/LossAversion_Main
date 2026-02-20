negLLprospect_lambda1 <- function(parameters, choiceset, choices) {
  # Negative log likelihood function with lambda fixed to 1.
  
  eps = .Machine$double.eps;
  
  # Note: no need to source choice_probability, you can define a new one or pass lambda=1
  choiceP = choice_probability_lambda1(choiceset, parameters);
  
  choiceP[choiceP==1] = 1 - eps;
  choiceP[choiceP==0] = eps;
  
  nll <- -sum(choices * log(choiceP) + (1 - choices) * log(1 - choiceP));
  return(nll)
}