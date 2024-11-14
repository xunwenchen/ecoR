# Build own functions here ----

# calculate chi-squared statistic
# Oi is observed vector, and Ei is expected proportional value

cal_x2 <- function(Oi, Ei_prop){
  # Oi is the observed value
  # Ei is expected proportional value
  t <- sum(Oi) # total number of observation
  r <- Ei_prop/sum(Ei_prop) # ratio/probability
  exp <- r*t # expected value
  chi_squared <- sum((Oi - exp)^2 / exp)# Calculate chi-squared statistic
  
  # Print the result
  return(chi_squared)
}
