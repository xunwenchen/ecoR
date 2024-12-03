

# Exercise: chi-squared statistic calculation integrated with Yates' correction
cal_x2_full <- function(Oi, Ei_prop){
  
  # only do correction when df = 1
  k = length(Oi)# obtain number of observation
  df = k - 1 # calculate degree of freedom
  if(df >= 2){
    
  # use the original equation
  # Oi is the observed value
  # Ei is expected proportional value
  t <- sum(Oi) # total number of observation
  r <- Ei_prop/sum(Ei_prop) # ratio/probability
  exp <- r*t # expected value
  chi_squared <- sum((Oi - exp)^2 / exp)# Calculate chi-squared statistic
  
  # Print the result
  return(chi_squared)
  
  } else {
    
    # use Yates' correction
    t <- sum(Oi) # total number of observation
    r <- Ei_prop/sum(Ei_prop) # ratio/probability
    exp <- r*t # expected value
    chi_squared_c <- sum((abs(Oi - exp)-0.5)^2 / exp)# Calculate chi-squared statistic
    
    # Print the result
    return(chi_squared_c)
    }
}



# 练习：试试用data.xlsx (sheet = 2)做参数之间的相关性分析
ww_data <- read_xlsx('data/data.xlsx', sheet = 2)
View(ww_data)
colnames(ww_data)

str(ww_data)
ww_data <- ww_data[ , -1] # remove first column: date

ww.pearson <- cor(ww_data) # default method = "pearson"
round(ww.pearson, 2)
# Reorder the variables prior to plotting
ww.pearson.o <- order.single(ww.pearson)

# pairs() is a function to plot a matrix of bivariate scatter 
# plots.panelutils.R is a set of functions that add useful 
# features topairs().
pairs(ww_data[ ,ww.pearson.o], 
      lower.panel = panel.smooth, 
      upper.panel = panel.cor,
      diag.panel = panel.hist, 
      main = "Pearson Correlation Matrix")

