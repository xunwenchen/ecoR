# Introduction to R ----

# case sensitive

# new variable covers older ones

# Clear the environment
rm(list = ls())

# Clear the console
cat("\014")  # This is equivalent to pressing Ctrl+L

# Clear all plots
if (!is.null(dev.list())) dev.off()

# Garbage Collection - to free up memory
gc()


# Basic ----
# ~~ Simple calculations ----
1 + 2 #addition
99 - 28 #subtraction
67 * 5.64 + 99 #multiplication and addition
(3+5)/(2+3) #division 

# ~~ Assign values to variables ----
x <- 2 #assign value 2 to variable x
x 

x <- 3 #it rewrites on top of the old value
x

y <- x + 10
y

z <- x - 2
i <- exp(2) #function exponential

rm(i) #remove variable

# check how to use exp() function
?exp()

# ~~ Logical questions ----
x <- 8
y <- 3
x > y
x < y
x == y
x != y
x >= y
x <= y

# 1. Basic Data Types ----
# ~~ Numeric (number) ----
x
class(x)
is.integer(x)

# ~~ Integer (number with no decimal) ----
y = as.integer(3) 

class(y)
is.integer(y)

as.integer(3.14)
as.integer("5.27")

z <- 2L
class(z)

as.integer("coffee or tea") # it will return NA

as.integer(TRUE) # the numeric value of TRUE 

as.integer(FALSE) # the numeric value of FALSE 

# ~~ Characters (text) ----
a <- wukong # it will return an error
a <- "Wukong"
class(a)
a <- c('Wukong', 'Bajie', 'Wujing', 'Sanzang')
class(a)

# ~~ Complex ----
z = 1 + 2i     # create a complex number 
z              # print the value of z 
class(z)       # print the class name of z 

sqrt(-1)       # square root of -1 
sqrt(-1+0i)    # square root of -1+0i 
sqrt(as.complex(-1)) 

# ~~~~ little game with complex numbers ----
# plot Mandelbrot set
# Load required libraries
library(ggplot2)

# Define the Mandelbrot function
mandelbrot <- function(x, y, max_iter) {
  z <- complex(real = x, imaginary = y)
  c <- z
  for (i in 1:max_iter) {
    z <- z^2 + c
    if (Mod(z) > 2) return(i)
  }
  return(max_iter)
}


# Ask AI to explain the code

# Generate a grid of points
a <- seq(-2.5, 1.5, length.out = 1000)
b <- seq(-1.5, 1.5, length.out = 1000)
grid <- expand.grid(x = a, y = b)

# Apply the Mandelbrot function to each point
grid$iter <- mapply(mandelbrot, grid$x, grid$y, MoreArgs = list(max_iter = 100))

# Plot the Mandelbrot set
ggplot(grid, aes(x = x, y = y, fill = iter)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "white") +
  theme_minimal() +
  labs(title = "Mandelbrot Set", x = "Re", y = "Im")

# ~~~~ using prime numbers to plot spirals ----
# Load required libraries
library(ggplot2)
library(dplyr)

# Generate prime numbers
generate_primes <- function(n) {
  primes <- 2:n
  for (i in 2:sqrt(n)) {
    primes <- primes[primes == i | primes %% i != 0]
  }
  return(primes)
}

# Generate prime numbers up to 1000
primes <- generate_primes(999)

# Create a data frame for plotting
spiral_data <- data.frame(
  n = primes,
  angle = primes * pi / 180,
  radius = sqrt(primes)
)

# Calculate x and y coordinates
spiral_data <- spiral_data %>%
  mutate(x = radius * cos(angle),
         y = radius * sin(angle))

# Plot the spiral
ggplot(spiral_data, aes(x = x, y = y)) +
  geom_point(color = "indianred") +
  theme_minimal() +
  labs(title = "Prime Number Spiral", x = "X", y = "Y")




# 2. Vectors ----
v1 <- c(1, 2, 3)
v2 <- c(4, 5, 6, 7, 8)
v3 <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun') # not shown as a vector
v4 <- c(1, 2, 3, 'Mon', 'Tue', 'Wed') # not shown as a vector
v5 <- rep(1.5, 10) # repeat 1.5 ten times
v6 <- seq(1, 10, by = 2) # sequence from 1 to 10 by 2
v7 <- seq(1, 10, length = 5) # sequence from 1 to 10 with 5 elements
v8 <- seq(1, 10, length.out = 5) # sequence from 1 to 10 with 5 elements
v9 <- 1:100 # sequence from 1 to 100
v1 <- seq(from = 1, to = 100, by = 7)
v1[c(10, 12)] # select the 10th and 12th elements
v1[10, 12] # NOT working
v1[10:12] # select the 10th to 12th elements
v2 <- v1[-3] # remove the 3rd element
v3 <- sample(100:1000, 99) # randomly select 99 numbers from 100 to 1000
hist(v3) # plot a histogram of v3
v4 <- rnorm(50, 60, 10) # generate 50 random numbers with mean 60 and standard deviation 10
hist(v4) # plot a histogram of v4 
v4 <- rnorm(999, 100, 100) # generate 50 random numbers with mean 100 and standard deviation 100
hist(v4) # plot a histogram of the new v4

# summary
v1
length(v1) # length of v1 
summary(v1) # summary of v1
mean(v1) # mean of v1
sum(v1) # sum of v1
sd(v1) # standard deviation of v1
min(v1) # minimum of v1
max(v1) # maximum of v1
range(v1) # range of v1
quantile(v1) # quantiles of v1
median(v1) # median of v1
var(v1) # variance of v1
sort(v1) # sort v1
order(v1) # order of v1
rev(v1) # reverse v1

# 3. Matrices ----
# data structures in two dimensions, numbers only
length(b)
mt1 <- matrix(b, nrow = 100, ncol = 10) # create a matrix with 100 rows and 10 columns
mt1 
View(mt1) # view the matrix

# reaction a vector containing 100 names, from sample1 to sample100
sample_id <- paste0('sample_', 1:100)
rownames(mt1) <- sample_id # assign the names to the columns
mt1

# create a vector containing 10 soil properties
soil_prop <- c('pH', 'N', 'P', 'K', 'Ca', 'Mg', 'S', 'Fe', 'Mn', 'Zn')
colnames(mt1) <- soil_prop # assign the names to the rows
head(mt1)

hist(as.data.frame(mt1)$pH) # plot a histogram of pH

# nagvigate the matrix
mt1[39, 4]
mt1[39, 'K']
mt1[39, c('K', 'Ca')]
mt1[39, c(4, 5)]
mt1[39, ] # select the 39th row
mt1[, 4] # select the 4th column
mt1[, 'K'] # select the K column
mt1[, c('K', 'Ca')] # select the K and Ca columns
mt1[, c(4, 5)] # select the 4th and 5th columns
mt1[1:10, 1:5] # select the first 10 rows and first 5 columns

# 4. Data frame ----
# data that combines matrix with different data types
data() # check pre-installed data sets
df <- iris # load the iris data set
dim(df) # check the iris data set
head(df) # view the first 6 rows of the iris data set
class(df) # check the class of the iris

summary(df) # check the summary of the iris
colnames(df) # check the column names of the iris
rownames(df) # check the row names of the iris

df[, 1] # select the first column
df[4, ] # select the 4th row
df$Petal.Width # select the species column
hist(df$Petal.Width) # plot a histogram of the Petal.Width column


# 5. List ----
# a vector concatenating different data objects of diff. kinds
list1 <- list(v1, v2, v3, v4, mt1, df) # create a list
class(list1) # check the class of the list
length(list1) # check the length of the list
list1[1] # select the first element
list1[[1]] # select the first element
list1[[6]] # select the 6th element
list1[[6]]$Petal.Width # select the Petal.Width column in the 6th element
list1[[6]][4, ] # select the 4th row in the 6th element
list1[[6]][, 4] # select the 4th column in the 6th element

str(list1)

# 6. Exercises ----
# https://www.w3resource.com/r-programming-exercises/


# 1. Write an R Program for “Hello Geeks”.
print("Hello Geeks")

# 2. Find the Sum, Mean and Product of the Vector
v <- c(1, 2, 3, 4, 5)
sum(v)
mean(v)
prod(v) # product

# 3. Create a sequence of numbers from 200 to 1000 incrementing by 3
seq(200, 1000, by = 3)

# 4. sum of numbers from 51 to 91
v <- c(51:91)
sum(51:91)
sum(v)

# 3. Create an R Program to Take Input From the User
# 4. How to Generate Random Numbers from Standard Distributions in R
# 5. R Program to Sample from a Population
# 6. Create an R Program to Find the Minimum and Maximum
# 7. R Program to Sort a Vector
# 8. How to Find the Factorial of a Number
# 9. How to create R Multiplication Table
# 10. Write an R Program to Check Prime Number
# 11. R Program to check Armstrong Number
# 12. R Program to Print the Fibonacci Sequence
# 13. R Program to Check for Leap Year
# 14. Check if a Number is Odd or Even in R Programming

# Using R ----
# 1. Import data ----
getwd()

# ways to install packages
# ~~ search within RStudio
# ~~ comment line 

# install package 'readxl'
# install.packages('readxl')

# load the package
library(readxl)


# Yan L, Hermans SM, Totsche KU, Lehmann R, Herrmann M, Küsel K (2021) Groundwater bacterial communities evolve over time in response to recharge. Water Res 201:117290. https://doi.org/10.1016/j.watres.2021.117290

# Wang Y, Ye J, Ju F, Liu L, Boyd JA, Deng Y, Parks DH, Jiang X, Yin X, Woodcroft BJ, Tyson GW, Hugenholtz P, Polz MF, Zhang T (2021) Successional dynamics and alternative stable states in a saline activated sludge microbial community over 9 years. Microbiome 9:199. https://doi.org/10.1186/s40168-021-01151-5


# import data in excel format and take the second sheet
df <- read_excel('data/data.xlsx', sheet = 2)
# df <- readxl::read_excel('data/data.xlsx', sheet = 2)
# df <- read.csv('data/data.csv')

# Note: data structure before importing? --> PPT slide ~33

# view data
View(df)


# find something in a file, and replace
# Ctrl + f

# find something in all files, and replaces
# ctrl + shift + f


# add a column time_step as the first column and the values are starting from 1 to the number of rows
df <- data.frame(time_step = 1:nrow(df), df)

dim(df)

View(df)

unique(df$pH)

colnames(df)

df$Temperature

hist(df$Temperature)

summary(df)
skimr::skim(df) # call skimr package without charing it into the Environment

# plot line chart of df$time_step and df$Rainfall
plot(df$Date, df$Rainfall, type = 'l', col = 'indianred', xlab = 'Time Step', ylab = 'Rainfall')


# use ggplot2 to plot the line chart with data points
library(ggplot2)
plot_rf <- ggplot(df, aes(x = time_step, y = Rainfall)) +
  geom_line(color = 'indianred') +
  labs(title = 'Rainfall over Time', x = 'Time Step', y = 'Rainfall')+
  # theme black and white without grid
  theme_bw() +
  theme(panel.grid = element_blank())+
  # add points
  geom_point(color = 'grey', alpha = 0.5, size = 2)
plot_rf
save_plot <- TRUE
if(save_plot){
  ggsave('out/rainfall.jpg', plot = plot_rf, width = 10, height = 6)
}

df_temp_op <- subset(df, df$Temperature > 25) # subset df when temp > 25

hist(df_temp_op$Temperature)

# 2. Data manipulation ----
head(df)
head(df, 10)
tail(df)
tail(df, 10)

colnames(df)
unique(df$DO)
length(unique(df$DO))
dim(df)
df[1, ]
df[ , 4]
df[ , 'DO']
df[ , c('DO', 'pH')]
df$Date
hist(df$DO)


sort(df$Rainfall)
cut_off <- 100
dry <- df[df$Rainfall < cut_off, ]
dim(dry)
dim(df)
hist(dry$Rainfall)


# create a new column with $
data()
npk <- npk

npk$all_fert<-paste(npk$N,npk$P,npk$K,sep='_')

# Exercise: Subset the data by N, just take the N fertilized.
fertN<- subset(npk, npk$N == 1)


# paste df with vectors
rbind(); cbind()

all_fert2 <- paste(npk$N, npk$P, npk$K, sep = '_')

npk_new <- cbind(npk, all_fert2)

# order a table
order(npk_new$yield)
sort(npk_new$yield)
sort(npk_new$yield, decreasing = TRUE)

# 3. Create basic models and statistics ----

# linear model
m1<- lm(yield ~ N +P+K, data = npk) # apply log(N)
m1

hist(resid(m1), nclass = 10)
summary(m1)

boxplot(yield~N, data = npk)

sum <- summary(m1)

sum$coefficients

# generalized linear model (glm)
m2 <- glm(yield ~ N+P+K, data = npk, family = Gamma (link = 'inverse')) # apply log(N)

summary(m2)

# Exercise: are there any differences in yield due to its block?

m3<-lm(yield~block,data=npk)
summary(m3)

# the answer is no. 

# linear mixed models - Back to course PPT ~ slide 37 about experimental designs of two contexts.

library(nlme)
library(lmerTest)

mm1 <- lme(yield ~ N+P+K, random = ~1 | block, data = npk, na.action = na.omit) # add na.omit

hist(resid(mm1), nclass = 8)

summary(mm1)

a <- summary(mm1)

a$varFix

save(mm1, 'out/mm1.R')

write.csv(mm1, 'out/mm1.csv')

# explanation: 
# fits a linear mixed-effects model to the npk dataset, with yield as the response variable and N, P, and K as predictor variables. 
# includes block as a random effect. 
# The na.omit argument ensures that any rows with missing values are omitted.  plots a histogram of the residuals to check their distribution.


# Hacking mode getting more interesting ----
# 1. Loops ----
# for
time <- 1:30
for(i in time){
  print(paste('The sampling time is', i))
}

# while
i <- 1
while (i < 10){
  print(i)
  i = i+1
}

# if else
x <- 10
if(x > 2){
  print('ok')
} else{
  print('not ok')
}

# 2. Functions ----
square <- function(x){

# square the argument and assign to y
y <- x^2 

# return the value
return(y)

}

square(9)



# multiple arguments and default values 
center1<-function(x, fun="mean", p=0.5 ){

if(fun == "mean"){x0 <- x-mean(x)}
else if(fun == "median"){x0 <- x-median(x)}
else if(fun == "quantile"){x0 <- x-quantile(x, probs=p)} # meaning of p is the quantile
else {stop("unrecognized function name")}
return(x0)

}

# define a vector for testing the function
y <- (1:9)^2
mean(y)
center1()
center1(y)# default is mean
center1(y, "median") # median
center1(y, "quantile") # default p = 0.5, return x - 50th quantile = median
center1(y, "quantile", 0.25) # return x - 25th quantile
center1(y, "quantile", 0.75) # return x- 75th quantile

# ~~ how many chicken and rabbits problem ----

how_many <- function(head.no, feet.no)
{
  chick.no = (4*head.no - feet.no)/2
  rabbit.no = head.no - chick.no
  return(list(Chick_number = chick.no, Rabbit_number = rabbit.no))
  # ans <- list(chick.no, rabbit.no)
  # return(ans)
  # 
}

# if the number of heads is 35 and the number of feet is 94

how_many(35, 94)

# if numbers of heads and feet are not ok, like 3 heads and more than 12 feet
how_many(3, 13) # the function will return non-integers 



# ~~ percentage increase / decrease calculator ----
# Back to course PPT ~ slide 37
perc_increase <- function(initial_val, final_val){
  diff = final_val-initial_val
  change = diff/abs(initial_val)
  return(change)
}

perc_decrease <- function(initial_val, final_val){
  diff = initial_val-final_val
  change = diff/abs(initial_val)
  return(change)
}

perc_increase(25, 42)
perc_decrease(42, 25)

# ~~ source fun.R ----
# write functions in fun.R first 
source('code/fun.R')

# data input
Oi <- c(23, 45, 34, 56, 37)
Ei_prop <- c(3, 5, 4, 5, 4)

x2 <- cal_x2(Oi, Ei_prop)
x2
# checking using existing function in another package
# Calculate expected probabilities
p <- Ei_prop/sum(Ei_prop)

# Perform chi-squared test using chisq.test()
x2_t <- chisq.test(Oi, p = p)$statistic
x2_t

if (x2 == x2_t) {
  print('Congratulations! Your manual calculation is correct.')
} else {
  print('Please check your calculation again. No worries, it is almost there.')
}

# Exercise: add Yates' correction in the cal_x2() function. 


# Compare the results
identical(data2$statistic, x21)
# 3. Good manners and tricks 好习惯和小技巧 ----

# format according to the subprocess 
# same policy with spaces
# search tool
# cmd/ctrl + shift + c <- convert all selected text to annotation and vice versa
# click on a key. It will highlight in bold the closing one
# you hide whole section by just clicking on the arrow next to its title
# gc() for cleaning memory

# 4. Practical examples ----

library(nlme)
npk <- npk
elem = colnames(npk[2:4])
data = npk
i = 1
mixed_m <- function(elem,data){
  out <- data.frame() # create an empty data frame
  for(i in 1:length(elem)){
    form <- paste("lme(yield ~",elem[i],", random = ~1 | block, data = data, na.action = na.omit)") # create model expression
    m4 <-eval(parse(text = form)) # run model expression
    sum<-summary(m4) # save model summary
    sum_mods <- data.frame(sum$tTable, sum$AIC, elem[i])# select and format ir
    out <- rbind(out, sum_mods) # add the results in the object "out" previously created.
    }
  colnames(out)[6:7] <- c('AIC', 'element') # modify some column names
  return(out) # define the output
}

treatment <-mixed_m(elem = colnames(npk[2:4]), data = npk)
rownames(treatment) <- c("Control N", "N fert", "Control P", "P fert", "Control K", "K fert")
treatment

options(scipen = 999) # sci notation to decimal style
write.csv(treatment, file = 'out/trt_example.csv')


  


