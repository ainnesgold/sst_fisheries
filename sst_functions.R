#growth rate plot
#0.3 commonly used as growth rate
#r = 0 when 9 degree difference - nonlethal critical thermal maxima was 36 degrees
#try this function with different optimal temps (that will be calculated in the anomaly part of sst_cmip6.R)
quad_r <- function (x, a=0.3, b=0, c=-0.0037) {
  a + b*x + c*x^2
}

curve(quad_r, from = -10, to = 10, xlab="SST deviation (°C)", ylab = "")
title(ylab= "Intrinsic growth rate", line=2, cex.lab=1.2)








##Piecewise linear function for K
# Define the linear equation
linear_equation <- function(x) {
  -2.1408 * x + 160.253358
}

# Create a function to restrict y between 10 and 101.3
#10 g/m2 is maintained so matter how warm it gets (Darling et al. 2017)
#max carrying capacity set to 101.3 g/m2 (MacNeil et al. 2015)

restricted_y <- function(x) {
  y <- linear_equation(x)
  if (y < 10) {
    y <- 10
  } else if (y > 101.3) {
    y <- 101.3
  }
  return(y)
}

# Test the restricted_y function with a value of x
x_value <- 28
result <- restricted_y(x_value)
cat("For x =", x_value, ", y =", result, "\n")

curve(linear_equation, from = 27, to = 35, xlab = "Mean SST (°C)", ylab = "")
title(ylab= bquote('Carrying capacity'~(g/m^2)), line=2, cex.lab=1.2)


##logarithmic - not sure this one makes sense. might leave out.
#log_equation <- function(x) {
 # - log(x) + 101.3
#}
#curve(log_equation, from = 0, to = 50, xlab = "Mean SST (°C)", ylab = "")


##quadratic - same points used to establish linear
#max is "current" temp and max K
# 3 degree deviation (2040 temp) associated with 6% decline in max
quad_K <- function (x, a=101.3, b=0, c=-.7) {
  y = a + b*x + c*x^2
 # if (y < 10) {
  #  y <- 10
  #}
  return(y)
}

quad_K(-3)

curve(quad_K, from = -3, to = 3, xlab="SST deviation (°C)", ylab = "")

