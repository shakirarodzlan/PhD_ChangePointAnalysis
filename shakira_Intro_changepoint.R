
## Introduction to change point analysis (Youtube, Dr Rebecca Killick)
    
      ## https://www.youtube.com/watch?v=UfGrLJ7S3sc


# install & load package

install.packages("changepoint")
library(changepoint)

# set working directory

setwd("C:/Users/shakirarodzlan/Desktop/PhD/PhD Chapters/Chapter 4-Joint point regression/Intro change point using R/Dr Rabecca Intro Change Point _ youtube")


#...................................................................

#1) Single change in mean (cpt.mean), example code

## using cpt.mean function

# cpt.mean(data, penalty = "MBIC", pen.value = 0, method = "AMOC",
        # Q = 5, test.stat = "normal", class = TRUE, param.estimates = TRUE,
        # minseglen = 1)


## scale the data if variance is not 1 (because cpt.mean assumed that the variance of a time series is 1)

set.seed(1)
m1 = c (rnorm(100, 0, 1), rnorm(100, 5, 1))
m1.amos = cpt.mean(m1)
cpts(m1.amos)


## plot

plot(m1.amos)


# Hands on other data - download this data from github Dr Rebecca
  ## (https://github.com/rkillick/intro-changepoint-course/blob/master/GPvisitsWeekNov1718.Rdata)

  
load("GPvisitsWeekNov1718.Rdata")   #GP appointment data
ts.plot(GPvisitsWeekNov1718)

## to get change point

GP.default = cpt.mean(GPvisitsWeekNov1718)
cpts(GP.default)

param.est(GP.default)

## plot

plot(GP.default)

## Scale the data - because variance is not equal to 1 
   # note that cpt.mean assumed that variance og time series is 1.

GP.scale = cpt.mean(as.vector(scale(GPvisitsWeekNov1718)))
cpts(GP.scale)

plot(GP.scale)


#.................................................

# 2) Multiple changepoint

## using cpt.var function

set.seed(1)
v1 = c(rnorm(100, 0, 1), rnorm(100, 0, 2), rnorm(100, 0, 10),
       rnorm(100,0,9))
v1.man = cpt.var(v1, method = "PELT", penalty = "Manual", 
                 pen.value = "2*log(n)")
cpts(v1.man)
param.est(v1.man)

plot(v1.man, cpt.width=3)


# using cpt.meanvar function

set.seed(1)
mv1 = c(rexp(50, rate = 1), rexp(50,5), rexp (50,2),
        rexp (50,7))
mv1.pelt = cpt.meanvar(mv1, test.stat = "Exponential",
                      method = "BinSeg", Q=10, penalty = "SIC")
cpts(mv1.pelt)

param.est(mv1.pelt)

plot (mv1.pelt, cpt.width = 3, cpt.col = "blue")


### Try other example data

## https://github.com/rkillick/changepoint/blob/master/data/HC1.RData

  ##load HC1 data 

data(HC1)
ts.plot(HC1)

hc1.pelt = cpt.meanvar(HC1, method = "PELT", penalty = "Manual",
                       pen.value = 14)
ncpts(hc1.pelt)

plot (hc1.pelt, ylab= "G+C content", cpt.width = 3)

  ## using CROP function (change for range of penalties)

v1.crops = cpt.var (v1, method = "PELT",
                    penalty = "CROPS", pen.value = c(5, 500))

cpts.full(v1.crops)

pen.value.full(v1.crops) # to get specific penalty values

plot(v1.crops, ncpts= 5)

plot(v1.crops, diagnostic = TRUE)



#........................................................................................


# Non-paramatric test using cpt.np function

install.packages("changepoint.np")
library(changepoint.np)

set.seed(12)
J <- function(x) {(1 + sign(x))/2} #function
n <- 1000 # number of data/observation
tau <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 
         0.65, 0.76, 0.78, 0.81) * n        # this is change point
h <- c(2.01, -2.51, 1.51, -2.01, 2.51, -2.11, 
       1.05, 2.16, -1.56, 2.56, -2.11)      # this is parameter to control
sigma <- 0.5
t <- seq(0, 1, length.out =n)
 # to create the data
data <- array()                            
for (i in 1: n){
  data [i] <- sum (h*J(n*t [i] - tau)) + (sigma * rnorm (1))
  }

out <- cpt.np (data, method="PELT", minseglen= 2,
               nquantiles = 4* log(length(data)))
cpts(out)
plot(out)


## other example for non parametric
    ## used data HeartRate from the changepoint.np package.
    ## use one of non-parametric functions to see if there is evidence for change in heart rate

data ("HeartRate")
HR.pelt = cpt.np (HeartRate, method = "PELT",
                  nquantiles = 4 * log(length(HeartRate)))               
ncpts(HR.pelt)
plot(HR.pelt)
# use CROP
HR.crops = cpt.np (HeartRate, penalty ="CROPS", 
                   pen.value = c(5, 200), method = "PELT",
                    minseglen = 2, nquantiles = 4 * log(length(HeartRate))) 
## diagnostic plot
plot(HR.crops, diagnostic = TRUE)
abline (v= 11, col = "red")  
abline (v= 15, col= "red")

plot(HR.crops, ncpts= 11) # u can choose no of changepoint
plot(HR.crops, ncpts= 15)



#.................................................................................


# Checking assumption 
  ## The main assumption for normal likelihood ratio test for change in mean are
     ### Independent data point
     ### Normal distribution point pre and post
     ### Constant variance across the data

ts.plot(m1)

 ## check normal distribution

  #Histogram
hist(m1)
  # shapiro test
shapiro.test(m1)
  # kolmogorov test
ks.test(m1, pnorm, mean = mean(m1), sd = sd (m1))
  # auto correlation function acf
acf(m1)

   ## all test showed not normally distributed, so we need to check the residual

 ## Check residual

means = param.est(m1.amos)$mean
m1.resid = m1-rep (means, seg.len(m1.amos))
   # Shapiro and Kolmogorov test
shapiro.test(m1.resid)
ks.test(m1.resid, pnorm, mean= mean(m1.resid), 
        sd = sd (m1.resid))
   # Q and Q plot
qqnorm(m1.resid)
qqline(m1.resid)


## Check assumptions for all simulated data using residual check

#1)  GP appointment
meansGP = param.est(GPvisitsWeekNov1718)$mean
GP.resid = GP.default-rep (means, seg.len(GP.scale))
 # Shapiro and Kolmogorov test
shapiro.test(GP.resid)
ks.test(m1.resid, pnorm, mean= mean(m1.resid), 
        sd = sd (m1.resid))
# Q and Q plot
qqnorm(m1.resid)
qqline(m1.resid)
 