###################################################################################
# run this code only once 
library(MASS)   
###################################################################################
#
###################################################################################
# Problem 1: Testing the Multinomial Model with Equal Probabilities
id<-20456458
set.seed(id)
k<-sample(5:9,1)   # randomly choose number of categories for Multinomial data
p<-sample(1:9,k,replace=TRUE)
p<-p/sum(p)       # choose random probabilities which must sum to one
y<-rmultinom(1,150,p)    # generate random data
e<- rep(150/k, k)  # calculate expected frequencies assuming equal probabilities for each category
# print table of observed and expected frequencies
cat("Table of Observed and Expected Frequencies ")
print(data.frame("Category" = rbind(y[,1],e), row.names = c("Observed", "Expected")),digits=4) 
# observed values of likelihood ratio test statistic and Goodness of Fit test statistic
# and corresponding p-values
df<-k-1      # degrees of freedom for the Chi-squared distribution
lambda<-2*sum(y*log(y/e))
pvalue<-1-pchisq(lambda,df)
cat("Observed value of likelihood ratio statistic = ", lambda)
cat("with p-value = ",pvalue, "and degrees of freedom = ",df)
pearson<-sum(((y-e)^2)/e) 
pvalue<-1-pchisq(pearson,df)  
cat("Observed value of Goodness of Fit statistic = ", pearson)
cat("with p-value = ", pvalue, "and degrees of freedom = ",df)
###################################################################################
#
###################################################################################
# Problem 2: Testing the Goodness of Fit of a Poisson Model
set.seed(id)
model<-sample(c(1:4),1)
cat("Model = ", model)
# Data are randomly generated from one of four different models all with mean 4
# Model=1: Poisson(4) distribution 
# Model=2:  Negative Binomial(3,3/7)
# Model=3:  G(4,1) distribution and discretized 
# Model=4:  Gamma(3,4/3) distribution and discretized 
if (model==1) {
  y<-rpois(150,4)    # 150 observations from Poisson(4)
} else if (model==2) {
  y<-rnbinom(150,3,3/7)   # 150 observations from NB(3,3/7)
} else if (model==3) {
  y<-round( rnorm(150,4,1))   # 150 observations from G(4,1) rounded
  y[y<0]<-0   # convert any negative observations to 0
} else if (model==4) {
  y<-round(rgamma(150,3,3/4))    # 150 observations from Gamma(3,4/3) rounded
} 
ymin<-min(y)
ymax<-max(y)
# determine categories and frequencies for the data
data<-table(c(y, ymin:ymax))-1      # Done to ensure all categories are accounted for
f<-as.numeric(data)     # frequencies
cat<-as.numeric(names(data))
# determine the maximum likelihood estimate of theta which is the sample mean calculated 
# from the frequency table
thetahat<-sum(cat*f)/150 
# determine the expected frequencies 
e<-dpois(cat,thetahat)*150   #expected frequencies for Poisson data
#frequency for ymin must be sum for y<=ymin
e[1]<-ppois(ymin,thetahat)*150
ncat<-length(e)
# frequency for ymax must be sum of frequencies for y>=ymax
e[ncat]<- ppois(ymax- 1,thetahat, lower = F)*150     
# Table of Observed and expected frequencies
data<-rbind("y" = ymin:ymax, "observed" = f, "expected" = e)    
# print table of observed and expected frequencies
cat("Table of Observed and Expected Frequencies ")
print(data,digits=4) 
# Expected frequencies must all be at least 5 to apply tests. Collapse categories if necessary.
nbins<-ncol(data)
while(data[3, nbins] < 5){
  data[2:3, nbins - 1]<-data[2:3, nbins - 1] + data[2:3, nbins]
  data<-data[, -nbins]
  nbins<-nbins - 1
}
nbins<-1
while(data[3, nbins] < 5){
  data[2:3, nbins + 1]<-data[2:3, nbins + 1] + data[2:3, nbins]
  data<-data[, -nbins]
}
cat("Table of Observed and Expected Frequencies ")
# print table of observed and expected frequencies
print(data,digits=4) 
# observed values of likelihood ratio test statistic and Goodness of Fit test statistic
# and corresponding p-values
df = ncol(data)-2      # degress of freedom for the Chi-squared distribution
f<-data[2,]
e<-data[3,]
lambda<-2*sum(f*log(f/e))
pvalue<-1-pchisq(lambda,df)
cat("Observed value of likelihood ratio statistic = ", lambda)
cat("with p-value = ",pvalue, "and degrees of freedom = ",df)
pearson<-sum(((f-e)^2)/e) 
pvalue<-1-pchisq(pearson,df)  
cat("Observed value of Goodness of Fit statistic = ", pearson)
cat("with p-value = ", pvalue, "and degrees of freedom = ",df)
###################################################################################
#
###################################################################################
# Problem 3: Testing for Independence in Two Way Tables
set.seed(id)
# generate data for a two way table by first simulating bivariate data  
# from the Bivariate Normal distribution and then discretize the data
# Random uniform between -0.75 and 0.75
corrCoef<-runif(1, -0.75, 0.75)
sigma<-max(id %% 10, 1)
# Last digit of UWID using modulo, minimum value of 1.
mu1<-max(id %% 100 - id %% 10, 20)
# (Second last digit*10) is extracted here, minimum value of 20
mu2<-max(id %% 1000 - id %% 100, 30)
# (Third last digit*100) is extracted here, minimum value of 30
VarCovar<-cbind(c(sigma^2, corrCoef*sigma^2), c(corrCoef*sigma^2, sigma^2))
# Simulate data from a bivariate Normal
n<-sample(c(100:200),1)    # n =  sample size
cat("Number of observations = ",n)
data2<-mvrnorm(n, mu = c(mu1, mu2), Sigma = VarCovar)
# Create smoker/non-smoker variable by mapping 1 to smoker and 2 to non-smoker
data3 = as.data.frame(data2)
data3[, 1]<-ifelse(data3[,1] < median(data3[,1]), 1, 2)
data3[, 1]<-c("Smoker", "Non-smoker") [data3[,1]]
# Create tall/avg/short height variable by mapping 1 to tall, 2 to average and 3 to short
data3[, 2]<-floor((rank(data3[, 2])-0.1)/nrow(data3)*3) + 1
data3[, 2]<-c("Tall", "Average", "Short")[data3[, 2]]
data3[, 1]<-factor(data3[, 1])
data3[, 2]<-factor(data3[, 2])
colnames(data3)<-c("Smoker Indicator", "Height Indicator")
f<-table(data3)
cat("Table of Observed Frequencies:")
f
r<-margin.table(f,1)     # row totals
c<-margin.table(f,2)     # column totals
e<-outer(r,c)/sum(f)   # matrix of expected frequencies
cat("Table of Expected Frequencies:")
print(e,digits=4)
lambda<-2*sum(f*log(f/e))  # observed value of likelihood ratio statistic
df<-(length(r)-1)*(length(c)-1)  # degrees of freedom
pvalue<-1-pchisq(lambda,df)
cat("Observed value of likelihood ratio statistic = ", lambda)
cat("with p-value = ",pvalue, "and degrees of freedom = ",df)
pearson<-sum(((f-e)^2)/e) 
pvalue<-1-pchisq(pearson,df)  
cat("Observed value of Goodness of Fit statistic = ", pearson)
cat("with p-value = ", pvalue, "and degrees of freedom = ",df)
###################################################################################


