

# https://datasciencegenie.com/computing-the-portfolio-var-using-copulas/

# Serilerin grafiðini çizme I. yöntem
plot.ts(BIST100)
plot.ts(USD)
tm<-cbind(BIST100,USD)
plot.ts(tm, col="blue")


plot(BIST100,USD,pch='o')
abline(lm(USD~BIST100),col='red',lwd=1)
cor(BIST100,USD,method='spearman')
library(VineCopula)





# Serilerin grafiðini çizme II. yöntem
plot(BIST100)
lines(BIST100, type = "o", col = "blue")
plot(USD)
lines(USD, type = "o", col = "blue")

# Histogram grafiðini çiz
hist(BIST100, breaks = 20, col = "green", density = 20)
hist(USD, breaks = 20, col = "green", density = 20)



par(mfrow=c(1,2))
hist(BIST100,main="Daily Returns of BIST100 Histogram",xlab="BÝST100 Daily Returns",breaks = 60, col = "green", density = 10)
hist(USD,main="Daily Returns of USD Histogram",xlab="USD Daily Returns",breaks = 60, col = "green", density = 10)
par(mfrow=c(1,1))



#Correlation
library(psych)
df <- cbind(BIST100,USD)
cor(df,method='spearman')
pairs.panels(df)
u <- pnorm(df)
pairs.panels(u)
library(rgl)
plot3d(u[,1],u[,2],pch=20,col='blue')
plot3d(u[,1],u[,2],u[,3],pch=20,col='navyblue')

pairs.panels(df)
cor(df,meth='spearman')


shapiro.test(BIST100)
shapiro.test(USD)



# Descriptives of BIST100 Returns

library(e1071) #package for skewness and kurtosis functions

mean(BIST100)
sd(BIST100)
skewness(BIST100)
kurtosis(BIST100)

# Descriptives of USD Returns

mean(USD)
sd(USD)
skewness(USD)
kurtosis(USD)



library(fitdistrplus)


par(mfrow=c(1,2))
descdist(BIST100, discrete=FALSE, boot=100)
descdist(USD, discrete=FALSE, boot=100)
par(mfrow=c(1,1))




# Fitting the t distribution to BIST100 Returns

library(QRM)

tfit1 <- fit.st(BIST100)
#tfit1
tpars1<-tfit1[[2]] #vector of parameter estimates

nu1<-tpars1[1]
mu1<-tpars1[2]
sigma1<-tpars1[3]

par(mfrow=c(1,2))
z<-seq(-0.1,0.1,0.001)
hist(BIST100,nclass=20,probability=TRUE,main = "Histogram of BIST100 Returns",xlab="BIST100 Returns",breaks = 60, col = "green", density = 10) # Ben ekledim.
a<-Vectorize(function(z) dnorm(z, mean=mean(BIST100), sd=sd(BIST100)))
curve(a, col="red", lwd=3,  add=TRUE)
b<-Vectorize(function (z) dt((z-mu1)/sigma1,df=nu1)/sigma1)
curve(b, col="blue", lwd=3,  add=TRUE)


# Fitting the t distribution to USD Returns

tfit2 <- fit.st(USD)
#tfit2
tpars2 <- tfit2[[2]] #vector of parameter estimates

nu2<-tpars2[1]
mu2<-tpars2[2]
sigma2<-tpars2[3]

z<-seq(-0.1,0.1,0.001)
hist(USD,nclass=20,probability=TRUE,main = "Histogram of USD Returns",xlab="USD Returns",breaks = 60, col = "green", density = 10)
a<-Vectorize(function(z) dnorm(z, mean=mean(USD), sd=sd(USD)))
curve(a, col="red", lwd=3,  add=TRUE)
b<-Vectorize(function (z) dt((z-mu2)/sigma2,df=nu2)/sigma2)
curve(b, col="blue", lwd=3,  add=TRUE)
par(mfrow=c(1,1))





#Kolmogorov-Smirnov Test for BIST100

ks.test(BIST100, "pnorm",mean=mean(BIST100),sd=sd(BIST100))

shiftedt<-function(k,nu,mu,sigma){pt((k-mu)/sigma,df=nu)}
ks.test(BIST100, "shiftedt",nu=nu1,mu=mu1,sigma=sigma1)





#Kolmogorov-Smirnov Test for USD
ks.test(USD, "pnorm",mean=mean(BIST100),sd=sd(BIST100))

shiftedt<-function(k,nu,mu,sigma){pt((k-mu)/sigma,df=nu)}
ks.test(USD, "shiftedt",nu=nu2,mu=mu2,sigma=sigma2)




par(mfrow=c(1,2))
plot(BIST100,USD,xlab="BIST100 Daily Returns",ylab="USD Daily Returns",main="BIST100 and USD Daily Returns Scatterplot")

BIST100trans<-ecdf(BIST100)(BIST100)
USDtrans<-ecdf(USD)(USD)
plot(BIST100trans,USDtrans,xlab="CDF BIST100 Daily Returns",ylab="CDF USD Daily Returns",main="BIST100 and USD Returns CDFs")
par(mfrow=c(1,1))



library(VineCopula)# Kullanýlmak istenen paketin(VineCopula) yüklenmesi.
FittedCopula<-BiCopSelect(BIST100trans,USDtrans,familyset=NA)#Ýki veri serisi arasýndaki en uygun copula fonksiyonunun belirlenmesi.
summary(FittedCopula)#Copula fonksiyonunun özelliklerinin gösterilmesi.

# Gerçek günlük getirilerin CDF grafiðinin çizilmesi.
plot(BIST100trans,USDtrans,xlab="CDF BIST100 Daily Returns",ylab="CDF USD Daily Returns",main="Actual and Sampled Returns CDFs")

CopulaSim<-BiCopSim(length(BIST100), FittedCopula)#Copula fonksiyonunun özelliklerine göre yeni bir veri üretilmesi.
points(CopulaSim, col = rgb(0,0,1,alpha=0.3), pch=16) #Yeni üretilen verinin CDF grafiðinin çizilmesi.






#Portfolio Weights
w1<-0.5 
w2<-0.5

library(copula)
library(VC2copula)

n<-1000000
Percentiles<-c()
for(i in 1:20)
  sim<-rCopula(n,copulaFromFamilyIndex(FittedCopula,FittedCopula))
sim2<-cbind(qt(sim[,1],df=nu1)*sigma1+mu1,qt(sim[,2],df=nu2)*sigma2+mu2)
Percentiles<-c(Percentiles,sort(sim2 %*% c(w1,w2))[0.01*n])

VaR<-mean(Percentiles)
VaR


#Computing VaR with the Multivariate Normal Distribution
n<-1000000
Percentiles<-c() 
for(i in 1:20)
sim<-rCopula(n,normalCopula(param = cor(BIST100,USD), dim = 2))
sim2<-cbind(qnorm(sim[,1],mean=mean(BIST100),sd=sd(BIST100)),qnorm(sim[,2],mean=mean(USD),sd=sd(USD)))
Percentiles<-c(Percentiles,sort(sim2 %*% c(w1,w2))[0.01*n])

VaR<-mean(Percentiles)
VaR







#t-Copula

rho<-0.7
df<-2

t_Cop<- mvdc(tCopula(rho, dim = 2,df=df), margins=c("norm","norm"),paramMargins=list(list(mean=0, sd=1),list(mean=0, sd=1)))

x1 <- x2 <- seq(-3, 3, length= 200)

v<-c()
for (i in x1)
{for (j in x2){v<-c(v,dMvdc(c(i,j), t_Cop))}}
f<-t(matrix(v,nrow=200,byrow=TRUE))

plot_ly(x=x1,y=x2,z=f,type = "contour")%>%layout(xaxis=list(title="x1"),yaxis=list(title="x2"),title=paste("t-Copula with ??=",as.character(rho)," and ",as.character(df),"df"))

