#-----WE REFER TO ALEXIOS R PACKGAGES RUGARCH AND RMGARCH -----------------------------
#http://www.unstarched.net/2013/01/03/the-garch-dcc-model-and-2-stage-dccmvt-estimation/
#--------------------------------------------------------------------------------------
library(rmgarch)
library(rugarch)
library(xts)
library(parallel)
# Import data is already done 
# y is vector with 2 log returns SEEPEX and HUPX
# while X is the same vector, just index are dates
#--------------------------------------------------------------------------------------
TT = dim(y)[1]
y <- y[1:730,]
Dat = X
# define a DCCspec object: 2 stage estimation should usually always use

xspec = ugarchspec(mean.model = list(armaOrder = c(1, 1)), variance.model = list(garchOrder = c(1,1), model = 'eGARCH'), distribution.model = 'norm')
uspec = multispec(replicate(2, xspec))
spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
spec1a = dccspec(uspec = uspec, dccOrder = c(1, 1), model='aDCC', distribution = 'mvnorm')
spec2 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvlaplace')
spec2a = dccspec(uspec = uspec, dccOrder = c(1, 1), model='aDCC', distribution = 'mvlaplace')


cl = makePSOCKcluster(2)
multf = multifit(uspec, Dat, cluster = cl)

fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
fit1a = dccfit(spec1a, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
fit2 = dccfit(spec2, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
fit2a = dccfit(spec2a, data = Dat, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)


spec = gogarchspec(mean.model = list(armaOrder = c(0, 0), 
                                     include.mean =FALSE),
                   variance.model = list(model = "sGARCH", 
                                         garchOrder = c(1,1)) , 
                   distribution.model =  "mvnorm"
)
# First Estimate a QML first stage model (multf already estimated). Then
# estimate the second stage shape parameter.
spec3 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvt')
fit3 = dccfit(spec3, data = Dat, fit.control = list(eval.se = FALSE), fit = multf)
# obtain the multivariate shape parameter:
mvt.shape = rshape(fit3)
# Plug that into a fixed first stage model and iterate :
mvt.l = rep(0, 6)
mvt.s = rep(0, 6)
mvt.l[1] = likelihood(fit3)
mvt.s[1] = mvt.shape
for (i in 1:5) {
  xspec = ugarchspec(mean.model = list(armaOrder = c(1, 1)), variance.model = list(garchOrder = c(1,1), model = 'eGARCH'), distribution.model = 'std', fixed.pars = list(shape = mvt.shape))
  spec3 = dccspec(uspec = multispec(replicate(2, xspec)), dccOrder = c(1,1), distribution = 'mvt')
  fit3 = dccfit(spec3, data = Dat, solver = 'solnp', fit.control = list(eval.se = FALSE))
  mvt.shape = rshape(fit3)
  mvt.l[i + 1] = likelihood(fit3)
  mvt.s[i + 1] = mvt.shape
}
# Finally, once more, fixing the second stage shape parameter, and
# evaluating the standard errors
xspec = ugarchspec(mean.model = list(armaOrder = c(1, 1)), variance.model = list(garchOrder = c(1,1), model = 'eGARCH'), distribution.model = 'std', fixed.pars = list(shape = mvt.shape))
spec3 = dccspec(uspec = multispec(replicate(2, xspec)), dccOrder = c(1, 1), distribution = 'mvt', fixed.pars = list(shape = mvt.shape))
fit3 = dccfit(spec3, data = Dat, solver = 'solnp', fit.control = list(eval.se = TRUE), cluster = cl)


#----asymmetric DCC (MVT) model

xspec = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(garchOrder = c(1,1), model = "eGARCH"),  distribution.model = "norm")
spec3a  = dccspec(uspec = multispec( replicate(2, xspec) ), dccOrder = c(1,1), distribution = "mvt", model="aDCC")
fit3a = dccfit(spec3a, data = Dat, fit.control = list(eval.se=FALSE), fit = multf)
# obtain the multivariate shape parameter:
mvtx.shape = rshape(fit3a)
# Plug that into a fixed first stage model and iterate :
mvtx.l = rep(0, 6)
mvtx.s = rep(0, 6)
mvtx.l[1] = likelihood(fit3)
mvtx.s[1] = mvtx.shape
for(i in 1:5){
  xspec = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(garchOrder = c(1,1), model = "eGARCH"),  distribution.model = "std", fixed.pars = list(shape=mvtx.shape))
  spec3a = dccspec(uspec = multispec( replicate(2, xspec) ), dccOrder = c(1,1), model="aDCC", distribution = "mvt")
  fit3a = dccfit(spec3a, data = Dat, solver = "solnp", fit.control = list(eval.se=FALSE))
  mvtx.shape = rshape(fit3a)
  mvtx.l[i+1] = likelihood(fit3a)
  mvtx.s[i+1] = mvtx.shape
}
# Finally, once more, fixing the second stage shaoe parameter, and evaluating the standard errors
xspec = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(garchOrder = c(1,1), model = "eGARCH"),  distribution.model = "std", fixed.pars = list(shape=mvtx.shape))
spec3a = dccspec(uspec = multispec( replicate(2, xspec) ), dccOrder = c(1,1), model="aDCC", distribution = "mvt", fixed.pars=list(shape=mvtx.shape))
fit3a = dccfit(spec3a, data = Dat, solver = "solnp", fit.control = list(eval.se=TRUE), cluster = cl)



#STOP CLUSTER PARALELE
stopCluster(cl)

# FORECAST
fcast1 <- dccforecast(fit1, n.ahead=1)
fcast1.cor <- rcor(fcast1)

# PLOTING
library(timeSeries)
R1 = rcor(fit1)
R2 = rcor(fit2)
R3 = rcor(fit3)
D = rownames(X)
# D = y1$datum1[1:730]
colx = c(colors()[24], colors()[33], colors()[139])

par(mfrow = c(1,1))
RR = timeSeries(cbind(R1[2,1,],R2[2,1,],R3[2,1,]), D)
plot(RR[,1], ylab = "cor", col = colx[1], lty=1 ,lwd=1)
for(i in 2:3) lines(RR[,i], col = colx[i], lty = i, lwd=1+i/10)
title("DCC HUPX-SEEPEX Returns")
legend("bottomright", c("DCC(Normal)", "DCC(Laplace)", "DCC(Student T)"), col = colx, lty=1:3, bty="n")
# We further continue in DCC2 file

#--------------------------------------------------------------------------------------

#This is from Jon Danielsson book 

fit = gogarchfit(spec = spec, data = y)

show(fit) 
xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE))
uspec = multispec(replicate(2, xspec))
spec = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
res = dccfit(spec, data = y)
H=res@mfit$H
DCCrho=vector(length=dim(y)[1])
for(i in 1:dim(y)[1]){
  DCCrho[i] =  H[1,2,i]/sqrt(H[1,1,i]*H[2,2,i])
}
matplot(DCCrho,type='l',las=1,lty=1,col=2:3,ylab="")
mtext("Corr",side=2,line=0.5,at=1,las=1,cex=0.9)
legend("Log returns","DCC",lty=1,col=2:3,bty="n",cex=0.9)
acf(y[,1])
forecast = gogarchforecast(fit, n.ahead = 1)
#--------------------------------------------------------------------------------------