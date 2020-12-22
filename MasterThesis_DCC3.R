#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##
#################################################################################

#################################################################################
# DCC model
#################################################################################

# We have already defined X in IMPORT DATA

	Dat = X

	
	cnames = colnames(Dat)
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), variance.model = list(garchOrder = c(1,1), model = "eGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(2, uspec) ), dccOrder = c(1,1), 
	                distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se=FALSE))
	
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), 
	                   variance.model = list(garchOrder = c(1,1), model = "eGARCH"), 
	                   distribution.model = "norm")
	vspec = vector(mode = "list", length = 2)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:2){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), 
	                dccOrder = c(1,1),  distribution = "mvnorm",
	                fixed.pars = dccfix)
	
	presigma = tail( sigma(fit1 ), 2 )
	preresiduals = tail( residuals(fit1), 2 )
	prereturns = tail( as.matrix(Dat), 2 )
	sim1 = dccsim(fitORspec = fit1, n.sim = 1000, n.start = 100, m.sim = 2, startMethod = "unconditional", 
	              presigma = presigma, preresiduals = preresiduals, prereturns = prereturns, 
	              preQ = last(rcor(fit1, type = "Q"))[,,1], Qbar = fit1@mfit$Qbar, 
	              preZ = tail(fit1@mfit$stdresid, 1),
	              rseed = c(100, 200), mexsimdata = NULL, vexsimdata = NULL)
	
	
	sim2 = dccsim(fitORspec = spec2, n.sim = 1000, n.start = 100, m.sim = 2,
	              presigma = presigma, preresiduals = preresiduals, prereturns = prereturns, 
	              preQ = last(rcor(fit1, type = "Q"))[,,1], Qbar = fit1@mfit$Qbar, 
	              preZ = tail(fit1@mfit$stdresid, 1),
	              rseed = c(100, 200), mexsimdata = NULL, vexsimdata = NULL)
	
	sim3 = dccsim(fitORspec = fit1, n.sim = 1000, n.start = 100, m.sim = 2, startMethod = "sample", 
	              rseed = c(100, 200), mexsimdata = NULL, vexsimdata = NULL)
	
	
	rc1 = rcor(sim1, sim = 1)
	rc2 = rcor(sim2, sim = 1)
	rc3 = rcor(sim3, sim = 1)
	rh1 = rcov(sim1, sim = 2)
	rh2 = rcov(sim2, sim = 2)
	rh3 = rcov(sim2, sim = 2)
	
	options(width = 120)
	zz <- file("Output/test2e1.txt", open="wt")
	sink(zz)
	print( all.equal(first(rc1)[,,1], first(rc2)[,,1], first(rc3)[,,1], check.attributes = FALSE))
	print( all.equal(last(rc1)[,,1], last(rc2)[,,1], last(rc3)[,,1], check.attributes = FALSE))
	print( all.equal(first(rh1)[,,1], first(rh2)[,,1], first(rh3)[,,1], check.attributes = FALSE))
	print( all.equal(last(rh1)[,,1], last(rh2)[,,1], last(rh3)[,,1], check.attributes = FALSE))
	print( all.equal(head(fitted(sim1)), head(fitted(sim2)), head(fitted(sim3)), check.attributes = FALSE))
	print( all.equal(head(fitted(sim1, sim=2)), head(fitted(sim2, sim=2)), head(fitted(sim3, sim=2)), check.attributes = FALSE))
	sink(type="message")
	sink()
	close(zz)
	
	# Now some 1-ahead rolling tests
	
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), 
	                   variance.model = list(garchOrder = c(1,1), model = "apARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(3, uspec) ), dccOrder = c(1,1), 
	                distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, out.sample = 100, fit.control = list(eval.se=FALSE))
	
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), 
	                   variance.model = list(garchOrder = c(1,1), model = "apARCH"), 
	                   distribution.model = "norm")
	vspec = vector(mode = "list", length = 3)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:3){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), 
	                dccOrder = c(1,1),  distribution = "mvnorm",
	                fixed.pars = dccfix)
	
	rcovfilt = rcorfilt = rcovsim = rcorsim = array(NA, dim = c(3,3,100))
	T = dim(Dat)[1] - 100
	p = 2
	preQ = last(rcor(fit1, type = "Q"))[,,1]
	presigma = tail( sigma(fit1 ), 2 )
	preresiduals = tail( residuals(fit1), 2 )
	prereturns = tail( as.matrix(Dat[1:T,]), 2 )
	
	filt = dccfilter(spec2, data = Dat[1:T,], filter.control = list(n.old = T))
	for(i in 1:100){
	  preQ = last(rcor(filt, type = "Q"))[,,1]
	  presigma = tail( sigma(filt), 2 )
	  preresiduals = tail( residuals(filt), 2 )
	  prereturns = tail( as.matrix(Dat[1:(T+i-1),]), 2 )
	  preZ = tail(filt@mfilter$stdresid, 1)
	  # when (i == 1) equivalent to fit1
	  if(i == 1){
	    print( all.equal(last(rcor(filt))[,,1], last(rcor(fit1))[,,1]) )
	    print( all.equal(first(rcov(filt))[,,1], first(rcov(fit1))[,,1]) )
	    print( all.equal(tail(residuals(filt)), tail(residuals(fit1))) )
	  }
	  sim = dccsim(fitORspec = spec2, n.sim = 1, n.start = 0, m.sim = 100, 
	               presigma = presigma, preresiduals = preresiduals, 
	               prereturns = prereturns, preQ = preQ, preZ = preZ, Qbar = fit1@mfit$Qbar)
	  filt = dccfilter(spec2, data = Dat[1:(T+i),], filter.control = list(n.old = T))
	  tmp = matrix(0, 3, 3)
	  for(j in 1:100) tmp = tmp + sim@msim$simH[[j]][,,1]
	  tmp = tmp/100
	  rcovsim[,,i] = tmp
	  rcorsim[,,i] = cov2cor(tmp)
	  rcovfilt[,,i] = last(rcov(filt), 1)[,,1]
	  rcorfilt[,,i] = last(rcor(filt), 1)[,,1]
	  # check upto point T that it agress with fit:
	  #print(all.equal(last(rcov(filt), 2)[,,1], last(rcov(fit1), 1)[,,1]))
	  #print(all.equal(first(rcov(filt), 1)[,,1],first(rcov(fit1), 1)[,,1]))
	  # check upto point T that it agress with fit:
	  print(i)
	}
	
	postscript("Output/test2e1.eps", width = 16, height = 14)
	par(mfrow = c(2,2))
	plot(rcovfilt[1,2,], type = "l", main = paste(cnames[1], "-", cnames[2], sep = ""),
	     ylab = "covariance", xlab = "periods [out of sample]")
	lines(rcovsim[1,2,], col = 2, lty = 2)
	legend("topleft", legend = c("Filtered Covariance", "Simulated Covariance"), col = 1:2, lty = 1:2,
	       bty = "n")
	
	plot(rcovfilt[2,3,], type = "l", main = paste(cnames[2], "-", cnames[3], sep = ""),
	     ylab = "covariance", xlab = "periods [out of sample]")
	lines(rcovsim[2,3,], col = 2, lty = 2)
	legend("topleft", legend = c("Filtered Covariance", "Simulated Covariance"), col = 1:2, lty = 1:2,
	       bty = "n")
	
	plot(rcorfilt[1,2,], type = "l", main = paste(cnames[1], "-", cnames[2], sep = ""),
	     ylab = "correlation", xlab = "periods [out of sample]")
	lines(rcorsim[1,2,], col = 2, lty = 2)
	legend("topleft", legend = c("Filtered Correlation", "Simulated Correlation"), col = 1:2, lty = 1:2,
	       bty = "n")
	plot(rcorfilt[2,3,], type = "l", main = paste(cnames[2], "-", cnames[3], sep = ""),
	     ylab = "correlation", xlab = "periods [out of sample]")
	lines(rcorsim[2,3,], col = 2, lty = 2)
	legend("topleft", legend = c("Filtered Correlation", "Simulated Correlation"), col = 1:2, lty = 1:2,
	       bty = "n")
	dev.off()
	
	
	
	# Check Sim by fit and Sim by Spec using VAR:
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(3, uspec) ), dccOrder = c(1,1), 
	                distribution = "mvnorm", VAR = TRUE, lag = 2)
	fit1 = dccfit(spec1, data = Dat, out.sample = 100, fit.control = list(eval.se=FALSE))
	
	T = dim(Dat)[1] - 100
	VAR.fit = varxfit(Dat[1:T,], p = 2, postpad = "constant")
	vspec = vector(mode = "list", length = 3)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:3){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), VAR = TRUE, lag = 2, 
	                dccOrder = c(1,1),  distribution = "mvnorm",
	                fixed.pars = dccfix)
	
	presigma = tail( sigma(fit1 ), 2 )
	prereturns = as.matrix(tail(Dat[1:T,], 2))
	preresiduals = as.matrix(tail(residuals(fit1), 2))
	sim1 = dccsim(fitORspec = fit1, n.sim = 1, n.start = 0, m.sim = 100, startMethod = "sample", 
	              presigma = presigma, prereturns = prereturns, rseed = 1:100)
	
	preQ = last(rcor(fit1, type = "Q"))[,,1]
	preZ = tail(fit1@mfit$stdresid, 1)
	sim2 = dccsim(fitORspec = spec2, n.sim = 1, n.start = 0, m.sim = 100, startMethod = "sample", 
	              presigma = presigma, Qbar = fit1@mfit$Qbar, preQ = preQ, preZ = preZ, 
	              prereturns = prereturns, preresiduals = preresiduals, rseed = 1:100, VAR.fit = VAR.fit)
	
	
	rc1 = rcor(sim1, sim = 1)
	rc2 = rcor(sim2, sim = 1)
	rh1 = rcov(sim1, sim = 2)
	rh2 = rcov(sim2, sim = 2)
	
	options(width = 120)
	zz <- file("Output/test2e2.txt", open="wt")
	sink(zz)
	print( all.equal(first(rc1)[,,1], first(rc2)[,,1], check.attributes=FALSE) )
	print( all.equal(last(rc1)[,,1], last(rc2)[,,1], check.attributes=FALSE) )
	print( all.equal(first(rh1)[,,1], first(rh2)[,,1], check.attributes=FALSE) )
	print( all.equal(last(rh1)[,,1], last(rh2)[,,1], check.attributes=FALSE) )
	print( all.equal(head(fitted(sim1)), head(fitted(sim2)), check.attributes=FALSE) )
	print( all.equal(head(fitted(sim1, sim=2)), head(fitted(sim2, sim=2)), check.attributes=FALSE) )
	sink(type="message")
	sink()
	close(zz)
	
	
	
	# Forecast Tests 
	
	
	#---------------------------------------------------------------------------------------------------------------
	# ARMA-GARCH-DCC
	cnames = colnames(Dat)
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), variance.model = list(garchOrder = c(1,1), model = "gjrGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(3, uspec) ), dccOrder = c(1,1), 
	                distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, out.sample = 100, fit.control = list(eval.se=FALSE))
	T = dim(Dat)[1]-100
	
	forc = dccforecast(fit1, n.ahead = 1, n.roll = 0)
	forc2 = dccforecast(fit1, n.ahead = 1, n.roll = 10)
	
	# forcx = dccforecast(fit1, n.ahead = 1, n.roll = 99)
	# REM: n.roll = 10 == 11 roll (1 roll is the standard + 10 beyond)
	# Now filter
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), 
	                   variance.model = list(garchOrder = c(1,1), model = "gjrGARCH"), 
	                   distribution.model = "norm")
	vspec = vector(mode = "list", length = 2)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:2){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), 
	                dccOrder = c(1,1),  distribution = "mvnorm",
	                fixed.pars = dccfix)
	
	# 1-step ahead forecast on out-sample data is equivalent to filtering
	# Note the use of n.old so that filtering assumptions are observed based on
	# the size of the original dataset
	filt1 = dccfilter(spec2, data = Dat[1:(T+100),], filter.control = list(n.old = T))
	
	# check mean and sigma
	options(width = 120)
	zz <- file("Output/test2f1-2.txt", open="wt")
	sink(zz)
	print(all.equal(fitted(forc)[1,,1], filt1@model$mu[T+1,], fitted(filt1)[T+1,], check.attributes = FALSE))	
	print(all.equal(sigma(forc)[1,,1], filt1@model$sigma[T+1,], sigma(filt1)[T+1,], check.attributes = FALSE))	
	print(all.equal(fitted(forc2)[1,,11], filt1@model$mu[T+11,], fitted(filt1)[T+11,], check.attributes = FALSE))	
	print(all.equal(rcor(forc)[[1]][,,1], rcor(filt1)[,,T+1], check.attributes = FALSE))	
	print(all.equal(rcor(forc2)[[11]][,,1], rcor(filt1)[,,T+11], check.attributes = FALSE))	
	print(all.equal(sigma(forc2)[,,11], filt1@model$sigma[T+11,], sigma(filt1)[T+11,], check.attributes = FALSE))	
	sink(type="message")
	sink()
	close(zz)
	
	
	#---------------------------------------------------------------------------------------------------------------
	# ARMA-GARCH-DCC-2
	cnames = colnames(Dat)
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(3, uspec) ), dccOrder = c(1,1), 
	                distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, out.sample = 500, fit.control = list(eval.se=FALSE))
	T = dim(Dat)[1]-500
	
	forc = dccforecast(fit1, n.ahead = 1, n.roll = 499)
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	vspec = vector(mode = "list", length = 3)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:3){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), 
	                dccOrder = c(1,1),  distribution = "mvnorm",
	                fixed.pars = dccfix)
	filt1 = dccfilter(spec2, data = Dat[1:(T+500),], filter.control = list(n.old = T))
	
	###########################################################################
	# DCC Filter vs rolling Forecast
	# The reason for the build up of the difference between the filter and rolling
	# forecast methods is that the filter uses a fixed estimate of the covariance
	# matrix initialization at time T (e.g. n.old), whilst the rolling forecast method
	# updates the initialization value at each iteration (T+1) to include the data upto 
	# point T.
	# plot(rcor(filt1)[1,2,(T+1):(T+500)], type = "l")
	# lines(sapply(rcor(forc), FUN = function(x) x[1,2,1]), col = 2)
	
	
	#---------------------------------------------------------------------------------------------------------------
	# VAR-GARCH-DCC
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(2, uspec) ), VAR = TRUE, lag = 2, dccOrder = c(1,1), 
	                distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, out.sample = 100, fit.control = list(eval.se=FALSE))
	T = dim(Dat)[1]-100
	
	forc = dccforecast(fit1, n.ahead = 1, n.roll = 10)
	# REM: n.roll = 10 == 11 roll (1 roll is the standard + 10 beyond)
	# Now filter (first the VAR so it is the same as in fit)
	VAR.filt = varxfilter(Dat[1:(T+100),], p = 2, Bcoef = fit1@model$varcoef, postpad = "constant")
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	vspec = vector(mode = "list", length = 2)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:2){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), VAR = TRUE, lag = 2, 
	                dccOrder = c(1,1),  distribution = "mvnorm",
	                fixed.pars = dccfix)
	
	
	filt1 = dccfilter(spec2, data = Dat[1:(T+100),], filter.control = list(n.old = T), varcoef = VAR.filt$Bcoef)
	
	# check mean and sigma
	options(width = 120)
	zz <- file("Output/test2f2.txt", open="wt")
	sink(zz)
	print(all.equal(as.numeric(fitted(forc)[,,1]), as.numeric(fitted(filt1)[T+1,])))
	print(all.equal(as.numeric(sigma(forc)[,,1]), as.numeric(sigma(filt1)[T+1,])))
	print(all.equal(as.numeric(fitted(forc)[,,11]), as.numeric(fitted(filt1)[T+11,])))
	print(all.equal(as.numeric(sigma(forc)[,,11]), as.numeric(sigma(filt1)[T+11,])))
	print(all.equal(rcor(forc)[[1]][,,1], rcor(filt1)[,,T+1]))
	# the difference is explained in the note above
	print(all.equal(rcor(forc)[[11]][,,1], rcor(filt1)[,,T+11]))
	
	sink(type="message")
	sink()
	close(zz)
	
	#---------------------------------------------------------------------------------------------------------------
	# VAR(1)-GARCH(1,1)-aDCC(1,2) n.ahead = 1 n.roll = 99
	# Dat = dji30retw[,1:10]
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(2, uspec) ), VAR = TRUE, lag = 1, 
	                dccOrder = c(1,2), model="aDCC", distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, out.sample = 100, solver = "solnp", fit.control = list(eval.se=FALSE))
	
	
	T = dim(Dat)[1]-100
	
	forc = dccforecast(fit1, n.ahead = 1, n.roll = 99)
	# REM: n.roll = 99 == 100 rolls (1 roll is the standard + 99 beyond)
	# Now filter (first the VAR so it is the same as in fit)
	VAR.filt = varxfilter(Dat[1:(T+100),], p = 1, Bcoef = fit1@model$varcoef, postpad = "constant")
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	vspec = vector(mode = "list", length = 2)
	midx = fit1@model$midx
	mpars = fit1@model$mpars
	for(i in 1:2){
	  vspec[[i]] = uspec
	  setfixed(vspec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit1, "dcc"))
	spec2 = dccspec(uspec = multispec( vspec ), VAR = TRUE, lag = 1, 
	                dccOrder = c(1,2), model="aDCC", distribution = "mvnorm",
	                fixed.pars = dccfix)
	
	filt1 = dccfilter(spec2, data = Dat[1:(T+100), ], filter.control = list(n.old = T), varcoef = VAR.filt$Bcoef)
	
	# check mean and sigma
	options(width = 120)
	zz <- file("Output/test2f3.txt", open="wt")
	sink(zz)
	print(all.equal(fitted(forc)[,,1], as.numeric(filt1@model$mu[T+1,]), check.attributes = FALSE))	
	print(all.equal(fitted(forc)[,,1], as.numeric(fitted(filt1)[T+1,]), check.attributes = FALSE))	
	print(all.equal(sigma(forc)[,,1], filt1@model$sigma[T+1,], check.attributes = FALSE))	
	print(all.equal(sigma(forc)[,,1], as.numeric(sigma(filt1)[T+1,]), check.attributes = FALSE))	
	print(all.equal(fitted(forc)[,,11], as.numeric(filt1@model$mu[T+11,]), check.attributes = FALSE))	
	print(all.equal(fitted(forc)[,,11], as.numeric(fitted(filt1)[T+11,]), check.attributes = FALSE))	
	print(all.equal(sigma(forc)[,,11], filt1@model$sigma[T+11,], check.attributes = FALSE))	
	print(all.equal(sigma(forc)[,,11], as.numeric(sigma(filt1)[T+11,]), check.attributes = FALSE))	
	print(all.equal(rcor(forc)[[1]][,,1], rcor(filt1)[,,T+1], check.attributes = FALSE))	
	print(all.equal(rcor(forc)[[11]][,,1], rcor(filt1)[,,T+11], check.attributes = FALSE))	
	
	sink(type="message")
	sink()
	close(zz)
	
	#################################
	# data(dji30ret)
	# Dat = dji30ret[,1:10]
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
	                   variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
	                   distribution.model = "norm")
	spec1 = dccspec(uspec = multispec( replicate(2, uspec) ), VAR = TRUE, lag = 1, 
	                dccOrder = c(1,2),  distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, out.sample = 100, solver = "solnp", fit.control = list(eval.se=FALSE))	
	
	forc1 = dccforecast(fit1, n.ahead = 1500, n.roll = 0)
	
	forc2 = dccforecast(fit1, n.ahead = 1, n.roll = 99)
	
	
	plot(forc1, which = 1)
	#plot(forc2, which = 1, series = 1:3)
	#plot(forc2, which = 4, series = 1:3)
	#plot(forc1, which = 4, series = 2:3)
	#plot(forc2, which = 5)
	#plot(forc1, which = 5)
	
	options(width = 120)
	zz <- file("Output/test2f4.txt", open="wt")
	sink(zz)
	# forecast R --> unconditional R (Rbar) as N --> large
	UQ = fit1@mfit$Qbar*(1-sum(coef(fit1, "dcc")))
	Rbar = UQ/(sqrt(diag(UQ)) %*% t(sqrt(diag(UQ))))
	print(all.equal(Rbar, forc1@mforecast$R[[1]][,,1500]))
	#############################################
	sink(type="message")
	sink()
	close(zz)
	
	
	
	# Roliing Tests 
	
	# data(dji30ret)	
	# Dat = dji30ret[, 1:5, drop = FALSE]
	Dat=X3
	uspec1 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "apARCH"), 
	                    distribution.model = "norm")
	uspec2 = ugarchspec(mean.model = list(armaOrder = c(2,0)), variance.model = list(model = "gjrGARCH"), 
	                    distribution.model = "norm")
	uspec3 = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec4 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(2, 1)), 
	                    distribution.model = "norm")
	uspec5 = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(model = "eGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "norm")
	
	uspec = c( uspec1, uspec5 )
	spec = dccspec(uspec = multispec( uspec ), dccOrder = c(1,1), distribution = "mvnorm")
	
	roll = dccroll(spec, data = Dat, n.ahead = 1, forecast.length = 30, refit.every = 2, 
	               refit.window = "moving", solver = "solnp", 
	               fit.control = list(eval.se = FALSE), solver.control = list(), 
	               save.fit = FALSE, save.wdir = NULL)
	
	options(width = 120)
	zz <- file("Output/test2g3-2.txt", open="wt")
	sink(zz)
	show(roll)
	cat("\nAverage Log-Likelihood Across Rolls:\n")
	print( round( likelihood(roll)/sapply(roll@model$rollind, FUN = function(x) length(x)), 2) )
	sink(type="message")
	sink()
	close(zz)
	
	postscript("Output/test2g3-2.eps", width = 12, height = 12)
	plot(roll, which = 4)
	dev.off()
	
	postscript("Output/test2g3-2.eps", width = 12, height = 12)
	plot(roll, which = 5)
	dev.off()
	
	postscript("Output/test2g1-3.eps", width = 12, height = 12)
	plot(roll, which = 3)
	dev.off()
	
	# simulation tests
	
	# Normal
	# Dat = dji30ret[, 1:5, drop = FALSE]
	uspec1 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec2 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec3 = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec4 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "norm")
	uspec5 = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "norm")
	
	uspec = c(  uspec2, uspec3 )
	spec1 = dccspec(uspec = multispec( uspec ), dccOrder = c(1,1), distribution = "mvnorm")
	
	fit1 = dccfit(spec1, data = Dat, solver = "solnp", fit.control = list(eval.se=FALSE), cluster = NULL)
	
	forc = dccforecast(fit1, n.ahead = 1)
	
	sim1 = dccsim(fit1, n.sim = 1, m.sim = 1000, startMethod = "sample",  rseed = 1:1000)
	
	uspec1 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "ged", fixed.pars=list(shape=1))
	uspec2 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "ged", fixed.pars=list(shape=1))
	uspec3 = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "ged", fixed.pars=list(shape=1))
	uspec4 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "ged", fixed.pars=list(shape=1))
	uspec5 = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "ged", fixed.pars=list(shape=1))
	spec2 = dccspec(uspec = multispec( uspec ), dccOrder = c(1,1), distribution = "mvlaplace")
	
	fit2 = dccfit(spec2, data = Dat, solver = "solnp", fit.control = list(eval.se=FALSE))
	
	sim2 = dccsim(fit2, n.sim = 1, m.sim = 1000, startMethod = "sample",  rseed = 1:1000)
	
	# gauge the uncertainty around 1-ahead
	forcM = fitted(forc)[,,1]
	
	sim1M = t(sapply(sim1@msim$simX, FUN = function(x) x))
	# equivalent: t(sapply(1:1000, FUN = function(i) fitted(sim1, sim=i)))
	sim2M = t(sapply(sim2@msim$simX, FUN = function(x) x))
	
	cnames = fit1@model$modeldata$asset.names
	mvnH = apply(sim1M, 2, FUN = function(x) hist(x, breaks="fd", plot=FALSE))
	mvlH = apply(sim2M, 2, FUN = function(x) hist(x, breaks="fd", plot=FALSE))
	minX = maxX = rep(0, 2)
	maxY = rep(0, 2)
	for(i in 1:2){
	  maxX[i] = max(c(mvnH[[i]]$breaks, mvlH[[i]]$breaks))
	  minX[i] = min(c(mvnH[[i]]$breaks, mvlH[[i]]$breaks))
	  maxY[i] = max(c(mvnH[[i]]$counts, mvlH[[i]]$counts))
	}
	
	postscript("Output/test2h1-1.eps", width = 16, height = 14)
	par(mfrow = c(2,3))
	for(i in 1:2){
	  plot(mvnH[[i]], xlim = c(minX[i], maxX[i]), ylim = c(0, maxY[i]), border = "steelblue", xlab = "", main = cnames[i])
	  plot(mvlH[[i]], border = "tomato", add = TRUE)
	  abline(v = forcM[i], col = "black", lwd=1)
	  legend("topleft", c("MVN", "MVL"), col = c("steelblue", "tomato"), bty="n", lty=c(1,1))
	  legend("topright", c("Forecast"), col = c("black"), bty="n", lty=1)
	}
	dev.off()
	
	# one more with Student (QML first stage)
	uspec1 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec2 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec3 = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "norm")
	uspec4 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "norm")
	uspec5 = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "norm")
	spec3 = dccspec(uspec = multispec( uspec ), dccOrder = c(1,1), distribution = "mvt")
	fit3 = dccfit(spec3, data = Dat, solver = "solnp", fit.control = list(eval.se=FALSE))
	
	# now go back and fix the shape
	shp = unname(coef(fit3, "dcc")["mshape"])
	uspec1 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "std", fixed.pars = list(shape=shp))
	uspec2 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "std", fixed.pars = list(shape=shp))
	uspec3 = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(model = "sGARCH"), 
	                    distribution.model = "std", fixed.pars = list(shape=shp))
	uspec4 = ugarchspec(mean.model = list(armaOrder = c(1,0)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "std", fixed.pars = list(shape=shp))
	uspec5 = ugarchspec(mean.model = list(armaOrder = c(1,1)), variance.model = list(model = "sGARCH",
	                                                                                 garchOrder = c(1, 1)), 
	                    distribution.model = "std", fixed.pars = list(shape=shp))
	spec3 = dccspec(uspec = multispec( uspec ), dccOrder = c(1,1), distribution = "mvt", fixed.pars=list(mshape=shp))
	fit3 = dccfit(spec3, data = Dat, solver = "solnp", fit.control = list(eval.se=FALSE))
	
	sim3 = dccsim(fit3, n.sim = 1, m.sim = 1000, startMethod = "sample",  rseed = 1:1000)
	
	sim3M = t(sapply(sim3@msim$simX, FUN = function(x) x))
	
	mvtH = apply(sim3M, 2, FUN = function(x) hist(x, breaks="fd", plot=FALSE))
	minX = maxX = rep(0, 2)
	maxY = rep(0, 2)
	for(i in 1:2){
	  maxX[i] = max(c(mvnH[[i]]$breaks, mvtH[[i]]$breaks))
	  minX[i] = min(c(mvnH[[i]]$breaks, mvtH[[i]]$breaks))
	  maxY[i] = max(c(mvnH[[i]]$counts, mvtH[[i]]$counts))
	}
	
	postscript("Output/test2h2-2.eps", width = 16, height = 14)
	par(mfrow = c(2,3))
	for(i in 1:2){
	  plot(mvtH[[i]], xlim = c(minX[i], maxX[i]), ylim = c(0, maxY[i]), border = "steelblue", xlab = "", main = cnames[i])
	  plot(mvnH[[i]], border = "tomato", add = TRUE)
	  abline(v = forcM[i], col = "black", lwd=1)
	  legend("topleft", c("MVT", "MVN"), col = c("steelblue", "tomato"), bty="n", lty=c(1,1))
	  legend("topright", c("Forecast"), col = c("black"), bty="n", lty=1)
	}
	dev.off()	
	




