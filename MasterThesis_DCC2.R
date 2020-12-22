#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##
#################################################################################

#################################################################################
# DCC model
#################################################################################
# Fit Tests
# DCC under different specifications
# We have already defined X in IMPORT DATA

	Dat = X
	uspec = ugarchspec(mean.model = list(armaOrder = c(2,1)), variance.model = list(garchOrder = c(3,1), model = "sGARCH"), 
			distribution.model = "norm")
	#Vodi racuna o ovom broju  multispec( replicate(3, uspec) ) (za HUPX-SEEPEX je 2 a kada koristim i OPCOM ili vise treba uzeti 3 ili koliko vise)
	spec1 = dccspec(uspec = multispec( replicate(2, uspec) ), dccOrder = c(1,1),  distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se=FALSE))
	
	specx1 = ugarchspec(mean.model = list(armaOrder = c(2,1)), variance.model = list(garchOrder = c(3,1), model = "sGARCH",
					variance.targeting=TRUE), 
			distribution.model = "norm")
	specx2 = ugarchspec(mean.model = list(armaOrder = c(0,1)), variance.model = list(garchOrder = c(1,1), model = "eGARCH",
					variance.targeting=F), 
			distribution.model = "norm")
	specx3 = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(garchOrder = c(1,1), model = "apARCH"), 
			distribution.model = "norm")
	spec2 = dccspec(uspec = multispec( list(specx1, specx2) ), dccOrder = c(1,1),  distribution = "mvnorm")
	fit2 = dccfit(spec2, data = Dat, fit.control = list(eval.se=T))
	
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0),  include.mean = FALSE), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
			distribution.model = "norm")
	spec3 = dccspec(uspec = multispec( replicate(2, uspec) ), VAR = TRUE, dccOrder = c(1,1),  distribution = "mvnorm")
	fit3 = dccfit(spec3, data = Dat, fit.control = list(eval.se=FALSE))
	
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0),  include.mean = FALSE), variance.model = list(garchOrder = c(1,1), model = "sGARCH"), 
			distribution.model = "norm")
	spec4 = dccspec(uspec = multispec( replicate(2, uspec) ), VAR = TRUE, lag = 3, dccOrder = c(1,1),  distribution = "mvnorm")
	fit4 = dccfit(spec4, data = Dat, fit.control = list(eval.se=FALSE))
	

	spec5 = dccspec(uspec = multispec( list(specx1, specx2) ), dccOrder = c(1,1), model="aDCC", distribution = "mvnorm")
	fit5 = dccfit(spec5, data = Dat, fit.control = list(eval.se=FALSE))
	
	np = c(
			length(fit1@mfit$matcoef[,1]),
			length(fit2@mfit$matcoef[,1]),
			length(fit3@mfit$matcoef[,1]),
			length(fit4@mfit$matcoef[,1]),
			length(fit5@mfit$matcoef[,1])) +  (3^2 - 3)/2
	aicmod = c(
			rugarch:::.information.test(fit1@mfit$llh, nObs = fit1@model$modeldata$T, 
					nPars = np[1])[[1]],
			rugarch:::.information.test(fit2@mfit$llh, nObs = fit2@model$modeldata$T, 
					nPars = np[2])[[1]],
			rugarch:::.information.test(fit3@mfit$llh, nObs = fit3@model$modeldata$T, 
					nPars = np[3])[[1]],
			rugarch:::.information.test(fit4@mfit$llh, nObs = fit4@model$modeldata$T, 
					nPars = np[4])[[1]],
			rugarch:::.information.test(fit5@mfit$llh, nObs = fit5@model$modeldata$T, 
					nPars = np[5])[[1]])
	tmp = data.frame(n.pars = np, AIC = aicmod)
	rownames(tmp) = paste("Model", 1:5, sep = "")
	
	options(width = 120)
	zz <- file("Output/test2a-AIC4Models.txt", open="wt")
	sink(zz)	
	print(tmp)
	sink(type="message")
	sink()
	close(zz)
	
	rc1 = rcor(fit1)
	rc2 = rcor(fit2)
	rc3 = rcor(fit3)
	rc4 = rcor(fit4)
	rc5 = rcor(fit5)
	
	postscript("Output/test2a.eps", width = 10, height = 8)
	plot(rc2[1,2, ], type = "l")
	lines(rc5[1,2, ], col = 2)
	legend("topleft", legend = c("DCC", "aDCC"), col= 1:2, fill = 1:2)
	dev.off()
	
#-------------------------------------------------------------------------------
	# We can either let dccfit fit everything (VAR and garch), else we can
	# supply the VAR.fit and garch (fit) objects which are pre-estimated (so
	# that we do not have to estimate them everytime we change 2-stage details).
	
	uspec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
			variance.model = list(garchOrder = c(1,1), model = "eGARCH"), 
			distribution.model = "norm")
	
	mspec = multispec( replicate(2, uspec) )
	
	VAR.fit = varxfit(Dat, p = 2, postpad = "constant")
	resx = VAR.fit$xresiduals
	multf = multifit(mspec, data = resx)
	
	spec1 = dccspec(uspec = mspec, dccOrder = c(1,1), VAR = TRUE, lag = 2, 
			 distribution = "mvnorm")
	fit1 = dccfit(spec1, data = Dat, fit.control = list(eval.se=TRUE), 
			VAR.fit = VAR.fit, fit = multf)
	
	# check that we have the same answers when not passing pre-fitted objects:
	# (will not work for the robust version which has simulation embedded).
	fitcheck = dccfit(spec1, data = Dat, fit.control = list(eval.se=FALSE))
	
	options(width = 120)
	zz <- file("Output/test2c1.txt", open="wt")
	sink(zz)
	print( all.equal( head(fitted(fitcheck)), head(fitted(fit1)) ) )
	print( all.equal( last(rcov(fitcheck), 1)[,,1], last(rcov(fitcheck), 1)[,,1] ) )
	print( all.equal( first(rcov(fitcheck), 1)[,,1], first(rcov(fitcheck), 1)[,,1] ) )
	sink(type="message")
	sink()
	close(zz)
	
	spec1a = dccspec(uspec = mspec, dccOrder = c(1,1), VAR = TRUE, lag = 2, 
			model="aDCC", distribution = "mvnorm")
	fit1a = dccfit(spec1a, data = Dat, fit.control = list(eval.se=TRUE), 
			VAR.fit = VAR.fit, fit = multf)
	
	spec2 = dccspec(uspec = mspec, dccOrder = c(1,1), VAR = TRUE, lag = 2,  
			distribution = "mvlaplace")
	fit2 = dccfit(spec2, data = Dat, fit.control = list(eval.se=TRUE), 
			VAR.fit = VAR.fit, fit = multf)
	
	spec2a = dccspec(uspec = mspec, dccOrder = c(1,1), VAR = TRUE, lag = 2, model="aDCC", 
			distribution = "mvlaplace")
	fit2a = dccfit(spec2a, data = Dat, fit.control = list(eval.se=TRUE), 
			VAR.fit = VAR.fit, fit = multf)
	
	# The proper way to estimate the DCC-Student Model in a 2-stage setup
	# here we only pass the VAR.fit since the garch model is different
	uspec = ugarchspec(
			mean.model = list(armaOrder = c(1,1)), 
			variance.model = list(garchOrder = c(1,1), model = "eGARCH"), 
			distribution.model = "std", fixed.pars = list(shape=5))
	
	spec3 = dccspec(uspec = multispec( replicate(2, uspec) ), dccOrder = c(1,1), 
			VAR = TRUE, lag = 2,  distribution = "mvt", fixed.pars = list(mshape = 5))
	fit3 = dccfit(spec3, data = Dat, solver = "solnp", 
			VAR.fit = VAR.fit, fit.control = list(eval.se=TRUE), parallel = parallel, 
			parallel.control = parallel.control)
	
	spec3a = dccspec(uspec = multispec( replicate(2, uspec) ), dccOrder = c(1,1), 
			VAR = TRUE, lag = 2, model="aDCC", distribution = "mvt", fixed.pars = list(mshape = 5))
	fit3a = dccfit(spec3a, data = Dat, solver = "solnp", fit.control = list(eval.se=TRUE),
			VAR.fit = VAR.fit)
	
	
	# create a table for the DCC coefficients
	dccf = matrix(NA, ncol = 6, nrow = 4)
	dccp = matrix(NA, ncol = 6, nrow = 3)
	dccf[1:2,1] = coef(fit1, "dcc")
	dccf[4, 1] = likelihood(fit1)
	dccp[1:2, 1] = c(fit1@mfit$matcoef["[Joint]dcca1", 4], fit1@mfit$matcoef["[Joint]dccb1", 4])
	dccf[1:3,2] = coef(fit1a, "dcc")
	dccf[4, 2]  = likelihood(fit1a)
	dccp[1:3, 2] = c(fit1a@mfit$matcoef["[Joint]dcca1", 4], 
			fit1a@mfit$matcoef["[Joint]dccb1", 4],
			fit1a@mfit$matcoef["[Joint]dccg1", 4])
	
	dccf[1:2,3] = coef(fit2, "dcc")
	dccf[4,3] = likelihood(fit2)
	dccp[1:2, 3] = c(fit2@mfit$matcoef["[Joint]dcca1", 4], fit2@mfit$matcoef["[Joint]dccb1", 4])
	
	dccf[1:3,4] = coef(fit2a, "dcc")
	dccf[4,4] = likelihood(fit2a)
	dccp[1:3, 4] = c(fit2a@mfit$matcoef["[Joint]dcca1", 4], 
			fit2a@mfit$matcoef["[Joint]dccb1", 4],
			fit2a@mfit$matcoef["[Joint]dccg1", 4])
	
	dccf[1:2,5] = coef(fit3, "dcc")[1:2]
	dccf[4,5] = likelihood(fit3)
	dccp[1:2, 5] = c(fit3@mfit$matcoef["[Joint]dcca1", 4], fit3@mfit$matcoef["[Joint]dccb1", 4])
	
	dccf[1:3,6] = coef(fit3a, "dcc")[1:3]
	dccf[4,6] = likelihood(fit3a)
	dccp[1:3, 6] = c(fit3a@mfit$matcoef["[Joint]dcca1", 4], 
			fit3a@mfit$matcoef["[Joint]dccb1", 4],
			fit3a@mfit$matcoef["[Joint]dccg1", 4])
	
	dccfdf = as.data.frame(dccf)
	starsdf = apply(dccp, 2, FUN = function(x) rugarch:::.stars(x, levels = c(0.01, 0.05, 0.1)))
	for(i in 1:6){
		for(j in 1:3){
			if(!is.na(dccf[j,i])) 
				dccfdf[j, i] = paste(as.character(round(dccf[j, i],5)), starsdf[j,i],sep="")
			else
				dccfdf[j, i] = ""
		}
		dccfdf[4,i] = as.character(round(dccf[4,i], 2))
	}
	colnames(dccfdf) = c("DCC-MVN", "aDCC-MVN", "DCC-MVL", "aDCC-MVL", "DCC-T[5]", "aDCC-T[5]")
	rownames(dccfdf) = c("a", "b", "g", "LL")
	
	options(width = 120)
	zz <- file("Output/test2c1.txt", open="wt")
	sink(zz)
	print(dccfdf)
	sink(type="message")
	sink()
	close(zz)
	
	rc1 = rcor(fit1)
	rc2 = rcor(fit2)
	rc3 = rcor(fit3)
	rc1a = rcor(fit1a)
	rc2a = rcor(fit2a)
	rc3a = rcor(fit3a)
	
	postscript("Output/test2c1-2.eps", width = 10, height = 8)
	par(mfrow = c(2,1))
	plot(rc1[1,2, ], type = "l", main = "DCC")
	lines(rc2[1,2, ], col = 3)
	lines(rc3[1,2, ], col = 4)
	legend("bottomright", legend = c("MVNORM", "MVLAPLACE", "MVT[5]"), col= c(1,3,4), fill = c(1,3,4), bty="n")
	
	plot(rc1a[1,2, ], type = "l", main = "aDCC")
	lines(rc2a[1,2, ], col = 3)
	lines(rc3a[1,2, ], col = 4)
	legend("bottomright", legend = c("MVNORM", "MVLAPLACE", "MVT[5]"), col= c(1,3,4), fill = c(1,3,4), bty="n")
	
	dev.off()
	
	# weighted margins (elliptical distributions)
	# fitted returns mu of size equal to dataset#
	port1 = wmargin("mvnorm", weights = rep(1/2,2), Sigma = rcov(fit1), mean = fitted(fit1))
	qport1 = qdist("norm", 0.01, port1[,1], port1[,2], 0, 0, 0)
	
	
	port2 = wmargin("mvlaplace", weights = rep(1/2,2), Sigma = rcov(fit2), mean = fitted(fit2))
	# REM: GED (with shape = 1) == Laplace
	# qport2 = apply(port2, 1, FUN = function(x) qdist("ged", 0.01, x[1], x[2], 0, 0, 1))
	qport2 = qdist("ged", 0.01, port2[,1], port2[,2], 0, 0, 1)
	
	port3 = wmargin("mvt", weights = rep(1/2,2), Sigma = rcov(fit3), mean = fitted(fit3), shape = 5, skew =0 )
	qport3 = qdist("std", 0.01, port3[,1], port3[,2], 0, 0, port3[,4])
	actual = apply(Dat, 1, "mean")
	
	postscript("Output/test2c2-1.eps", width = 10, height = 8)
	VaRplot(0.01, xts::as.xts(actual), xts::xts(qport1, as.POSIXct(rownames(Dat))))
	lines(xts::xts(qport3, as.POSIXct(rownames(Dat))), col = 3)
	lines(xts::xts(qport2, as.POSIXct(rownames(Dat))), col = 4)
	legend("topright", legend = c("MVNORM", "MVLAPLACE", "MVT[5]"), bty="n", col= c(1,3,4), fill = c(1,3,4))
	dev.off()
	
	options(width = 120)
	zz <- file("Output/test2c2.txt", open="wt")
	sink(zz)
	rugarch:::.VaRreport("Port", "DCC", "MVNORM", 0.01, actual, qport1)
	rugarch:::.VaRreport("Port", "DCC", "MVLAPLACE", 0.01, actual, qport2)
	rugarch:::.VaRreport("Port", "DCC", "MVT[5]", 0.01, actual, qport3)
	sink(type="message")
	sink()
	close(zz)
	




