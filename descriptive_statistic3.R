
library(moments)
library(tseries)
library(MASS)
library(stats)
library(corrplot)

momCOR = tibble()
for (i in 1:9){
  temp_val<- paste("H0",i, sep = "")
  temp_mom <- cor(DBHourIN(temp_val)$HU,DBHourIN(temp_val)$RS)
  momCOR[i,1] = temp_mom  
}

# BEcause of the Hours labeling H01, H02....have to do it in two loops
for (i in 10:24){
  temp_val <- paste("H",i, sep = "")
  temp_mom <- show_moments1(DBHourOUT(temp_val)$RS)
  temp_mom <- cor(DBHourIN(temp_val)$HU,DBHourIN(temp_val)$RS)
  momCOR[i,1] = temp_mom 
  
}
rm(temp_val,i,j,temp_mom)
write_csv(momCOR,file = "/Users/marijanrancic/Documents/IMQF/Master/Master teza TEXT/Slike/momCOR.csv")

y <- DBTS.SPREAD$base.spread
Box.test(y, lag = 20, type = c("Ljung-Box"))
Box.test(y^2, lag = 20, type = c("Ljung-Box"))
jarque.bera.test(y)
jarque.test(y)
adf.test(y)
pp.test(y)
kpss.test(y, null="Trend")


# M je cor ()
corrplot(M, type="upper", order="hclust", 
         p.mat = p.mat, sig.level = 0.01)


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(M, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
)


library(ggplot2)
library(ggcorrplot)
ggcorrplot(M, hc.order = TRUE, type = "lower",
           lab = TRUE)

ggcorrplot(corr, hc.order = TRUE,
           type = "lower", p.mat = p.mat)
library(GGally)

forecastHUPX = data.frame(testing$HU,forecast_arima_HU,forecast_prophet_HU$yhat,forecast_gbm_HU, forecast_rf_HU)


#----------------




