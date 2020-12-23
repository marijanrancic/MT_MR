
library(moments)
library(tseries)
library(MASS)
library(stats)

momenit = tibble()

for (i in 1:9){
  temp_val<- paste("H0",i, sep = "")
  temp_mom <- show_moments1(DBHourOUT(temp_val)$RS)
    for (j in 1:7){
      momenit[i,j] = temp_mom[j]
    }
}

# BEcause of the Hours labeling H01, H02....have to do it in two loops
for (i in 10:24){
  temp_val <- paste("H",i, sep = "")
  temp_mom <- show_moments1(DBHourOUT(temp_val)$RS)
  for (j in 1:7){
    momenit[i,j] = temp_mom[j]
  }
  
}
rm(temp_val,i,j,temp_mom)

write_csv(momenit,file = "/Users/marijanrancic/Documents/IMQF/Master/Master teza TEXT/Slike/momentiRSOUT.csv")

y <- DBHourIN("H12")$HU
Box.test(y, lag = 20, type = c("Ljung-Box"))
Box.test(y^2, lag = 20, type = c("Ljung-Box"))
