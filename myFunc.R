show_moments <-  function(vektor) {
  library(moments)
  DS_mean = mean (vektor)
  DS_max = max (vektor)
  DS_min = min(vektor)
  DS_sd = sd(vektor)       
  DS_skew = skewness(vektor)  
  DS_kurt = kurtosis(vektor)
  descriptive_statistic <- c(DS_mean,DS_max,DS_min,DS_sd,DS_skew,DS_kurt)
  naslov <- c("Mean", "Max", "Min", "SD", "Skewness", "Kurtosis")
  des_stat <- tibble(descriptive_statistic,naslov)
  return(des_stat)
  invisible(vektor)
}

#---------------------------------
DBHour <- function(sat){
  library(dplyr)
  DBH <- DBTS %>% dplyr::filter(DHour==sat)
  return(DBH)
}

DBHourIN <- function(sat){
  library(dplyr)
  DBH <- DBTS %>% filter(year(date.time)==2017 | year(date.time)==2018) %>% dplyr::filter(DHour==sat)
  return(DBH)
}
DBHourOUT <- function(sat){
  library(dplyr)
  DBH <- DBTS %>% filter(year(date.time)==2019 ) %>% dplyr::filter(DHour==sat)
  return(DBH)
}
#----------------------------------
mojastatistika <- function(vektor){
  library(tseries)
  library(moments)
  LJB = Box.test(vektor, lag = 20, type = c("Ljung-Box"))
  LJB1= Box.test(vektor^2, lag = 20, type = c("Ljung-Box"))
  JBT= jarque.bera.test(vektor)
  DIK=adf.test(vektor)
  Statistika = c(LJB,LJB1,JBT,DIK)
  return(Statistika)
}

show_moments1 <-  function(vektor1) {
  library(moments)
  DS_mean1 = mean (vektor1)
  DS_max1 = max (vektor1)
  DS_min1 = min(vektor1)
  DS_var1 = var(vektor1)      
  DS_sd1 = sd(vektor1)       
  DS_skew1 = skewness(vektor1)  
  DS_kurt1 = kurtosis(vektor1)
  descriptive_statistic <- c(DS_mean1,DS_max1,DS_min1,DS_var1,DS_sd1,DS_skew1,DS_kurt1)
  return(descriptive_statistic)
  invisible(vektor1)
}

n.colmeans = function(df, n = 24){
  aggregate(x = df,
            by = list(gl(ceiling(nrow(df)/n), n)[1:nrow(df)]),
            FUN = mean)
} 

#-----------------RETURNS
# Provide dataframe and name of the column and it will add two columns arithmetic and log returns

add.ret = function(data,name){
  pars <- as.list(match.call()[-1])
  vektor = data[,as.character(pars$name)]
  temp_val1<- paste("ret.art.",as.character(pars$name), sep = "")
  kurac <- data
data <- data %>% mutate(temp_val1 = c(0, diff(vektor,lag=1))/vektor)


}

show_moments1(prices$HU)



