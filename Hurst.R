
data<-DBTS

N<-as.numeric(dim(data)[1])
GBPEUR<-data$HU
SEKEUR<-data$RS
CADEUR<-data$WGenRO
hurst_vol<-function(X)
{
  returns<-diff(X)/X[-length(X)]
  D=data$datum[2:length(data$datum)]
  test=c(1,2,4,8,16,32)
  #this vector will be used to get to the number of subdivision at any time
  adjusted_returns<-matrix(0,nrow=129,ncol=N)
  #i'll build a matrix of adjusted returns
  #I initialise this matrix with zeros because i will sum those returns
  stdev<-c(0)
  moyenne<-matrix(0,nrow=6,ncol=62)
  #i build a matrix of means of returns
  for(i in seq(1,length(test))){
    splitlength=length(returns)/test[i]
    #splitlength is the number of elements in a subdivision
    for(j in seq(1,test[i]))
    {
      temp<-returns[((j-1)*splitlength+1):(j*splitlength)]
      #temp is a temporary vector filled with the returns of a subdivision
      moyenne[i,j]=mean(temp)
      adjusted_returns[test[i]:sum(test[1:i]),1:splitlength]=temp-moyenne[i,j]
      stdev[test[i]:sum(test[1:i])]<-sd(temp)
      #i compute subdivision standard deviation for the rescaled range calculation
    }
  }
  deviate_series<-matrix(0,nrow=129,ncol=N)
  #I build a deviation series of the adjusted returns
  for(i in 1:129)
  {
    deviate_series[i,]<-cumsum(adjusted_returns[i,])
  }
  widest_difference<-c(0)
  #I compute the widest difference in all the deviation series
  for(i in 1:63)
  {
    widest_difference[i]<-max(deviate_series[i,])-min(deviate_series[i,])
  }
  rescaled_range<-c(0)
  #i use the widest difference and the standard deviation of each subdivision to compute rescaled ranges
  for(i in 1:63)
  {
    rescaled_range[i]<-widest_difference[i]/stdev[i]
  }
  rescaled_range_m<-c(0)
  size<-c(0)
  #i compute the rescaled range's means and I create a subdivision's size vector
  for(i in 1:6)
  {
    rescaled_range_m[i]<-mean(rescaled_range[test[i]:sum(test[1:i])])
    size[i]<-N/test[i]
  }
  #i put some log on the two vectors i just compute
  rescaled_range_m_l<-log10(rescaled_range_m)
  size_l<-log10(size)
  plot(size_l,rescaled_range_m_l,type='l',xlab='log10(size of subdivision)',ylab='log10(meaned rescaled range)',main='logarithmically adjusted rescaled range \n function of the size of subdivision',sub='the slope is the Hurst exponent of the time series')
  linear_regression<-lm(rescaled_range_m_l~size_l)
  #i make a linear regression to get the Hurst exponent
  Hurst<-as.numeric(linear_regression$coefficients[2])
  
  vol=sd(X)
  
  annualized_vol=vol*(length(X)**(Hurst))
  cat('Hurst exponent : ',Hurst)
  cat(sep='\n')
  cat('Annualized volatility : ',vol*sqrt(length(X)),'%')
  cat(sep='\n')
  cat('Annualized volatility using Hurst exponent: ',annualized_vol,'%')
}

hurst_vol(X=GBPEUR)
hurst_vol(CADEUR)
hurst_vol(SEKEUR)
