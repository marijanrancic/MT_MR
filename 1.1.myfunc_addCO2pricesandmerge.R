library(data.table)
EUCO2 <- EUCO2 %>%
  mutate(Date= ymd(Date)) %>%
  as_tsibble(index =Date)


EUCO2 <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/eua-price-2.xlsx", sheet = 1) #stavi ih na WWW
EUCO2 <- EUCO2 %>%
  mutate(Date= ymd(Date)) %>%
  as_tsibble(index =Date)
DBTS.SPREAD <- left_join(DBTS.SPREAD  , EUCO2  ,
                         by = c("datum" = "Date")) 


# This is way around to fill NA values with last value (so for Saturday and Sunday we have Friday values)
d <- data.table(x=DBTS.SPREAD$datum, y=DBTS.SPREAD$Price)
setnafill(d, type = "locf", cols = "y") #
d1<- d %>%
  mutate(x= ymd(x)) %>%
  as_tsibble(index =x)
DBTS.SPREAD <- left_join(DBTS.SPREAD  , d1  ,
                         by = c("datum" = "x")) 

write_excel_csv(DBTS.SPREAD, file = "/Users/marijanrancic/Documents/IMQF/Master/Data/BAZA.csv")
