dji_raw <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/DJI2017-2019.xlsx")
dji.ts <- dji_raw %>%
  mutate(datum= ymd(datum)) %>%
  as_tsibble(index =datum)

peak <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/FRM/DB_FRM.xlsx", sheet = 1)
peak.ts <- peak %>%
  mutate(datum= ymd(datum)) %>%
  as_tsibble(index =datum)

ts_com <- left_join(dji.ts , peak.ts ,
                    by = c("datum" = "datum")) %>% filter(!is.na(Close))