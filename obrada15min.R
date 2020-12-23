
#---------SREDJIVANJE GENERATION fajla za 2017 i 2018 godinu
# Procitam fajl koji u excelu ima moving average za houly i data podatke
# onda citam svaki 4 i svaki 96. za satne, tj. dnevne srednje vrednosti respektivno
# onda ih spojim......malo je seljacki ali to je jedini nacin da izvucem podatke koji su na sajtu ENTSO-E Transparency
# dati sa 15min sampling time za HUPX i satni za SEEPEX

daily2017DB <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/ENTSOE/Ovo koristi/Hungary/GenerationForecast2017.xlsx", sheet = 1)
daily2018DB <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/ENTSOE/Ovo koristi/Hungary/GenerationForecast2018.xlsx", sheet = 1)

#dailyload2017f <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/ENTSOE/Ovo koristi/Hungary/LoadForecast2017.xlsx", sheet = 1)
#dailyload2018f <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/ENTSOE/Ovo koristi/Hungary/LoadForecast2018.xlsx", sheet = 1)
#dailyload2017a <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/ENTSOE/Ovo koristi/Hungary/LoadActual2017.xlsx", sheet = 1)
#dailyload2018a <- read_excel("/Users/marijanrancic/Documents/IMQF/Master/Data/ENTSOE/Ovo koristi/Hungary/LoadActual2018.xlsx", sheet = 1)


#----------------GENERATION------------------------------------------------------
# ---- citam samo 4. i 96. case
hourly2017 <- daily2017DB %>% filter(row_number() %% 4 == 1)
daily2017 <- daily2017DB %>% filter(row_number() %% 96 == 1)
# ---- citam samo 4. i 96. case
hourly2018 <- daily2018DB %>% filter(row_number() %% 4 == 1)
daily2018 <- daily2018DB %>% filter(row_number() %% 96 == 1)

# --sad obradjujem samo dnevne vrednosti --------
#sve ovo je malo seljacki ali nemam vremena da pisem funkcije za dva seta

# ------GENERATION DAILY-------------------------
avgrawdailygen2017 <- select(avgdailygen2017, c(-1,-2,-3))
avgrawdailygen2018 <- select(avgdailygen2018, c(-1,-2,-3))
# GEN HOURLY
avgrawhourlygen2017 <- select(avghourlygen2017, c(-1,-2,-4))
avgrawhourlygen2018 <- select(avghourlygen2018, c(-1,-2,-4))




#-------------------SADA UKLAPAM---------------------
# povezivanje godina 2017 i 2018

dailydata1 <- bind_rows(avgrawdailygen2017,avgrawdailygen2018)
daily_hupx_gen <- bind_cols(daily_hupx,dailydata1)


theme_set(theme_bw())  # pre-set the bw theme.
g <- ggplot(daily_hupx_gen, aes(x = daily, y = avg))
# Scatterplot
g + geom_point() + 
  geom_smooth(method="lm", se=F) +
  labs(y="Spot Price", 
       x="Generation", 
       title="Scatterplot with overlapping points")