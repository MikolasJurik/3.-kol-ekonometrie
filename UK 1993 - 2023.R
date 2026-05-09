
library(tidyquant) 
library(dplyr)
library(tidyr)
library(lubridate)
library(mFilter)
library(ggplot2)

# 1. STAŽENÍ DAT Z FRED (1993 - 2023)
fred_symbols <- c(
  "GBRGDPDEFQISMEI", # 1. Deflátor HDP
  "NGDPRSAXDCGBQ",   # 2. Reálné HDP 
  "ULQEUL01GBQ661S", # 3. Nominální ULC 
  "IRLTLT01GBQ156N", # 4. Dlouhodobá úroková míra (10Y)
  "IR3TIB01GBQ156N"  # 4. Krátkodobá úroková míra (3M)
)

fred_data_long <- tq_get(fred_symbols, 
                         get = "economic.data",
                         from = "1993-01-01",
                         to = "2023-04-01")

# 2. PŘEVOD DAT A VÝPOČTY PROMĚNNÝCH
uk_data_final <- fred_data_long %>%
  select(date, symbol, price) %>%
  pivot_wider(names_from = symbol, values_from = price) %>%
  rename(
    gdp_deflator = GBRGDPDEFQISMEI,
    real_gdp = NGDPRSAXDCGBQ,
    nominal_ulc = ULQEUL01GBQ661S,
    rate_10y = IRLTLT01GBQ156N,
    rate_3m = IR3TIB01GBQ156N
  ) %>%
  arrange(date) %>%
  drop_na() %>% # Odstranění chybějících hodnot před výpočty
  mutate(
    # 1. Míra inflace (mezikvartální procentní změna deflátoru)
    inflation = (log(gdp_deflator) - log(lag(gdp_deflator))) * 100,
    
    # 2. Labor Income Share (Reálné ULC = Nominální ULC / Deflátor)
    real_ulc = (nominal_ulc / gdp_deflator) * 100,
    
    # 3. Úrokové rozpětí 
    interest_spread = rate_10y - rate_3m,
    
    # 4. Mzdová inflace 
    wage_inflation = (log(nominal_ulc) - log(lag(nominal_ulc))) * 100
  ) %>%
  drop_na() # Odstranění prvního řádku (NA kvůli funkci lag)

# ==============================================================================
# 2.1 STACIONARIZACE DAT A VIZUALIZACE



# Výpočet stacionárních řad (procentních odchylek)
uk_data_model <- uk_data_final %>%
  mutate(
    # Krok A: Logaritmování a násobení 100 pro získání procentního měřítka
    log_gdp = log(real_gdp) * 100,
    log_ulc = log(real_ulc) * 100
  )

# Krok B: Aplikace HP filtru
# Pro čtvrtletní makroekonomická data se standardně používá parametr lambda = 1600
hp_gdp <- hpfilter(uk_data_model$log_gdp, freq = 1600)
hp_ulc <- hpfilter(uk_data_model$log_ulc, freq = 1600)

# Krok C: Extrakce cyklické složky (cycle) zpět do našeho datasetu
uk_data_model <- uk_data_model %>%
  mutate(
    output_gap = hp_gdp$cycle, # Mezera výstupu (značeno jako x_t)
    ulc_gap = hp_ulc$cycle     # Odchylka reálných mezních nákladů (značeno jako s_t)
  )



# Graf 1: Mezera výstupu (Output Gap)
plot_output_gap <- ggplot(uk_data_model, aes(x = date, y = output_gap)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Mezera výstupu (Output Gap) - Velká Británie",
    subtitle = "Procentní odchylka reálného HDP od trendu (HP filtr, lambda = 1600)",
    x = "Rok",
    y = "Odchylka (%)"
  ) +
  theme_minimal()

# Graf 2: Odchylka podílu produktu z práce (Labor Income Share Gap)
plot_ulc_gap <- ggplot(uk_data_model, aes(x = date, y = ulc_gap)) +
  geom_line(color = "forestgreen", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Odchylka Labor Income Share - Velká Británie",
    subtitle = "Procentní odchylka reálných ULC od trendu (HP filtr, lambda = 1600)",
    x = "Rok",
    y = "Odchylka (%)"
  ) +
  theme_minimal()

# Vykreslení odchylek
print(plot_output_gap)
print(plot_ulc_gap)
