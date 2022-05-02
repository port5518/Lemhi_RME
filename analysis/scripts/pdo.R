# Author: Kevin See
# Purpose: PDO data for use with survival from LLRTP to LGR
# Created: 6/14/21
# Last Modified: 5/2/22
# Notes:

rm(list = ls())

#---------------------
# load packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(here)
library(corrr)

#---------------------
# read in PDO data from NOAA
pdo_df = read_table("https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/pdo.long.data",
                    skip = 1,
                    col_names = c("year",
                                  "Jan",
                                  "Feb",
                                  "Mar",
                                  "Apr",
                                  "May",
                                  "Jun",
                                  "Jul",
                                  "Aug",
                                  "Sep",
                                  "Oct",
                                  "Nov",
                                  "Dec"),
                    n_max = 122) %>%
  mutate(across(everything(),
                as.numeric)) %>%
  pivot_longer(cols = -year,
               names_to = "month",
               values_to = "pdo") %>%
  mutate(date = ym(paste(year, month)) + days(14)) %>%
  mutate(type = if_else(pdo > 0, "Pos", "Neg"))

# for plotting, interpolate between points when the PDO crosses 0
pdo_df2 = pdo_df %>%
  mutate(grp = "orig") %>%
  bind_rows(pdo_df %>%
              mutate(lead_pdo = lead(pdo),
                     lead_date = lead(date),
                     sign_pdo = sign(pdo),
                     sign_lead = sign(lead_pdo),
                     switch = sign_pdo != sign_lead) %>%
              mutate(lin_extrap = map2_dbl(pdo, lead_pdo,
                                           .f = function(a, b) {
                                             tibble(x = c(0, 1),
                                                    y = c(a, b)) %>%
                                               lm(x ~ y,
                                                  data = .) %>%
                                               coef() %>%
                                               magrittr::extract(1)
                                           }),
                     lin_extrap = if_else(switch, lin_extrap, NA_real_)) %>%
              filter(!is.na(lin_extrap)) %>%
              mutate(x_date = difftime(lead_date, date, units = "days"),
                     x_date = as.numeric(x_date),
                     x_date = x_date * lin_extrap,
                     x_date = date + days(round(x_date))) %>%
              mutate(y_pdo = 0) %>%
              select(year, month,
                     pdo = y_pdo,
                     date = x_date) %>%
              mutate(grp = "new")) %>%
  arrange(date) %>%
  tidyr::fill(type, .direction = "up") %>%
  # filter errant value
  filter(pdo < 5)

#---------------------
# plot
yr = min(pdo_df2$year)
pdo_p = pdo_df2 %>%
  filter(year >= yr) %>%
  filter(pdo < 5) %>%
  ggplot(aes(x = date,
             y = pdo)) +
  geom_hline(yintercept = 0) +
  geom_line(color = "black") +
  # geom_line(aes(color = type,
  #               group = 1)) +
  geom_area(data = pdo_df2 %>%
              filter(year >= yr) %>%
              filter(pdo >= 0),
            aes(fill = "Pos")) +
  geom_area(data = pdo_df2 %>%
              filter(year >= yr) %>%
              filter(pdo <= 0),
            aes(fill = "Neg")) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1",
                    name = "Sign") +
  scale_color_brewer(palette = "Set1",
                     name = "Sign") +
  labs(y = "PDO",
       x = "Date")
pdo_p

# save plot
ggsave(file = here("analysis/figures/pdo.png"),
       pdo_p)

#---------------------
# read in components of PDO
for(nm in list.files(here("analysis/data/raw_data/PDO"))) {
  df = read_csv(here("analysis/data/raw_data/PDO", nm))
  names(df)[2] = str_split(names(df)[2], "https", simplify = T)[1] %>%
    str_remove("^PDO from ") %>%
    str_squish()
  if(nm == list.files(here("analysis/data/raw_data/PDO"))[1]) {
    pdo_comp_df = df
  } else {
    pdo_comp_df %<>%
      full_join(df)
  }
  rm(df)
}

#---------------------
# look at correlation between metrics
pdo_comp_df %>%
  select(-Date) %>%
  corrr::correlate() %>%
  corrr::shave() %>%
  corrr::stretch() %>%
  filter(!is.na(r)) %>%
  arrange(desc(r))

#---------------------
# plot time series of PDO components
pdo_comp_p = pdo_comp_df %>%
  pivot_longer(cols = -Date,
               names_to = "metric",
               values_to = "value") %>%
  filter(value != -9999) %>%
  ggplot(aes(x = Date,
             y = value,
             color = metric)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  geom_line()
pdo_comp_p
