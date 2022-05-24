#Length Freq for radio tags.
#Nick Porter 5.24.22


rm(list = ls()) # clear the working environment


#---------------------------
# load packages
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(esquisse)
# library(dbplyr)
# library(readxl)

#---------------------------
# read in screw trap data and save to NAS
# data collection events
lemhi_rst_tags_cleaned <- readRDS("~/Documents/GitHub/Lemhi_RME/analysis/data/derived_data/lemhi_rst_tags_cleaned.rds")

sum_by_yr=lemhi_rst_tags_cleaned %>%
  group_by(brood_year,emig_stage) %>%
  summarise(n = n(), "Min length" = min(length), "Mean length"= mean(length), "Max length" = max(length))

glimpse(lemhi_rst_tags_cleaned)


esquisser(lemhi_rst_tags_cleaned)
# histogram for all captured and PIT tagged fish
(yearly_len_freq = ggplot(lemhi_rst_tags_cleaned) +
  aes(x = length) +
  geom_histogram(bins = 20L) +
  theme_classic() +
  facet_wrap(~brood_year) +
  labs(x = "Length (mm)", y = "Frequency"))

setwd("~/Documents/GitHub/Lemhi_RME/analysis/figures")

ggsave('RT_all_years_length_freq.png',
       yearly_len_freq,
       width = 8,
       height = 8)


sum_by_yr=tagged_fish %>%
  group_by(season) %>%
  summarise(n= n())

min(tagged_fish$length)


