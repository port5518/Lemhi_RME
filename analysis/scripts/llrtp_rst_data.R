# Author: Kevin See and Mike Ackerman
# Purpose: prep screw trap data from Lower Lemhi RST (LLRTP)
# Created: 6/9/21
# Last Modified: 4/29/22
# Notes:

#---------------------------
# clear environment
rm(list = ls())

#---------------------------
# load packages
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
# library(dbplyr)
# library(readxl)

#---------------------------
# read in screw trap data and save to NAS
# data collection events
dce = read_csv(here("analysis/data/raw_data/fld_DataCollectionEvent.csv"),
               guess_max = 1e6) %>%
  filter(dceType == "Rotary Screw Trap") %>%
  # pull out DCEs related to the lower Lemhi RST
  filter(SiteName %in% c("LLRTP",
                         "L3AO",
                         "L3A"))

# fish observations
fish_obs = read_csv("S:/main/data/fish/lemhi_rst/fld_FishObservation.csv",
                    guess_max = 1e6) %>%
  mutate(PitTagCaptureType = recode(PitTagCaptureType,
                                    "Non-efficiency Recapture" = "Non-Efficiency Recapture")) %>%
  filter(dceID %in% unique(dce$dceID))

# what species are present?
tabyl(fish_obs, Species) %>%
  arrange(desc(n)) %>%
  adorn_pct_formatting()

# trap operations
llrtp = read_csv(here("analysis/data/raw_data/LLRTP Trap Operations.csv")) %>%
  clean_names()

#---------------------------
# merge them together
rst_df = dce %>%
  select(dceID, dceType, SurveyDateTime, SiteName, dceName,
         Equipment,
         starts_with("RiverLocation"),
         ends_with("udeDD"),
         TrapLocation,
         matches("Release")) %>%
  # let's grab juv Chinook
  inner_join(fish_obs %>%
               filter(Species == "Chinook",
                      FishLifeStage != "Adult")) %>%
  # remove dead fish
  filter(FishStatus == "Live") %>%
  # assign life stages based on trapping date
  mutate(life_stage = if_else(between(month(SurveyDateTime), 7, 8),
                              "Parr",
                              if_else(month(SurveyDateTime) > 8,
                                      "Presmolt",
                                      if_else(FishLifeStage == "Fry",
                                              "Fry",
                                              "Smolt")))) %>%
  # add brood year
  mutate(brood_year = if_else(life_stage == "Smolt",
                              year(SurveyDateTime) - 2,
                              year(SurveyDateTime) - 1))

# life stage by brood year
tabyl(rst_df, brood_year, life_stage)

# life stage by brood year, new tags versus non-effeciency recaptures
rst_df %>%
  filter(PitTagCaptureType != "Efficiency Recapture") %>%
  group_by(life_stage, brood_year, PitTagCaptureType) %>%
  summarise(across(FishCount,
                   sum),
            .groups = "drop") %>%
  pivot_wider(names_from = brood_year,
              values_from = FishCount,
              names_sort = T)

# marks and recaptures by brood year and life stage
rst_df %>%
  group_by(life_stage, brood_year) %>%
  summarise(n_tot = sum(FishCount),
            n_M = sum(FishCount[PitTagCaptureType %in% c("New Tag", "Non-Efficiency Recapture")], na.rm = T),
            n_R = sum(FishCount[PitTagCaptureType == "Efficiency Recapture"], na.rm = T)) %>%
  pivot_longer(cols = starts_with("n_"),
               names_to = "type",
               values_to = "count") %>%
  pivot_wider(names_from = brood_year,
              values_from = count,
              names_sort = T)

# write rst_df to NAS
write_csv(rst_df,
          "S:/main/data/fish/lemhi_rst/llrtp_rst_df.csv")

# END SCRIPT
