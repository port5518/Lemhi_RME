# Author: Kevin See and Mike Ackerman
# Purpose: grab PIT tags from Lemhi RST
# Created: May 27, 2021
# Last Modified: April 22, 2022

# Notes:

rm(list = ls())

#------------------------
# load packages
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(ggpubr)

#------------------------
# read in tagging details
# ran PTAGIS Tagging Details report, filtered on mark subbasin is Lemhi, capture method is screw trap, species is Chinook and mark year is 1986-2020.
lem_mrk_tags = read_csv(here("analysis/data/raw_data", "Lemhi_RST_Chnk_Tagging_Detail_1986-2020.csv")) %>%
  clean_names() %>%
  mutate(across(ends_with("mmddyyyy"),
                mdy)) %>%
  filter(mark_site_code_value %in% c("LEMHIR",
                                     "LLRTP"),
         release_site_code_value != "LEMHIW")

# ran PTAGIS Recapture Details report, filtered on mark subbasin is Lemhi, capture method is screw trap, species is Chinook and recap year is 1986-2020.
lem_recap_tags = read_csv(here("analysis/data/raw_data", "Lemhi_RST_Chnk_Recapture_Detail_1986-2020.csv")) %>%
  clean_names() %>%
  mutate(across(ends_with("mmddyyyy"),
                mdy)) %>%
  filter(recap_site_code_value %in% c("LEMHIR",
                                      "LLRTP"),
         recap_species_code == 1,
         recap_release_site_info_site_code != "LEMHIW")

# put them together
lem_tags_all = lem_mrk_tags %>%
  select(tag_code,
         site_name = mark_site_name,
         site_code = mark_site_code_value,
         release_site = release_site_code_value,
         release_rkm = release_site_rkm_value,
         date = mark_date_mmddyyyy,
         length = length_mm,
         weight = weight_g) %>%
  distinct() %>%
  mutate(source = "mark_details") %>%
  bind_rows(lem_recap_tags %>%
              select(tag_code,
                     site_name = recap_site_name,
                     site_code = recap_site_code_value,
                     release_site = recap_release_site_info_site_code,
                     release_rkm = recap_release_site_rkm_value,
                     date = recap_date_mmddyyyy,
                     length = recap_length_mm,
                     weight = recap_weight_g) %>%
              distinct() %>%
              mutate(source = "recap_details")) %>%
  arrange(tag_code, date) %>%
  # assign juveniles to a life stage, standard classification
  mutate(
    emig_stage = "nada",
    emig_stage = if_else(month(date) > 8, "Presmolt", emig_stage),
    emig_stage = if_else(month(date) < 7, "Smolt", emig_stage),
    emig_stage = if_else(between(month(date), 7, 8), "Parr", emig_stage)
  ) %>%
  # make some corrections based on size & date cutoffs
  mutate(
    emig_stage = if_else(yday(date) >  105 & yday(date) <= 120 & length <= 75,  "Parr", emig_stage),
    emig_stage = if_else(yday(date) >= 121 & yday(date) <= 135 & length <= 85,  "Parr", emig_stage),
    emig_stage = if_else(yday(date) >= 136 & yday(date) <= 151 & length <= 90,  "Parr", emig_stage),
    emig_stage = if_else(yday(date) >= 152 & yday(date) <= 158 & length <= 95,  "Parr", emig_stage),
    emig_stage = if_else(yday(date) >= 159 & yday(date) <= 166 & length <= 105, "Parr", emig_stage),
    emig_stage = if_else(yday(date) >= 167 & yday(date) <= 182 & length <= 115, "Parr", emig_stage),
  ) %>%
  mutate(brood_year = if_else(emig_stage == "Smolt",
                              year(date) - 2,
                              year(date) - 1)) %>%
  mutate(strategy = case_when(emig_stage == "Smolt" ~ "NRR",
                              emig_stage == "Parr" ~ "DSR",
                              emig_stage == "Presmolt" ~ "DSR"))

# plot life history strategies
strategy_p = lem_tags_all %>%
  select(tag_code,
         date,
         length,
         emig_stage,
         strategy,
         brood_year) %>%
  filter(length < 200 & length > 30) %>%
  ggplot(aes(x = as.Date(yday(date), "1900-01-01"),
             y = length)) +
  geom_point(aes(
    #shape = factor(strategy),
    color = factor(strategy))) +
  # scale_shape_manual(name = "Strategy",
  #                    values = c(18, 16)) +
  geom_vline(xintercept = as.Date("1900-03-01")) +
  geom_vline(xintercept = as.Date("1900-05-31")) +
  geom_vline(xintercept = as.Date("1900-09-15"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("1900-11-30"), linetype = "dashed") +
  scale_color_manual(name = "Strategy",
                     values = c("orangered2", "dodgerblue2")) +
  scale_x_date(date_breaks = "months",
               date_labels = "%b") +
  theme_classic() +
  labs(x = "Date Tagged",
       y = "Fork Length (mm)") +
  theme(axis.title.x = element_text(size = 11,
                                    color = "black"),
        axis.title.y = element_text(size = 11,
                                    color = "black"),
        axis.text.x = element_text(size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11,
                                   color = "black"),
        #legend.position = c(0.5, 0.9)
        legend.position = "top")
strategy_p

# plot by smolt, parr, presmolt
es_p = lem_tags_all %>%
  select(tag_code,
         date,
         length,
         emig_stage,
         strategy,
         brood_year) %>%
  filter(length < 200 & length > 30) %>%
  ggplot(aes(x = as.Date(yday(date), "1900-01-01"),
             y = length)) +
  geom_point(aes(color = factor(emig_stage))) +
  geom_vline(xintercept = as.Date("1900-03-01")) +
  geom_vline(xintercept = as.Date("1900-05-31")) +
  geom_vline(xintercept = as.Date("1900-09-15"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("1900-11-30"), linetype = "dashed") +
  scale_color_manual(name = "Emigrant Stage",
                     values = c("seagreen2", "slateblue2", "deepskyblue2")) +
  scale_x_date(date_breaks = "months",
               date_labels = "%b") +
  theme_classic() +
  labs(x = "Date Tagged",
       y = "Fork Length (mm)") +
  theme(axis.title.x = element_text(size = 11,
                                    color = "black"),
        axis.title.y = element_text(size = 11,
                                    color = "black"),
        axis.text.x = element_text(size = 11,
                                   color = "black"),
        axis.text.y = element_text(size = 11,
                                   color = "black"),
        legend.title = element_text(size = 11,
                                    color = "black"),
        legend.text = element_text(size = 11,
                                   color = "black"),
        legend.position = "top")
es_p

rst_p = ggarrange(strategy_p,
                  es_p,
                  ncol = 1)
rst_p

ggsave(file = here("analysis/figures/lem_rst_size_timing.png"),
       rst_p)

# deal with tags that have more than one brood year, or more than one emigration stage
# get the minimum and maximum detection dates for each tag
min_max_tags = lem_tags_all %>%
  group_by(tag_code) %>%
  mutate(min_det = min(date),
         max_det = max(date),
         time_diff = difftime(max_det, min_det, units = "days")) %>%
  mutate(across(time_diff,
                as.numeric)) %>%
  ungroup()

# for tags recaptured within a week, use the first capture date
lt_week_tags = min_max_tags %>%
  filter(time_diff <= 7) %>%
  group_by(tag_code) %>%
  mutate(across(c(length, weight),
                max,
                na.rm = T)) %>%
  filter(date == min_det) %>%
  slice(1) %>%
  ungroup() %>%
  select(any_of(names(lem_tags_all)))

# for tags recaptured longer than a week later, use the last capture date
gt_week_tags = min_max_tags %>%
  filter(time_diff > 7) %>%
  group_by(tag_code) %>%
  mutate(across(c(length, weight),
                max,
                na.rm = T)) %>%
  filter(date == max_det) %>%
  slice(1) %>%
  ungroup() %>%
  select(any_of(names(lem_tags_all)))

# put them together
lem_tags_fix = lt_week_tags %>%
  bind_rows(gt_week_tags)

lem_tags_fix %>%
  filter(is.na(length))
lem_tags_fix %>%
  filter(is.na(weight))

lem_tags_fix %>%
  filter(tag_code %in% tag_code[duplicated(tag_code)])

identical(n_distinct(lem_tags_all$tag_code), n_distinct(lem_tags_fix$tag_code))

# grab all tag codes and save as a .txt file to upload to PTAGIS
# lem_tags_all %>%
#   select(tag_code, brood_year) %>%
#   distinct() %>%
#   filter(tag_code %in% tag_code[duplicated(tag_code)]) %>%
#   left_join(lem_tags_all) %>%
#   arrange(tag_code, date) %>%
#   group_by(tag_code) %>%
#   arrange(date) %>%
#   mutate(id = 1:n()) %>%
#   ungroup() %>%
#   arrange(tag_code, id)
#
# lem_tags_all %>%
#   select(tag_code, brood_year) %>%
#   distinct() %>%
#   filter(tag_code %in% tag_code[duplicated(tag_code)]) %>%
#   left_join(lem_tags_all) %>%
#   group_by(tag_code, brood_year) %>%
#   filter(date == max(date)) %>%
#   ungroup() %>%
#   select(tag_code:date) %>%
#   distinct() %>%
#   pivot_wider(names_from = brood_year,
#               values_from = date,
#               names_sort = T)

# split by life-stage
lem_tags_fix %>%
  split(list(.$emig_stage)) %>%
  map(.f = function(x) {
    x %>%
      select(tag_code) %>%
      distinct() %>%
      write_delim(file = here('analysis/data/raw_data/tag_lists',
                              paste0('Lemhi_', unique(x$emig_stage), '_tags.txt')),
                  delim = "\n",
                  col_names = F)
  })

# split by brood year
lem_tags_fix %>%
  filter(brood_year > 2003) %>%
  split(list(.$brood_year)) %>%
  map(.f = function(x) {
    x %>%
      select(tag_code) %>%
      distinct() %>%
      write_delim(file = here('analysis/data/raw_data/tag_lists',
                              paste0('Lemhi_BY', unique(x$brood_year), '_tags.txt')),
                  delim = "\n",
                  col_names = F)
  })

# save details
write_rds(lem_tags_all,
          file = here('analysis/data/derived_data/lemhi_rst_tags.rds'))
write_rds(lem_tags_fix,
          file = here('analysis/data/derived_data/lemhi_rst_tags_cleaned.rds'))
