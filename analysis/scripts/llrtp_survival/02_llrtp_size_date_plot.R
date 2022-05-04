# Author: Mike Ackerman
# Purpose: Re-create Figure 2 from Copeland et al. 2014
# Created: 2/23/22
# Last Modified: 4/25/22
# Notes:

#-------------------
# clear environment
rm(list = ls())

#-------------------
# load packages
library(tidyverse)
library(here)
library(janitor)
library(lubridate)
library(ggpubr)

#-------------------
# read in LLRTP tag data
lem_rst_tags = read_rds(here("analysis/data/derived_data/lemhi_rst_tags_cleaned.rds"))

#-------------------
# explore a bit
tabyl(lem_rst_tags, brood_year, site_name)
tabyl(lem_rst_tags, brood_year, source)
tabyl(lem_rst_tags, brood_year, emig_stage)
tabyl(lem_rst_tags, brood_year, strategy)

lem_df = lem_rst_tags %>%
  select(tag_code,
         brood_year,
         date,
         length,
         emig_stage,
         strategy)

tabyl(lem_df, strategy, emig_stage)

#-------------------
# plot by strategy (DSR, NRR)
strategy_p = lem_df %>%
  ggplot(aes(x = as.Date(yday(date), "1900-01-01"),
             y = length)) +
  geom_point(aes(
    #shape = factor(strategy),
    color = factor(strategy))) +
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

#-------------------
# plot by life stage (smolt, parr, presmolt)
ls_p = lem_df %>%
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
ls_p

# arrange the plots
rst_p = ggarrange(strategy_p,
                  ls_p,
                  ncol = 1)
rst_p

ggsave(file = here("analysis/figures/llrtp_size_timing.png"),
       rst_p)

# END SCRIPT
