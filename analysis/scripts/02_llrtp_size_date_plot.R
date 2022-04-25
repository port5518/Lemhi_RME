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

#-------------------
# read in LLRTP tag data
lem_rst_tags = read_rds(here("analysis/data/derived_data/lemhi_rst_tags_cleaned.rds"))
