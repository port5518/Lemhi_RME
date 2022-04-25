# Author: Kevin See
# Purpose: Survival from LLRTP to LGR
# Created: 6/14/21
# Last Modified: 4/25/22
# Notes: Reviewed by Mike A in April 2022

# clear environment
rm(list = ls())

#-----------------------------------------------------------------
# load packages
library(tidyverse)
library(here)
# remotes::install_github("mackerman44/PITcleanr@develop")
library(PITcleanr)
library(magrittr)
library(marked)
library(lubridate)
library(janitor)

theme_set(theme_bw())

#-----------------------------------------------------------------
# read in info about tags at lower Lemhi RST
lem_tags = read_rds(here('analysis/data/derived_data/lemhi_rst_tags_cleaned.rds'))

#-----------------------------------------------------------------
# separate juvenile and adult detection nodes at each dam
#-----------------------------------------------------------------

# build a configuration file that combines a few sites
# all detections at MRR site LGR are also seen at GRJ, so we're ignoring LGR as a site
config_file = buildConfig() %>%
  mutate(node = site_code) %>%
  # mutate(node = if_else(site_code == "LEMHIR",
  #                       "LLRTP",
  #                       node)) %>%
  # combine 2 sites at Lower Granite
  mutate(node = if_else(site_code %in% c("GRS"),
                        "GRJ",
                        node),
         node = if_else(site_code %in% c("LGRLDR"),
                        "GRA",
                        node),
         node = if_else(site_code %in% c("IHA",
                                         "IHR"),
                        "ICH",
                        node),
         node = if_else(site_code %in% c("MC1",
                                         "MC2"),
                        "MCN",
                        node),
         node = if_else(site_code %in% c("JDALD2",
                                         "JDALD1"),
                        "JO1",
                        node),
         node = if_else(site_code %in% c("TD2",
                                         "TDA"),
                        "TD1",
                        node),
         # combine all juvenile sites at Bonneville
         node = if_else(site_code %in% c("B2J",
                                         "BCC",
                                         "B1J",
                                         "BVX"),
                        "BOJ",
                        node),
         # combine all adult sites at Bonneville
         node = if_else(site_code %in% c("BO1",
                                         "BO2",
                                         "BO3",
                                         "BO4",
                                         "BONAFF"),
                        "BON",
                        node))

# create parent-child table
parent_child = tribble(~parent, ~child,
                       "LLRTP", "LLR",
                       "LLR", "GRJ",
                       "GRJ", "GOJ",
                       "GOJ", "LMJ",
                       "LMJ", "ICH",
                       "ICH", "MCJ",
                       "MCJ", "JDJ",
                       # "JDJ", "TD1",
                       # "TD1", "BOJ",
                       "JDJ", "BOJ",
                       "BOJ", "BON",
                       "BON", "TD1",
                       "TD1", "JO1",
                       "JO1", "MCN",
                       "MCN", "ICH",
                       "ICH", "LMA",
                       "LMA", "GOA",
                       "GOA", "GRA",
                       "GRA", "USE")

plotNodes(parent_child = parent_child)

#-----------------------------------------------------------------
# simpler version of parent_child
# parent_child = tribble(~parent, ~child,
#                        "LLRTP", "GRJ",
#                        "GRJ", "BOJ",
#                        "BOJ", "BON",
#                        "BON", "GRA",
#                        "GRA", "USE")
#
# plotNodes(parent_child = parent_child)

#-----------------------------------------------------------------
# get observation data from all years
obs_df = tibble(brood_year = 2004:2019) %>%
  mutate(ptagis_raw = map(brood_year,
                          .f = function(byr) {
                            readCTH(here("analysis/data/raw_data/PTAGIS",
                                         paste0("Lemhi_Chnk_Juv_BY", byr, ".csv")))
                          })) %>%
  # compress PTAGIS detections
  mutate(comp = map(ptagis_raw,
                    .f = function(x) {
                      compress(x,
                               configuration = config_file,
                               max_minutes = 60 * 24 * 10,
                               units = "days",
                               ignore_event_vs_release = T)
                    })) %>%
  # reset all capture histories to start at LLRTP
  mutate(llrtp_start = map(comp,
                           .f = function(comp_obs) {
                             comp_obs %>%
                               left_join(lem_tags %>%
                                           select(tag_code, start_date = date, emig_stage),
                                         by = "tag_code") %>%
                               # left_join(comp_obs %>%
                               #             filter(node %in% c("LLRTP", "LEMHIR"),
                               #                    event_type_name %in% c("Mark", "Recapture")) %>%
                               #             group_by(tag_code) %>%
                               #             filter(min_det == max(min_det)) %>%
                               #             slice(1) %>%
                               #             summarise(start_date = min_det,
                               #                       .groups = "drop"),
                               #           by = "tag_code") %>%
                               filter(min_det >= start_date | is.na(start_date)) %>%
                               group_by(tag_code) %>%
                               mutate(slot = slot - min(slot) + 1) %>%
                               ungroup() %>%
                               mutate(node = if_else(slot == 1,
                                                     "LLRTP",
                                                     node)) %>%
                               mutate(llrtp_time = difftime(min_det, start_date, units = "days")) %>%
                               mutate(across(llrtp_time,
                                             as.numeric)) %>%
                               mutate(across(emig_stage,
                                             factor,
                                             levels = c("Parr",
                                                        "Presmolt",
                                                        "Smolt"))) %>%
                               mutate(life_stage = if_else(llrtp_time < 350 | node %in% (filter(config_file, site_type_name == "AvianColony") %>% pull(node)),
                                                           "juv",
                                                           "adult"))
                           }))

# add capture histories
model_df = obs_df %>%
  mutate(cjs_df = map(llrtp_start,
                      .f = function(x) {
                        ch_df = x %>%
                          filter(emig_stage != "Parr") %>%
                          mutate(across(emig_stage,
                                        fct_drop)) %>%
                          left_join(config_file %>%
                                      arrange(site_code, node, rkm) %>%
                                      group_by(node) %>%
                                      filter(!is.na(rkm),
                                             rkm != "*") %>%
                                      slice(1) %>%
                                      ungroup() %>%
                                      select(node, rkm) %>%
                                      separate(rkm,
                                               into = paste("rkm", 1:2, sep = "_")) %>%
                                      mutate(across(starts_with('rkm'),
                                                    as.numeric)),
                                    by = "node") %>%
                          mutate(node = if_else(rkm_1 == 522 & rkm_2 == 173 & life_stage == "adult",
                                                "GRA",
                                                if_else(rkm_1 == 522 & rkm_2 > 173 & life_stage == "adult",
                                                        "above_GRA",
                                                        node))) %>%
                          # mutate(node = if_else((rkm_1 < 522 | (rkm_1 == 522 & rkm_2 < 173)) & life_stage == "juvenile",
                          #                       "BOJ",
                          #                       if_else((rkm_1 < 522 | (rkm_1 == 522 & rkm_2 < 173)) & life_stage == "adult",
                          #                               "BON",
                          #                               node))) %>%
                          filter(node %in% c("LLRTP", "GRJ", "BOJ", "BON", "GRA", "above_GRA")) %>%
                          mutate(node = factor(node,
                                               levels = c("LLRTP", "GRJ", "BOJ", "BON", "GRA", "above_GRA"))) %>%
                          select(tag_code, emig_stage, node) %>%
                          distinct() %>%
                          mutate(seen = 1) %>%
                          pivot_wider(names_from = node,
                                      values_from = seen,
                                      values_fill = 0,
                                      names_sort = T)

                        # for(node in c("BON", "GRA", "above_GRA")) {
                        #   if(!node %in% names(ch_df)) {
                        #     ch_df %>%
                        #       add_column(glue::glue({node}),
                        #                  .after = Inf)
                        #   }
                        # }

                        ch_df %<>%
                          unite(ch,
                                -c(tag_code:emig_stage),
                                sep = "",
                                remove = F)
                        return(ch_df)
                      })) %>%
  # determine if that brood year can estimate survival to Lower Granite, or Granite-to-Granite SAR
  mutate(incl_GRJ = map_lgl(cjs_df,
                            .f = function(x) {
                              "GRJ" %in% names(x)
                            }),
         incl_GRA = map_lgl(cjs_df,
                            .f = function(x) {
                              "GRA" %in% names(x)
                            }))


model_df %<>%
  mutate(cjs_mod = map(cjs_df,
                       .f = function(cjs_df) {
                         if(n_distinct(cjs_df$emig_stage) == 1) {
                           cat("Only one life-stage\n")
                           cjs_proc = cjs_df %>%
                             select(-tag_code) %>%
                             as.data.frame() %>%
                             process.data()

                           cjs_ddl = make.design.data(cjs_proc)

                           mod = try(crm(data = cjs_proc,
                                         ddl = cjs_ddl,
                                         model.parameters = list(Phi = list(formula = ~ time),
                                                                 p = list(formula = ~ time)),
                                         accumulate = F,
                                         itnmax = 1e4) %>%
                                       cjs.hessian())
                         }

                         cjs_proc = cjs_df %>%
                           select(-tag_code) %>%
                           as.data.frame() %>%
                           process.data(group = c("emig_stage"))

                         cjs_ddl = make.design.data(cjs_proc)

                         mod = try(crm(data = cjs_proc,
                                       ddl = cjs_ddl,
                                       model.parameters = list(Phi = list(formula = ~ emig_stage * time),
                                                               p = list(formula = ~ emig_stage * time)),
                                       accumulate = F,
                                       itnmax = 1e4) %>%
                                     cjs.hessian())

                         if(class(mod)[1] == "try-error") {
                           mod = try(crm(data = cjs_proc,
                                         ddl = cjs_ddl,
                                         model.parameters = list(Phi = list(formula = ~ time),
                                                                 p = list(formula = ~ time)),
                                         accumulate = F,
                                         itnmax = 1e4) %>%
                                       cjs.hessian())
                         }
                         return(mod)
                       })) %>%
  mutate(pred = map(cjs_mod,
                    .f = function(mod) {
                      predict(mod,
                              newdata = data.frame(emig_stage = factor(c("Presmolt", "Smolt"))),
                              se = TRUE) %>%
                        map_df(.id = "param",
                               .f = identity)
                    })) %>%
  mutate(bon_sar = map(pred,
                       .f = function(x) {
                         try(x %>%
                               filter(param == "Phi",
                                      occ <= 3) %>%
                               group_by(emig_stage) %>%
                               summarize(across(c(estimate, lcl, ucl),
                                                prod)))
                       }),
         gra_sar = map(pred,
                       .f = function(x) {
                         try(x %>%
                               filter(param == "Phi",
                                      occ <= 4) %>%
                               group_by(emig_stage) %>%
                               summarize(across(c(estimate, lcl, ucl),
                                                prod)))
                       }))

# detection probability and survival to LGR
dodge_width = 0.4
model_df %>%
  filter(incl_GRJ) %>%
  select(brood_year,
         pred) %>%
  unnest(pred) %>%
  mutate(emig_stage = if_else(is.na(emig_stage) & brood_year == 2010,
                              "Presmolt",
                              if_else(is.na(emig_stage) & brood_year == 2011,
                                      "Smolt",
                                      as.character(emig_stage)))) %>%
  mutate(across(emig_stage,
                as_factor)) %>%
  filter((param == "Phi" & occ == 1) |
           (param == "p" & occ == 2)) %>%
  mutate(param_nm = recode(param,
                           "p" = "Detection at GRJ",
                           "Phi" = "Survival to GRJ")) %>%
  ggplot(aes(x = brood_year,
             y = estimate,
             color = emig_stage,
             fill = emig_stage)) +
  geom_errorbar(aes(ymin = lcl,
                    ymax = ucl),
                width = 0.3,
                position = position_dodge(width = dodge_width)) +
  # geom_crossbar(aes(ymin = estimate - se,
  #                   ymax = estimate + se),
  #               alpha = 0.5) +
  geom_point(size = 4,
             position = position_dodge(width = dodge_width)) +
  scale_color_brewer(palette = "Set1",
                     name = "Life stage") +
  scale_fill_brewer(palette = "Set1",
                    name = "Life stage") +
  facet_wrap(~ param_nm,
             ncol = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Brood Year",
       y = "Probability")

# look at travel times to Lower Granite
travel_time = model_df %>%
  select(brood_year, llrtp_start) %>%
  unnest(llrtp_start) %>%
  group_by(tag_code) %>%
  mutate(grj_det = if_else(sum(node == "GRJ") > 0, T, F)) %>%
  ungroup() %>%
  filter(grj_det) %>%
  group_by(brood_year, tag_code, emig_stage) %>%
  summarize(grj_time = difftime(min_det[node == "GRJ"], max_det[node=="LLRTP"], units = "days"),
            .groups = "drop") %>%
  mutate(across(grj_time,
                as.numeric)) %>%
  group_by(brood_year, emig_stage) %>%
  summarise(across(grj_time,
                   list(mean = mean,
                        median = median)),
            .groups = "drop")

# daily survival rates
model_df %>%
  select(brood_year,
         pred) %>%
  unnest(pred) %>%
  mutate(emig_stage = if_else(is.na(emig_stage) & brood_year == 2010,
                              "Presmolt",
                              if_else(is.na(emig_stage) & brood_year == 2011,
                                      "Smolt",
                                      as.character(emig_stage)))) %>%
  mutate(across(emig_stage,
                as_factor)) %>%
  filter((param == "Phi" & occ == 1)) %>%
  select(-se) %>%
  left_join(travel_time) %>%
  mutate(across(c(estimate, lcl, ucl),
                ~ . ^ (1/grj_time_mean))) %>%
  ggplot(aes(x = brood_year,
             y = estimate,
             color = emig_stage,
             fill = emig_stage)) +
  geom_errorbar(aes(ymin = lcl,
                    ymax = ucl),
                width = 0.3,
                position = position_dodge(width = dodge_width)) +
  geom_point(size = 4,
             position = position_dodge(width = dodge_width)) +
  scale_color_brewer(palette = "Set1",
                     name = "Life stage") +
  scale_fill_brewer(palette = "Set1",
                    name = "Life stage") +
  theme_bw() +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0.9, 1)) +
  labs(x = "Brood Year",
       y = "Probability",
       title = "Daily Survival Rate")

#---------------------------------------------------------
# save some things
save(lem_tags, config_file, obs_df, model_df, travel_time,
     file = here("analysis/data/derived_data/CJS_model_fits.rda"))

#---------------------------------------------------------
# combine all years, and fit a CJS model with a year effect
# cjs_df = model_df %>%
#   select(brood_year, cjs_df) %>%
#   unnest(cjs_df) %>%
#   mutate(across(-c(brood_year:emig_stage),
#                 replace_na,
#                 0)) %>%
#   select(-ch) %>%
#   unite(ch,
#         -c(brood_year:emig_stage),
#         sep = "",
#         remove = F)
#
# cjs_proc = cjs_df %>%
#   select(-tag_code) %>%
#   mutate(across(brood_year,
#                 as.factor)) %>%
#   as.data.frame() %>%
#   process.data(group = c("brood_year", "emig_stage"))
#
# cjs_ddl = make.design.data(cjs_proc)
#
# mod = try(crm(data = cjs_proc,
#               ddl = cjs_ddl,
#               model.parameters = list(Phi = list(formula = ~ (brood_year + emig_stage) * time),
#                                       p = list(formula = ~ (brood_year + emig_stage) * time)),
#               accumulate = F,
#               itnmax = 1e4) %>%
#             cjs.hessian())
#
# mod_pred = tidyr::expand(cjs_df, brood_year, emig_stage) %>%
#   mutate(across(brood_year,
#                 as_factor)) %>%
#   filter(emig_stage != "Parr") %>%
#   mutate(across(emig_stage,
#                 fct_drop)) %>%
#   # filter(brood_year )
#   predict(mod,
#           newdata = .,
#           se = TRUE) %>%
#   map_df(.id = "param",
#          .f = identity)

#---------------------------------------------------------
# look at run timing statistics
obs_df %>%
  select(brood_year,
         llrtp_start) %>%
  unnest(llrtp_start) %>%
  filter(emig_stage != "Parr") %>%
  filter(node %in% c("GRJ", "BOJ", "BON")) %>%
  mutate(node = recode(node,
                       "GRJ" = "Lower Granite",
                       "BOJ" = "Bonneville - Juv.",
                       "BON" = "Bonneville - Adult"),
         node = factor(node,
                       levels = c("Lower Granite",
                                  "Bonneville - Juv.",
                                  "Bonneville - Adult"))) %>%
  mutate(min_jday = yday(min_det)) %>%
  ggplot(aes(x = min_jday,
             y = fct_rev(as.factor(brood_year)),
             fill = emig_stage)) +
  geom_boxplot() +
  facet_wrap(~ node) +
  theme(legend.position = "bottom") +
  labs(x = "Julian Day of Arrival",
       y = "Brood Year",
       fill = "Emigrant Stage",
       title = "DSR/NRR Run Timing")

#---------------------------------------------------------
# look at proportion of fish by life history that come back to BON
# as 3 vs. 4 vs. 5 year olds
ocean_df = obs_df %>%
  select(brood_year,
         llrtp_start) %>%
  unnest(llrtp_start) %>%
  filter(node == "BON",
         life_stage == "adult") %>%
  mutate(return_age = if_else(emig_stage == "Smolt",
                              if_else(llrtp_time < 500,
                                      3,
                                      if_else(between(llrtp_time, 500, 900),
                                              4,
                                              5)),
                              if_else(emig_stage == "Presmolt",
                                      if_else(llrtp_time < 700,
                                              3,
                                              if_else(between(llrtp_time, 700, 1200),
                                                      4, 5)),
                                      if_else(emig_stage == "Parr",
                                              if_else(llrtp_time < 700,
                                                      3,
                                                      if_else(between(llrtp_time, 700, 1200),
                                                              4, 5)),
                                              NA_real_))))

ocean_df %>%
  tabyl(emig_stage, return_age) %>%
  adorn_percentages(denominator = "row") %>%
  # adorn_percentages(denominator = "col") %>%
  adorn_pct_formatting()

ocean_df %>%
  ggplot(aes(x = llrtp_time,
             fill = emig_stage,
             color = as.factor(return_age))) +
  scale_color_brewer(palette = "Set1",
                     name = "Return Age") +
  scale_fill_brewer(palette = "Set2",
                    name = "Emig. Stage") +
  geom_histogram() +
  labs(x = "Days Since LLRTP")

ocean_df %>%
  group_by(emig_stage, return_age) %>%
  summarize(n_tags = n(),
            .groups = "drop") %>%
  group_by(emig_stage) %>%
  mutate(prop = n_tags / sum(n_tags)) %>%
  ggplot(aes(x = emig_stage,
             y = prop,
             fill = as.factor(return_age))) +
  geom_col(position = "stack") +
  scale_fill_brewer(palette = "Set1",
                    name = "Return Age") +
  labs(x = "Emigration Stage",
       y = "Proportion")

ocean_df %>%
  filter(emig_stage != "Parr",
         brood_year < 2018) %>%
  group_by(emig_stage, brood_year, return_age) %>%
  summarize(n_tags = n(),
            .groups = "drop") %>%
  group_by(emig_stage, brood_year) %>%
  mutate(prop = n_tags / sum(n_tags)) %>%
  ggplot(aes(x = emig_stage,
             y = prop,
             fill = as.factor(return_age))) +
  geom_col(position = "stack") +
  facet_wrap(~ brood_year) +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Return Age") +
  labs(x = "Emigration Stage",
       y = "Proportion")

ocean_df %>%
  group_by(brood_year, emig_stage, return_age) %>%
  summarize(n_tags = n(),
            .groups = "drop") %>%
  group_by(brood_year, return_age) %>%
  mutate(tot_tags = sum(n_tags),
         prop = n_tags / tot_tags) %>%
  group_by(brood_year) %>%
  mutate(n_lh = n_distinct(emig_stage)) %>%
  ungroup() %>%
  filter(n_lh > 1) %>%
  # filter(n_lh == 1)
  select(-n_tags) %>%
  pivot_wider(names_from = emig_stage,
              values_from = prop,
              values_fill = 0) %>%
  filter(return_age == 3)

model_df %>%
  select(brood_year,
         cjs_df) %>%
  unnest(cjs_df) %>%
  filter(emig_stage != "Parr") %>%
  mutate(across(c(BOJ:above_GRA),
                replace_na,
                0)) %>%
  pivot_longer(cols = LLRTP:above_GRA,
               names_to = "node",
               values_to = "seen") %>%
  filter(seen == 1) %>%
  mutate(across(node,
                factor,
                levels = c("LLRTP",
                           "GRJ",
                           "BOJ",
                           "BON",
                           "GRA",
                           "above_GRA"))) %>%
  group_by(brood_year, emig_stage, node) %>%
  summarise(n_obs = n_distinct(tag_code),
            .groups = "drop") %>%
  pivot_wider(names_from = "node",
              values_from = n_obs,
              values_fill = 0,
              names_sort = T)

