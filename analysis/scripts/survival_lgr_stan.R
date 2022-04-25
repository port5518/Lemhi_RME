# Author: Kevin See
# Purpose: Survival from LLRTP to LGR - Bayesian style with Stan
# Created: 6/29/21
# Last Modified: 4/25/22
# Notes:

# clear environment
rm(list = ls())

#-----------------------------------------------------------------
# load packages
library(tidyverse)
library(here)
library(magrittr)
library(janitor)
library(lubridate)
# remotes::install_github("mackerman44/PITcleanr@develop")
library(PITcleanr)
library(rstan)
library(shinystan)

# set some rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#-----------------------------------------------------------------
# load prepped data
load(here("analysis/data/derived_data/CJS_model_fits.rda"))

# create some inputs
# pull together capture histories
ch_df = model_df %>%
  select(brood_year, cjs_df) %>%
  unnest(cjs_df) %>%
  filter(emig_stage != "Parr") %>%
  mutate(across(emig_stage,
                fct_drop)) %>%
  mutate(across(c(BON, GRA, above_GRA),
                replace_na,
                0)) %>%
  select(-ch) %>%
  unite(ch,
        -c(brood_year:emig_stage),
        sep = "",
        remove = F)

#-----------------------------------------------------------------
# fit model to all years separately
# set parameters to track
param_names <- c("phi",
                 "p",
                 "beta",
                 "delta_llrtp_grj",
                 "delta_grj_boj",
                 "delta_llrtp_boj",
                 "delta_llrtp_bon",
                 "delta_llrtp_gra",
                 "delta_sar_bon",
                 "delta_sar_gra")

# what is the name of the stan model file
# model_file = here("analysis/data/jags/cjs.stan")
model_file = here("analysis/data/jags/cjs_lh.stan")

# fit models and extract some summaries
stan_df = ch_df %>%
  # filter(!brood_year %in% c(2010, 2011)) %>%
  nest(data = -brood_year) %>%
  mutate(stan_data = map(data,
                         .f = function(ch_df) {
                           y = ch_df %>%
                             select(-c(tag_code:ch)) %>%
                             as.matrix()

                           # occasions
                           T = ncol(y)
                           # num. individuals
                           I = nrow(y)

                           # life-history type (DSR or NRR)
                           lh = ch_df %>%
                             mutate(across(emig_stage,
                                           fct_drop)) %>%
                             pull(emig_stage) %>%
                             as.numeric()

                           n_lh = n_distinct(lh)

                           list(y = y,
                                T = T,
                                I = I,
                                lh = lh,
                                n_lh = n_lh) %>%
                             return()
                         })) %>%
  mutate(stan_fit = map(stan_data,
                        .f = function(x) {
                          fit = stan(file = model_file,
                                     data = x,
                                     pars = param_names,
                                     chains = 3,
                                     cores = 3)
                          return(fit)
                        }))

# summarize stan results
stan_df %<>%
  mutate(stan_summ = map(stan_fit,
                         .f = function(fit) {
                           fit_summ = try(summary(fit)$summary %>%
                                            as_tibble(rownames = "param") %>%
                                            add_column(emig_stage = NA_character_,
                                                       .after = "param") %>%
                                            mutate(emig_stage = if_else(grepl(",1\\]$", param),
                                                                        "Presmolt",
                                                                        if_else(grepl(",2\\]$", param),
                                                                                "Smolt",
                                                                                NA_character_))))
                           if(class(fit_summ)[1] == "try-error") {
                             return(NULL)
                           } else {
                             return(fit_summ)
                           }
                         }))

# for a couple years with only one life history, fix that in the summary
stan_df %<>%
  mutate(stan_summ = map2(stan_summ, brood_year,
                          .f = function(fit_summ, yr) {
                            if(yr == 2010) {
                              fit_summ %<>%
                                mutate(emig_stage = "Presmolt") %>%
                                return()
                            } else if(yr == 2011) {
                              fit_summ %<>%
                                mutate(emig_stage = "Smolt") %>%
                                return()
                            } else {
                              return(fit_summ)
                            }
                          }))


# save some things
save(ch_df, stan_df,
     file = here("analysis/data/derived_data/CJS_Stan_fits.rda"))

#-----------------------------------------------------------------
# diagnostics with shinystan
stan_df %>%
  filter(brood_year == 2009) %>%
  pull(stan_fit) %>%
  extract2(1) %>%
  shinystan::launch_shinystan()

#-----------------------------------------------------------------
# DSR/NRR comparison
comp_df = stan_df %>%
  filter(!brood_year %in% c(2010, 2011)) %>%
  select(brood_year,
         stan_summ) %>%
  filter(brood_year < 2018) %>%
  unnest(stan_summ) %>%
  filter(grepl('delta', param)) %>%
  mutate(param = fct_recode(param,
                            "LLRTP to GRJ" = "delta_llrtp_grj",
                            "GRJ to BOJ" = "delta_grj_boj",
                            "LLRTP to BOJ" = "delta_llrtp_boj",
                            "LLRTP to BON" = "delta_llrtp_bon",
                            "LLRTP to GRA" = "delta_llrtp_gra",
                            "SAR (GRJ-BON)" = "delta_sar_bon",
                            "SAR (GRJ-GRA)" = "delta_sar_gra")) %>%
  mutate(across(param,
                fct_relevel,
                c("LLRTP to GRJ"),
                after = 0)) %>%
  mutate(across(brood_year,
                as_factor)) %>%
  mutate(winner = if_else(mean > 0 & `50%` > 0,
                          "DSR",
                          if_else(mean < 0 & `50%` < 0,
                                  "NRR",
                                  NA_character_)))
delta_surv_p = comp_df %>%
  filter(!grepl("SAR", param)) %>%
  ggplot(aes(x = brood_year)) +
  geom_boxplot(aes(ymin = `2.5%`,
                   lower = `25%`,
                   middle = mean,
                   upper = `75%`,
                   ymax = `97.5%`,
                   fill = winner),
               stat = "identity") +
  facet_wrap(~ param,
             nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Better Choice") +
  geom_hline(yintercept = 0,
             linetype = 2) +
  labs(x = "Brood Year",
       y = "<--  Better for DSR ... Better for NRR  -->\n\n\nLog( DSR / NRR )")

ggsave(here("analysis/figures/log_surv_comp.pdf"),
       delta_surv_p,
       width = 12,
       height = 9)

delta_sar_p = comp_df %>%
  filter(grepl("SAR", param)) %>%
  ggplot(aes(x = brood_year)) +
  geom_boxplot(aes(ymin = `2.5%`,
                   lower = `25%`,
                   middle = mean,
                   upper = `75%`,
                   ymax = `97.5%`,
                   fill = winner),
               stat = "identity") +
  facet_wrap(~ param,
             ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Better Choice") +
  geom_hline(yintercept = 0,
             linetype = 2) +
  labs(x = "Brood Year",
       y = "<--  Better for DSR ... Better for NRR  -->\n\n\nLog( DSR / NRR )")

ggsave(here("analysis/figures/log_sar_comp.pdf"),
       delta_sar_p,
       width = 9,
       height = 9)

delta_fresh_p = comp_df %>%
  filter(param %in% c("LLRTP to GRJ",
                      "GRJ to BOJ",
                      "LLRTP to BOJ")) %>%
  ggplot(aes(x = brood_year)) +
  geom_boxplot(aes(ymin = `2.5%`,
                   lower = `25%`,
                   middle = mean,
                   upper = `75%`,
                   ymax = `97.5%`,
                   fill = winner),
               stat = "identity") +
  facet_wrap(~ param,
             ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Better Choice") +
  geom_hline(yintercept = 0,
             linetype = 2) +
  labs(x = "Brood Year",
       y = "<--  Better for DSR ... Better for NRR  -->\n\n\nLog( DSR / NRR )")

ggsave(here("analysis/figures/log_fresh_comp.pdf"),
       delta_fresh_p,
       width = 12,
       height = 9)

stan_df %>%
  select(brood_year,
         stan_summ) %>%
  # filter(brood_year < 2018) %>%
  unnest(stan_summ) %>%
  filter(grepl('phi\\[1,', param)) %>%
  mutate(across(brood_year,
                as_factor)) %>%
  ggplot(aes(x = brood_year)) +
  geom_boxplot(aes(ymin = `2.5%`,
                   lower = `25%`,
                   middle = mean,
                   upper = `75%`,
                   ymax = `97.5%`,
                   fill = emig_stage),
               stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Life History") +
  labs(x = "Brood Year",
       y = "Survival: LLRTP to GRJ")

stan_df %>%
  filter(brood_year < 2018) %>%
  mutate(stan_mcmc = map(stan_fit,
                         .f = function(x) {
                           as.data.frame(x) %>%
                             as_tibble(rownames = "iter")
                         })) %>%
  select(brood_year, stan_mcmc) %>%
  unnest(cols = stan_mcmc) %>%
  select(brood_year, iter,
         starts_with("p[")) %>%
  pivot_longer(cols = starts_with("p["),
               names_to = "param",
               values_to = "value") %>%
  mutate(emig_stage = if_else(str_detect(param, ",1]") & brood_year != 2011,
                              "Presmolt",
                              if_else(str_detect(param, ",2]") | brood_year == 2011,
                                      "Smolt",
                                      NA_character_))) %>%
  mutate(det_site = str_extract(param, "[:digit:]"),
         det_site = recode(det_site,
                           "2" = "GRJ",
                           "3" = "BOJ",
                           "4" = "BON",
                           "5" = "GRA")) %>%
  filter(!det_site %in% c("1", "6")) %>%
  mutate(across(c(brood_year, det_site),
                as_factor)) %>%
  ggplot(aes(x = brood_year,
             y = value,
             fill = emig_stage)) +
  geom_violin(draw_quantiles = c(0.5),
              scale = "area") +
  # geom_boxplot() +
  facet_wrap(~ det_site) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Life History") +
  labs(x = "Brood Year",
       y = "Detection Probability")

stan_df %>%
  filter(brood_year < 2018) %>%
  mutate(stan_mcmc = map(stan_fit,
                         .f = function(x) {
                           as.data.frame(x) %>%
                             as_tibble(rownames = "iter")
                         })) %>%
  select(brood_year, stan_mcmc) %>%
  unnest(cols = stan_mcmc) %>%
  select(brood_year, iter,
         starts_with("phi[")) %>%
  pivot_longer(cols = starts_with("phi["),
               names_to = "param",
               values_to = "value") %>%
  mutate(emig_stage = if_else(str_detect(param, ",1]") & brood_year != 2011,
                              "Presmolt",
                              if_else(str_detect(param, ",2]") | brood_year == 2011,
                                      "Smolt",
                                      NA_character_))) %>%
  mutate(surv_period = str_extract(param, "[:digit:]"),
         surv_period = recode(surv_period,
                              "1" = "LLRTP to GRJ",
                              "2" = "GRJ to BOJ",
                              "3" = "BOJ to BON",
                              "4" = "BON to GRA")) %>%
  filter(!surv_period %in% c("5")) %>%
  mutate(across(c(brood_year, surv_period),
                as_factor)) %>%
  ggplot(aes(x = brood_year,
             y = value,
             fill = emig_stage)) +
  # geom_violin(draw_quantiles = c(0.5),
  #             scale = "area") +
  geom_boxplot() +
  facet_wrap(~ surv_period) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1",
                    name = "Life History") +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Brood Year",
       y = "Survival Probability")

# END SCRIPT
