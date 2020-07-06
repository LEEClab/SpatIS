######################################################
#
# Linking individual specialization to space use
# for Sturnira lilium bats
# 
# Testing the relationship between food plant density
# and the use of foraging areas by individuals
#
# Bernardo Niebuhr - bernardo_brandaum at yahoo.com.br
# Patricia Rogeri - parogeri at gmail.com
#
# Jun 2020
# No copyrights - feel free to use, modify, and share
#######################################################

# Load pacakges
if(!require(tidyverse)) install.packages("tidyverse", dep=T); library(tidyverse)
if(!require(broom)) install.packages("broom", dep=T); library(broom)

# Load data
use.areas <- read.csv("data/proportion_use_foraging_sites.csv", header = T)
fruit.density <- read.csv("data/use_foraging_areas_density_fruits.csv", header = T)

# get values for all areas for individual (adding zeros)
inds <- unique(use.areas$Animal_ID)
areas <- sort(unique(use.areas$area))
inds.areas <- paste(use.areas$Animal_ID, use.areas$area)
for(i in inds){
  line <- use.areas %>% 
    dplyr::filter(Animal_ID == i) %>% 
    slice(1)
  for(j in areas) {
    if(!(paste(i, j) %in% inds.areas)) {
      line$area <- j
      line$success <- 0
      line$failure <- line$sample.size
      use.areas <- dplyr::bind_rows(use.areas, line)
    }
  }
}
use.areas <- use.areas %>% 
  dplyr::arrange(Animal_ID, areas)

# merge
use.density <- use.areas %>% 
  dplyr::left_join(fruit.density, 
                   by = c("area" = "Foraging_area")) %>% 
  dplyr::as_tibble() %>% 
  tidyr::pivot_longer(cols = c(Solanum, Piper, Cecropia),
                      names_to = "plant",
                      values_to = "density")

# plot
use.density %>% 
  ggplot(aes(x = density, y = success)) +
  facet_wrap(vars(Animal_ID, plant)) +
  geom_point()

g1 <- use.density %>% 
  dplyr::filter(!(density == 80 & !(Animal_ID %in% c(4, 5)))) %>% 
  ggplot(aes(x = density, y = success, color = plant)) +
  facet_wrap(~Animal_ID, scales = "free") +
  geom_point() + 
  theme_bw() +
  labs(x = "Fruit density",
       y = "Number of activity points", 
       color = "Food-plant species") + 
  geom_smooth(method = "lm", se = TRUE) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed")
g1

ggsave("use_density_plots.png", plot = g1, path = "output", 
       width = 20, height = 13, units = "cm", dpi = 600)

# fit

# one model per individual and plant species
fits <- use.density %>% 
  dplyr::group_by(Animal_ID, plant) %>% 
  tidyr::nest_legacy() %>% 
  dplyr::mutate(fit = purrr::map(data, ~ lm(.x$success ~ .x$density)),
                tidied = purrr::map(fit, broom::tidy),
                glanced = purrr::map(fit, broom::glance)) %>% 
  dplyr::select(-data, -fit) %>% 
  tidyr::unnest_legacy(glanced) %>% 
  dplyr::select(Animal_ID, plant, tidied, r.squared, df) %>% 
  tidyr::unnest_legacy(tidied)

(to.save <- fits %>% 
  dplyr::filter(!grepl("Int", term)) %>% 
  dplyr::select(-term) %>% 
  dplyr::select(-std.error, -statistic))

to.save %>% 
  dplyr::filter(p.value <= 0.05)

