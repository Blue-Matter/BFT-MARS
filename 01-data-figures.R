
library(tidyverse)

# Get data list from ABFTMSE package ----
dat <- local({
  load("data/MARS_input_March24.rda")
  MARSinput
})

## Figures ----
year_df <- data.frame(Year = 1:dat$ny, Real_year = seq(dat$years[1], dat$years[2]))
area_df <- data.frame(Area = 1:dat$nr, AreaName = dat$areas) %>%
  mutate(AreaName = factor(AreaName, levels = AreaName))
fleet_df <- data.frame(Fleet = 1:dat$nf) %>%
  mutate(FleetName = paste0(Fleet, "-", dat$Fleets$name),
         FleetName = factor(FleetName, levels = FleetName))


# Annual catch
Cobs <- dat$Cobs %>%
  as.data.frame() %>%
  left_join(year_df) %>%
  left_join(area_df) %>%
  left_join(fleet_df) %>%
  mutate(Catch = Catch/1e6)

g <- Cobs %>%
  summarize(Catch = sum(Catch), .by = c(Real_year, AreaName, FleetName)) %>%
  ggplot(aes(Real_year, Catch, fill = FleetName)) +
  geom_col(width = 1, colour = "grey40", linewidth = 0.25) +
  facet_wrap(vars(AreaName)) +
  theme(panel.spacing = unit(0, "in"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, title = "Fleet")) +
  labs(x = "Year", y = "Catch")
ggsave("figures/Cobs_annual.png", g, height = 6, width = 6)

g <- Cobs %>%
  summarize(Catch = sum(Catch), .by = c(Real_year, AreaName, FleetName)) %>%
  ggplot(aes(Real_year, Catch, fill = FleetName)) +
  geom_col(width = 1, colour = "grey40", linewidth = 0.25) +
  facet_wrap(vars(AreaName), scales = "free_y") +
  theme(panel.spacing = unit(0, "in"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, title = "Fleet")) +
  labs(x = "Year", y = "Catch")
ggsave("figures/Cobs_annual2.png", g, height = 6, width = 6)

# Seasonal catch
g <- Cobs %>%
  mutate(Quarter = paste("Season", Quarter)) %>%
  ggplot(aes(Real_year, Catch, fill = FleetName)) +
  geom_col(width = 1, colour = "grey40", linewidth = 0.1) +
  facet_grid(vars(AreaName), vars(Quarter)) +
  theme(panel.spacing = unit(0, "in"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, title = NULL)) +
  labs(x = "Year", y = "Catch")
ggsave("figures/Cobs_season.png", g, height = 6, width = 6)

g <- Cobs %>%
  mutate(Quarter = paste("Season", Quarter)) %>%
  ggplot(aes(Real_year, Catch, fill = FleetName)) +
  geom_col(width = 1, colour = "grey40", linewidth = 0.1) +
  facet_grid(vars(AreaName), vars(Quarter), scales = "free_y") +
  theme(panel.spacing = unit(0, "in"),
        legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 4, title = NULL)) +
  labs(x = "Year", y = "Catch")
ggsave("figures/Cobs_season2.png", g, height = 6, width = 6)

# Spool up catch
g <- dat$HCobs %>%
  structure(dimnames = list(Year = 1864:1964, Season = 1:4, Age = 1:35, Area = area_df$AreaName)) %>%
  reshape2::melt() %>%
  filter(Age %in% seq(1, 35, 5)) %>%
  mutate(Season = paste("Season", Season), Age = factor(Age)) %>%
  ggplot(aes(Year, value, colour = Age)) +
  facet_grid(vars(Area), vars(Season), scales = "free_y") +
  geom_line() +
  theme(panel.spacing = unit(0, "in"), legend.position = "bottom") +
  labs(y = "Catch")
ggsave("figures/catch_reconstruction.png", height = 5, width = 6)

# Indices
#Ft = Fleet selectivity to mirror
#wt = additional weighting to be able to turn things off or emphasize certain indices or parts of series
#lt = length type, the US rod and reel indices are set up in length brackets
#     for example 66cm - 114cm
#     so the selectivity is mirrored to the fleet, but its only a range of lengths in that fleets CPUE observations

Iname_df <- data.frame(Ino = 1:dat$nI, Iname = dat$Inames)
g <- dat$Iobs %>%
  as.data.frame() %>%
  left_join(year_df) %>%
  left_join(Iname_df) %>%
  left_join(area_df, by = c("area" = "Area")) %>%
  mutate(Real_year = Real_year + 0.25*(subyear - 1),
         Stock = ifelse(stock == 1, "EBFT", "WBFT")) %>%
  ggplot(aes(Real_year, index, colour = AreaName, shape = Stock)) +
  geom_line(linewidth = 0.1, linetype = 3, color = "grey40") +
  geom_linerange(aes(ymin = exp(log(index - 2 * CV)), ymax = exp(log(index + 2*CV)))) +
  geom_point() +
  facet_wrap(vars(Iname), scales = "free_y") +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 21)) +
  labs(x = "Year", y = "Stock-specific index", colour = NULL, shape = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/Index.png", g, height = 6, width = 6)


g <- dat$CPUEobs %>%
  as.data.frame() %>%
  left_join(year_df) %>%
  left_join(fleet_df) %>%
  left_join(area_df, by = c("Area" = "Area")) %>%
  mutate(Real_year = Real_year + 0.25*(Quarter - 1),
         CV = pmin(CV, 1)) %>%
  ggplot(aes(Real_year, Index, colour = AreaName, group = AreaName)) +
  geom_line(linewidth = 0.1, linetype = 3, color = "grey40") +
  geom_linerange(aes(ymin = exp(log(Index) - 2 * CV), ymax = exp(log(Index) + 2*CV))) +
  geom_point() +
  facet_wrap(vars(FleetName), ncol = 2, scales = "free_y") +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Fishery CPUE", colour = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/Index2.png", g, height = 9, width = 6)

# CAL ----
len_df <- data.frame(Length_category = 0:dat$nl+1, Bin = dat$lenbins)
g <- dat$CLobs %>%
  as.data.frame() %>%
  left_join(year_df) %>%
  left_join(area_df) %>%
  left_join(len_df) %>%
  left_join(fleet_df) %>%
  mutate(p = N/sum(N), .by = c(Real_year, Subyear, Area, Fleet)) %>%
  mutate(Year = Real_year + 0.25 * (Subyear - 1)) %>%
  ggplot(aes(Bin, N, group = Year, colour = Year)) +
  geom_line(alpha = 0.5) +
  facet_grid(vars(FleetName), vars(Area)) +
  theme(panel.spacing = unit(0, "in"))
ggsave("figures/CAL.png", g, height = 9, width = 6)

# SC ----
# Type 1 = microchemistry (do not use!)
# Type 2 = genetics
g <- dat$SOOobs %>%
  as.data.frame() %>%
  left_join(year_df, by = c("y" = "Year")) %>%
  left_join(area_df, by = c("r" = "Area")) %>%
  filter(Type == 1) %>%
  mutate(Year = Real_year + 0.25 * (s-1), a = paste("Age class", a)) %>%
  ggplot(aes(Year, plogis(probE))) +
  geom_col(width = 0.25, colour = "black", fill = "grey40") +
  facet_grid(vars(AreaName), vars(a)) +
  labs(x = "Year", y = "Probability EBFT") +
  theme(panel.spacing = unit(0, "in")) +
  ggtitle("Stock composition - microchemistry")
ggsave("figures/SC1.png", g, height = 4, width = 6)

g <- dat$SOOobs %>%
  as.data.frame() %>%
  left_join(year_df, by = c("y" = "Year")) %>%
  left_join(area_df, by = c("r" = "Area")) %>%
  filter(Type == 2, !is.na(r)) %>%
  mutate(Year = Real_year + 0.25 * (s-1), a = paste("Age class", a)) %>%
  ggplot(aes(Year, plogis(probE))) +
  geom_col(width = 0.25, colour = "black", fill = "grey40") +
  #geom_point() + geom_line() +
  facet_grid(vars(AreaName), vars(a)) +
  labs(x = "Year", y = "Probability EBFT") +
  theme(panel.spacing = unit(0, "in")) +
  ggtitle("Stock composition - genetics")
ggsave("figures/SC2.png", g, height = 3, width = 6)


# etag ----
gE <- dat$PSAT %>%
  as.data.frame() %>%
  #left_join(year_df, by = c("y" = "Year")) %>%
  left_join(area_df, by = c("fr" = "Area")) %>%
  left_join(area_df, by = c("tr" = "Area")) %>%
  rename(AreaFrom = AreaName.x, AreaTo = AreaName.y) %>%
  filter(p == 1) %>%
  mutate(s = paste("Season", s), a = paste("Age class", a)) %>%
  ggplot(aes(AreaTo, AreaFrom)) +
  geom_tile(aes(fill = FR), colour = "black") +
  geom_text(aes(label = round(FR, 2))) +
  facet_grid(vars(s), vars(a)) +
  ggtitle("EBFT") +
  coord_cartesian(expand = FALSE) +
  theme(panel.background = element_rect(fill = "grey90"),
        legend.position = "bottom",
        strip.background = element_blank()) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5) +
  labs(x = "Destination", y = "Origin", fill = "Proportion")
ggsave("figures/psat_east.png", gE, height = 5, width = 5)

gW <- dat$PSAT %>%
  as.data.frame() %>%
  #left_join(year_df, by = c("y" = "Year")) %>%
  left_join(area_df, by = c("fr" = "Area")) %>%
  left_join(area_df, by = c("tr" = "Area")) %>%
  rename(AreaFrom = AreaName.x, AreaTo = AreaName.y) %>%
  filter(p == 2) %>%
  mutate(s = paste("Season", s), a = paste("Age class", a)) %>%
  ggplot(aes(AreaTo, AreaFrom)) +
  geom_tile(aes(fill = FR), colour = "black") +
  geom_text(aes(label = round(FR, 2))) +
  facet_grid(vars(s), vars(a)) +
  ggtitle("EBFT") +
  coord_cartesian(expand = FALSE) +
  theme(panel.background = element_rect(fill = "grey90"),
        legend.position = "bottom",
        strip.background = element_blank()) +
  scale_fill_gradient2(low = "blue", high = "red", midpoint = 0.5) +
  labs(x = "Destination", y = "Origin", fill = "Proportion")
ggsave("figures/psat_west.png", gW, height = 5, width = 4)


