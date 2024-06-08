
key <- data.frame(
  filename = c("model_April2024", "modelannual_April2024", "modelE_April2024", "modelW_April2024", "modelWmult_April2024"),
  name = c("MS-seasonal", "MS-annual", "EATL-seasonal", "WATL-seasonal", "WATL-seasonal-multi"),
  EBFT = c(1, 1, 1, NA, NA),
  WBFT = c(2, 2, NA, 1, 1)
)

key <- data.frame(
  filename = c("model_April2024", "modelannual_April2024", "modelE_April2024", "modelW_April2024"),
  name = c("MS-seasonal", "MS-annual", "EATL-seasonal", "WATL-seasonal"),
  EBFT = c(1, 1, 1, NA),
  WBFT = c(2, 2, NA, 1)
)

fits <- lapply(key$filename, function(i) readRDS(paste0("fit/", i, ".rds")))

S_ys <- lapply(1:nrow(key), function(i) {

  year <- get_MARSdata(fits[[i]])@Dlabel@year
  if (!is.na(key$EBFT[i])) {

    s <- key$EBFT[i]

    EBFT <- data.frame(
      year = year,
      est = fits[[i]]@report$S_yrs[, , s, drop = FALSE] %>% apply(1, sum),
      ref = fits[[i]]@report$phi_s[s] * fits[[i]]@report$R0_s[s],
      stock = "EBFT"
    )

  } else {
    EBFT <- data.frame()
  }

  if (!is.na(key$WBFT[i])) {

    s <- key$WBFT[i]

    WBFT <- data.frame(
      year = year,
      est = fits[[i]]@report$S_yrs[, , s, drop = FALSE] %>% apply(1, sum),
      ref = fits[[i]]@report$phi_s[s] * fits[[i]]@report$R0_s[s],
      stock = "WBFT"
    )

  } else {
    WBFT <- data.frame()
  }

  rbind(EBFT, WBFT) %>%
    mutate(model = key$name[i], ratio = est/ref)
}) %>%
  bind_rows() %>%
  mutate(model = factor(model, levels = key$name))

p <- mutate(S_ys, prop = est/sum(est), .by = c(year, model))

g1 <- ggplot(S_ys, aes(year, est, colour = model)) +
  facet_wrap(vars(stock), scales = "free_y") +
  geom_line() +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Spawning biomass", colour = "Model") +
  theme(legend.position = "bottom")
ggsave("fit/compare_SB.png", g1, height = 3.5, width = 6)

g2 <- ggplot(S_ys, aes(year, ratio, colour = model)) +
  facet_wrap(vars(stock)) +
  geom_line() +
  expand_limits(y = 0) +
  labs(x = "Year", y = expression(B/B[0]), colour = "Model") +
  theme(legend.position = "bottom")
ggsave("fit/compare_dep.png", g2, height = 3.5, width = 6)
