na_remove <- function(data) {
  mutate_if(data, is.numeric, ~ ifelse(. %in% c(-99, -88, -77), NA, .))
}

na_delete <- function(data) {
  r_na <- apply(data, 1, function(x) {
    all(is.na(x))
  }) %>% which()
  c_na <- apply(data, 1, function(x) {
    all(is.na(x))
  }) %>% which()
  return(data[-r_na, -c_na])
}
