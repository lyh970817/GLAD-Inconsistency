questionnaires <- ls() %>% grep("[A-Z]+", ., v = T)
raw_dat_list <- mget(questionnaires)
cache("raw_dat_list")

# sheet_list <- GLAD_sheet(questionnaires)
cache("sheet_list")

item_list <- map2(raw_dat_list, sheet_list, function(d, s) {
  if ("score_key" %in% colnames(s)) {
    items <- s %>%
      filter(score_key %in% c(-1, 1)) %>%
      select(easyname) %>%
      pull() %>%
      unique() %>%
      subset(!is.na(.) & . %in% colnames(d))
    return(items)
  }
})
cache("item_list")

descr_list <- map2(item_list, sheet_list, ~ GLAD_getdescr(.x, .y))
cache("sheet_list")

dat_list <- map2(item_list, raw_dat_list, function(i, d) {
  num_items <- paste(i, "numeric", sep = "_")
  if (length(i) > 0) {
    d[c("ID", num_items)] %>%
      na_remove() %>%
      return()
  }
})

dat_list_lgl <- map_lgl(dat_list, ~ !is.null(.x))
dat_list <- dat_list %>% subset(dat_list_lgl)
cache("dat_list")

dat <- reduce(dat_list, full_join)
cache("dat")
