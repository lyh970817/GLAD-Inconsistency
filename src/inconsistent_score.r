setwd("../")
library(ProjectTemplate)
library(polycor)
library(igraph)
load.project()

dat_fct <- mutate_if(dat, is.numeric, factor)

# cor_mat <- hetcor(dat_fct[-1], use = "pairwise.complete.obs", std.err = FALSE)
cache("cor_mat")

cor_m <- cor_mat$correlations
cor_m[abs(cor_m) < 0.6] <- 0
cor_m[is.na(cor_m)] <- 0
diag(cor_m) <- 0

no_connect_row <- apply(cor_m, 1, function(r) {
  all(r == 0)
})

no_connect_col <- apply(cor_m, 2, function(r) {
  all(r == 0)
})

rm <- which(no_connect_row & no_connect_col)

cor_m <- cor_m[-rm, -rm]

item_names <- dimnames(cor_m)[[1]]

dimnames(cor_m) <- NULL

g_cor <- graph_from_adjacency_matrix(cor_m, mode = "undirected", diag = F, weighted = T)

match_cor <- blossom(g_cor, weighted = T)
pairs <- item_names[match_cor$matching]

first <- pairs[seq(1, length(pairs), by = 2)]
second <- pairs[seq(2, length(pairs), by = 2)]

n_levels <- map_dbl(pairs, ~ range(dat[[.x]], na.rm = T) %>% diff())
which(n_levels[seq(1, length(pairs), by = 2)] != n_levels[seq(2, length(pairs), by = 2)])

first[8]
n_levels[seq(1, length(pairs), by = 2)][8]
second[8]
n_levels[seq(2, length(pairs), by = 2)][8]

dat[[first[8]]] <- ifelse(dat[[first[8]]] == 1, 3, dat[[first[8]]])

first[11]
n_levels[seq(1, length(pairs), by = 2)][11]
second[11]
n_levels[seq(2, length(pairs), by = 2)][11]

dat[[first[11]]] <- ifelse(dat[[first[11]]] == 1, 3, dat[[first[11]]])

first[38]
n_levels[seq(1, length(pairs), by = 2)][38]
second[38]
n_levels[seq(2, length(pairs), by = 2)][38]

dat[[first[38]]] <- ifelse(dat[[first[38]]] %in% c(2, 3), 2, dat[[first[38]]])
dat[[first[38]]] <- ifelse(dat[[first[38]]] == 4, 3, dat[[first[38]]])
dat[[first[38]]] <- ifelse(dat[[first[38]]] == 5, 4, dat[[first[38]]])

# If we want to compare say between gad scoring highest and all, gad pairs should be
# removed because they will definitely be 0 (scoring highest means answered
# the same to all quetsions)
sub_lgl <- (grepl("phq", first) & grepl("phq", second)) | (grepl("gad", first) & grepl("gad", second))

first_sub <- first[!sub_lgl]
second_sub <- second[!sub_lgl]

diff_score <- map2(first_sub, second_sub, function(x, y) {
  diff <- (abs(dat[[x]] - dat[[y]]))
  norm_diff <- (diff - min(diff, na.rm = T)) / max(diff, na.rm = T)
  return(norm_diff)
}) %>%
  bind_cols() %>%
  rowSums(na.rm = T)

sum(table(diff_score) %>% .[as.numeric(dimnames(.)[[1]]) > 13])

# hist(diff_score)

diff_dat <- tibble(ID = dat[["ID"]], diff_score = diff_score)

dat_gad <- left_join(raw_dat_list[["GAD"]], diff_dat)
dat_gad <- GLAD_derive(dat_gad, sheet_list[["GAD"]])
max_score_gad <- max(dat_gad[["gad.total_score"]], na.rm = T)

max_gad <- dat_gad %>%
  filter(gad.total_score %in% max_score_gad) %>%
  select(ID) %>%
  pull()

t.test(
  filter(diff_dat, ID %in% max_gad)[["diff_score"]],
  filter(diff_dat, !ID %in% max_gad)[["diff_score"]]
)

dat_phq <- left_join(raw_dat_list[["PHQ"]], diff_dat)
dat_phq <- GLAD_derive(dat_phq, sheet_list[["PHQ"]])
max_score_phq <- max(dat_phq[["phq.new.total"]], na.rm = T)

max_phq <- dat_phq %>%
  filter(phq.new.total %in% max_score_phq) %>%
  select(ID) %>%
  pull()

t.test(
  filter(diff_dat, ID %in% max_phq)[["diff_score"]],
  filter(diff_dat, !ID %in% max_phq)[["diff_score"]]
)

alpha_subset <- function(dat_list, ids = NULL) {
  cor_list <- map(dat_list, function(d) {
    if (!is.null(ids)) {
      d <- d %>%
        filter(ID %in% ids)
    }

    cor <- hetcor(modify(d[-1], factor),
      use = "pairwise.complete.obs",
      std.err = FALSE
    )$correlations
    cor[is.na(cor)] <- 0
    return(cor)
  })

  alphas <- map_dbl(cor_list, ~ psych::alpha(.x, check.keys = T)$total$raw_alpha)
  return(alphas)
}

gad_alphas <- alpha_subset(dat_list, max_gad)
phq_alphas <- alpha_subset(dat_list, max_phq)
all_alphas <- alpha_subset(dat_list)

mean(gad_alphas[gad_alphas >= 0.5])
mean(phq_alphas[phq_alphas >= 0.5])
mean(all_alphas)
