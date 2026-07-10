# =====================================================================
# 00_export_from_db.R  —  从 cesd_analysis.db 导出 Phase-B 分析用数据集
# 与预注册 §3–§4 对齐。仅做整形,不含任何 DIF/差距估计,可在冻结前运行。
# =====================================================================
suppressMessages({library(DBI); library(RSQLite); library(dplyr)})

DB   <- Sys.getenv("CESD_DB", "D:/clinicdatabase/SQLitedatabase/cesd_analysis.db")
OUT  <- Sys.getenv("CESD_OUT", "./data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
con  <- dbConnect(SQLite(), DB)

COMMON6 <- c("depres","effort","sleepr","whappy","flone","going")
CH_ITEMS <- c("depres","effort","sleepr","whappy","flone","going",
              "bother","mindts","fhope","fear")   # CHARLS 10, 0-3
EN_ITEMS <- c("depres","effort","sleepr","whappy","flone","fsad","going","enlife") # 8, 0/1

charls <- dbReadTable(con, "charls_cesd_items_long")
elsa   <- dbReadTable(con, "elsa_cesd_items_long")
hrs    <- dbReadTable(con, "hrs_cesd_items_long")
dbDisconnect(con)

## ---- CHARLS: w4 主分析(10 条目 0-3, 抑郁方向) + w1-4 纵向 ----
ch_d <- paste0(CH_ITEMS, "_d")
charls <- charls %>% mutate(across(all_of(c(ch_d,"cesd","agey","ragender","raeducl")), as.numeric))
ch_ok  <- function(d) d %>%
  filter(agey >= 60, ragender %in% c(1,2), !is.na(raeducl),
         if_all(all_of(ch_d), ~ !is.na(.))) %>%
  mutate(female = as.integer(ragender == 2),
         edu_low = as.integer(raeducl <= 1))

charls_w4  <- ch_ok(filter(charls, wave == 4))
charls_all <- ch_ok(charls)   # w1-4 纵向
write.csv(charls_w4,  file.path(OUT,"charls_w4.csv"),  row.names = FALSE)
write.csv(charls_all, file.path(OUT,"charls_w1to4.csv"), row.names = FALSE)

## ---- 三库共同 6 条目二分(cohort x sex), 各库取 >=60 完整个案最多之波 ----
prep_common <- function(d, id, items_d = paste0(COMMON6,"_d")){
  d <- d %>% mutate(across(all_of(c(items_d,"agey","ragender")), as.numeric)) %>%
    filter(agey >= 60, ragender %in% c(1,2), if_all(all_of(items_d), ~ !is.na(.)))
  w <- d %>% count(wave) %>% slice_max(n, n = 1) %>% pull(wave)
  d <- d %>% filter(wave == w) %>% mutate(female = as.integer(ragender == 2))
  # 二分: 抑郁方向 >=1 视为症状出现(CHARLS 0-3 -> 0/1; 英美已 0/1)
  for (b in COMMON6) d[[b]] <- as.integer(d[[paste0(b,"_d")]] >= 1)
  d %>% transmute(id = as.character(.data[[id]]), cohort, wave = w, female, age = agey,
                  depres, effort, sleepr, whappy, flone, going)
}
pooled <- bind_rows(prep_common(charls,"ID"),
                    prep_common(elsa,"idauniq"),
                    prep_common(hrs,"hhidpn"))
pooled$group6 <- as.integer(factor(paste(pooled$cohort, pooled$female)))  # 1..6 for Mplus
write.csv(pooled, file.path(OUT,"pooled_common6_binary.csv"), row.names = FALSE)

## ---- HRS + CIDI(外部效度) ----
hrs_d <- paste0(EN_ITEMS,"_d")
hrs2 <- hrs %>% mutate(across(any_of(c(hrs_d,"cesd","agey","ragender","cidi_mde3",
                                       "cidi_mde5","cidi_symp","cidi_dep","cidi_anh")), as.numeric)) %>%
  filter(agey >= 60, ragender %in% c(1,2), if_all(all_of(hrs_d), ~ !is.na(.))) %>%
  mutate(female = as.integer(ragender == 2))
write.csv(hrs2, file.path(OUT,"hrs_with_cidi.csv"), row.names = FALSE)

## ---- ELSA + 客观效标(握力/CRP) ----
en_d <- paste0(EN_ITEMS,"_d")
elsa2 <- elsa %>% mutate(across(any_of(c(en_d,"cesd","agey","ragender","grip_max",
                                         "crp","fibrinogen")), as.numeric)) %>%
  filter(agey >= 60, ragender %in% c(1,2), if_all(all_of(en_d), ~ !is.na(.))) %>%
  mutate(female = as.integer(ragender == 2))
write.csv(elsa2, file.path(OUT,"elsa_with_anchors.csv"), row.names = FALSE)

cat("Exported to", OUT, ":\n",
    " charls_w4.csv (N=", nrow(charls_w4), ")\n",
    " charls_w1to4.csv (N=", nrow(charls_all), ")\n",
    " pooled_common6_binary.csv (N=", nrow(pooled), ")\n",
    " hrs_with_cidi.csv (N=", nrow(hrs2), ")\n",
    " elsa_with_anchors.csv (N=", nrow(elsa2), ")\n", sep="")
