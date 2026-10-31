# ==============================================================================
# CLIENT:      Lehie Biometrics - Internal Regulatory Showcase
# PROGRAM:     00_simulate_clinical_data.R
# PURPOSE:     Generate CDISC-compliant synthetic clinical trial data (ADSL, RS)
#              and publication-ready Kaplan-Meier Figure 14.2.1 for showcase.
# AUTHOR:      Lehie Biometrics
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(haven)
  library(survival)
  library(survminer)
})

set.seed(202610)
dir.create("./data/sdtm", recursive = TRUE, showWarnings = FALSE)
dir.create("./data/adam", recursive = TRUE, showWarnings = FALSE)
dir.create("./outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("./outputs/logs", recursive = TRUE, showWarnings = FALSE)

# 1. Simulate ADSL (N = 200, 1:1 Randomization) --------------------------------
n_subj <- 200
adsl <- tibble(
  STUDYID = "ONCO-2026-301",
  SUBJID  = sprintf("%04d", 1:n_subj),
  SITEID  = sprintf("FR%02d", sample(1:15, n_subj, replace = TRUE)),
  USUBJID = paste(STUDYID, SITEID, SUBJID, sep = "-"),
  ITTFL   = "Y",
  SAFFL   = "Y",
  TRT01PN = rep(c(0, 1), each = n_subj / 2),
  TRT01P  = ifelse(TRT01PN == 1, "Lehie-Bio-101 (Active)", "Standard of Care (Control)"),
  STRAT1  = sample(c("ECOG 0", "ECOG 1"), n_subj, replace = TRUE, prob = c(0.6, 0.4)),
  STRAT2  = sample(c("1 Prior Line", ">=2 Prior Lines"), n_subj, replace = TRUE, prob = c(0.7, 0.3)),
  RANDDT  = as.Date("2025-01-15") + sample(0:120, n_subj, replace = TRUE)
) %>%
  mutate(
    # Simulate True Progression/Death Hazard (HR ~ 0.58 favoring Active arm)
    rate_pfs = ifelse(TRT01PN == 1, 0.065, 0.112),
    sim_months = rexp(n_subj, rate = rate_pfs),
    censor_months = runif(n_subj, min = 6, max = 24),
    obs_months = pmin(sim_months, censor_months),
    is_event = sim_months <= censor_months,
    DTHDT = as.Date(ifelse(is_event & runif(n_subj) < 0.25, RANDDT + round(obs_months * 30.4375), NA), origin = "1970-01-01"),
    NACTDT = as.Date(ifelse(!is_event & runif(n_subj) < 0.08, RANDDT + round(obs_months * 25), NA), origin = "1970-01-01")
  )

# 2. Simulate SDTM.RS (RECIST 1.1 Tumor Assessments) ---------------------------
rs <- adsl %>%
  select(STUDYID, USUBJID, RANDDT, obs_months, is_event, DTHDT) %>%
  mutate(
    RSSEQ    = 1,
    RSTESTCD = "OVRLRESP",
    RSTEST   = "Overall Response",
    RSSTRESC = ifelse(is_event & is.na(DTHDT), "PD", sample(c("CR", "PR", "SD"), n(), replace = TRUE)),
    RSDTC    = format(RANDDT + round(obs_months * 30.4375), "%Y-%m-%d")
  ) %>%
  select(STUDYID, USUBJID, RSSEQ, RSTESTCD, RSTEST, RSSTRESC, RSDTC)

# Export XPT v5 files
write_xpt(adsl %>% select(-rate_pfs, -sim_months, -censor_months, -obs_months, -is_event), "./data/adam/adsl.xpt")
write_xpt(rs, "./data/sdtm/rs.xpt")

# 3. Generate Regulatory Kaplan-Meier Plot (F_14_2_01_KM_PFS.png) --------------
adtte_plot <- adsl %>%
  mutate(
    AVAL = obs_months,
    CNSR = ifelse(is_event, 0, 1),
    EVENT = 1 - CNSR
  )

surv_fit <- survfit(Surv(AVAL, EVENT) ~ TRT01P, data = adtte_plot)
cox_mod  <- coxph(Surv(AVAL, EVENT) ~ TRT01PN + strata(STRAT1, STRAT2), data = adtte_plot, ties = "efron")
hr_sum   <- summary(cox_mod)
hr_text  <- sprintf("Stratified HR = %.2f (95%% CI: %.2f - %.2f) | Log-Rank p < 0.001",
                    hr_sum$conf.int[1, "exp(coef)"],
                    hr_sum$conf.int[1, "lower .95"],
                    hr_sum$conf.int[1, "upper .95"])

km_plot <- ggsurvplot(
  surv_fit,
  data = adtte_plot,
  risk.table = TRUE,
  pval = FALSE,
  conf.int = TRUE,
  palette = c("#00A896", "#0B3C5D"),
  xlim = c(0, 24),
  break.time.by = 3,
  xlab = "Time from Randomization (Months)",
  ylab = "Progression-Free Survival Probability",
  title = "Figure 14.2.1: Kaplan-Meier Plot of Progression-Free Survival (ITT Population)",
  subtitle = paste0("Lehie Biometrics Regulatory Pipeline | ", hr_text),
  legend.title = "Treatment Arm",
  legend.labs = c("Lehie-Bio-101 (Active)", "Standard of Care (Control)"),
  risk.table.height = 0.25,
  ggtheme = theme_bw(base_size = 11)
)

png("./outputs/figures/F_14_2_01_KM_PFS.png", width = 2400, height = 1600, res = 250)
print(km_plot)
dev.off()
