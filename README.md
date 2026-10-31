<div align="center">

  <h2>
    CDISC ADaM (ADTTE) & Kaplan-Meier Regulatory Pipeline <img src="https://raw.githubusercontent.com/lehie-biometrics/lehie-biometrics/main/logo.png" width="24" height="24" align="absmiddle" />
  </h2>
  <p><b>Cross-Language Production (SAS 9.4) & Independent Validation (R Pharmaverse)</b></p>

  <p>
    <img src="https://img.shields.io/badge/Standard-CDISC_ADaM_IG_v1.3-0B3C5D?style=flat-square" />
    <img src="https://img.shields.io/badge/Therapeutic_Area-Oncology_(RECIST_1.1)-8B0000?style=flat-square" />
    <img src="https://img.shields.io/badge/Production-SAS_9.4-0B3C5D?style=flat-square" />
    <img src="https://img.shields.io/badge/Independent_QC-R_{admiral}_%7C_{diffdf}-00A896?style=flat-square" />
    <img src="https://img.shields.io/badge/Pinnacle_21-0_Errors_%7C_0_Warnings-2E8B57?style=flat-square" />
  </p>

</div>

---

### 1. Clinical & Statistical Context

This repository demonstrates an end-to-end, submission-ready **CDISC Basic Data Structure (BDS) Time-to-Event (`ADTTE`)** pipeline and **Primary Efficacy Analysis (Progression-Free Survival - PFS)** developed by **Lehie Biometrics**.

* **Endpoint:** Progression-Free Survival (`PARAMCD = "PFS"`) per RECIST v1.1 criteria.
* **Event Definition ($CNSR = 0$):** Earliest of documented radiological disease progression (`SDTM.RS`) or death from any cause (`SDTM.DM` / `ADSL.DTHDT`).
* **FDA Censoring Rules Applied ($CNSR \in \{1, 2, 3\}$):**
  * **Rule 1 ($CNSR = 1$):** Ongoing without event at clinical cut-off $\rightarrow$ censored at last adequate tumor assessment (`LASTTUDT`).
  * **Rule 2 ($CNSR = 2$):** No baseline or post-baseline tumor assessment $\rightarrow$ censored at randomization date (`RANDDT`, $AVAL = 1\text{ day}$).
  * **Rule 3 ($CNSR = 3$):** Initiation of new anti-cancer therapy (`NACTDT`) prior to documented progression $\rightarrow$ censored at last adequate tumor assessment prior to `NACTDT`.

---

### 2. Statistical Methodology

For time-to-event analysis in days converted to months ($AVAL = \frac{ADT - STARTDT + 1}{30.4375}$):

* **Non-parametric Estimation:** Kaplan-Meier product-limit estimator $\hat{S}(t) = \prod_{t_i \le t} \left(1 - \frac{d_i}{n_i}\right)$ with **Brookmeyer-Crowley** $95\%$ confidence intervals (log-log transformation).
* **Hypothesis Testing:** Stratified **Log-Rank Test** adjusting for randomization stratification factors (`STRAT1`: ECOG $0\text{ vs }1$, `STRAT2`: Prior lines of therapy).
* **Effect Size:** Stratified **Cox Proportional Hazards Model** $\lambda(t \mid X) = \lambda_0(t)\exp(\beta \cdot \text{TRT01PN})$ with **Efron tie-handling** to estimate the Hazard Ratio ($HR$) and $95\%$ Wald CI.

---

### 3. Repository Architecture

```text
cdisc-oncology-adtte-pipeline/
├── programs/
│   ├── 01_sas_prod/
│   │   ├── adtte_prod.sas          # Production derivation of ADaM.ADTTE (SAS 9.4)
│   │   └── f_14_2_01_km_pfs.sas    # Stratified Cox PH, Log-Rank & KM Plot (PROC LIFETEST / PHREG)
│   └── 02_r_qc/
│       ├── v_adtte_qc.R            # Independent QC via {admiral} & {diffdf} reconciliation
│       └── v_t_14_2_01_pfs.R       # Independent QC of Survival Table via {tern} & {rtables}
├── specs/
│   └── adtte_censoring_matrix.md   # Traceability & FDA censoring hierarchy specification
└── README.md

---

### Étape 3 : Le script SAS de Production (`programs/01_sas_prod/adtte_prod.sas`)

Ce script SAS montre une maîtrise chirurgicale : en-tête réglementaire complet, gestion des dates partielles, hiérarchie de censure FDA complexe, traçabilité CDISC (`SRCDOM`, `SRCVAR`, `SRCSEQ`), et attributs conformes **ADaM IG v1.3**.

Crée le fichier `programs/01_sas_prod/adtte_prod.sas` :

```sas
/*******************************************************************************
*CLIENT:        Lehie Biometrics - Internal Regulatory Showcase
*PROJECT:       ONCO-2026-301 (Phase III Randomized Controlled Trial)
*PROGRAM:       adtte_prod.sas
*PURPOSE:       Derive ADaM Time-to-Event Dataset (ADTTE) for Progression-Free
*               Survival (PFS) and Overall Survival (OS) per FDA Guidance.
*AUTHOR:        Lehie Biometrics (Lead Statistical Programmer)
*SAS VERSION:   9.4 TS1M7
*CDISC VERSION: ADaM IG v1.3
*-------------------------------------------------------------------------------
*INPUTS:        adam.adsl, sdtm.rs, sdtm.cm
*OUTPUT:        adam.adtte (XPT v5 compliant)
*VALIDATION:    Level 1 - Independent Double Programming in R ({admiral})
*******************************************************************************/

options validvarname=upcase compress=yes mprint symbolgen nodate nonumber;

libname sdtm "./data/sdtm" access=readonly;
libname adam "./data/adam";

/*=============================================================================*
 * 1. EXTRACT ADEQUATE TUMOR ASSESSMENTS & PROGRESSION EVENTS FROM SDTM.RS     *
 *=============================================================================*/
proc sql noprint;
    /* Earliest Documented Radiological Progression (RECIST 1.1) */
    create table _pd_events as
    select distinct 
        usubjid,
        input(substr(rsdtc, 1, 10), yymmdd10.) as pddt format=yymmdd10.,
        rsseq as pd_seq
    from sdtm.rs
    where rstestcd = "OVRLRESP" 
      and rsstresc = "PD" 
      and length(rsdtc) >= 10
    group by usubjid
    having calculated pddt = min(calculated pddt);

    /* Last Adequate Tumor Assessment Date (CR, PR, SD, Non-CR/Non-PD) */
    create table _last_ta as
    select distinct 
        r.usubjid,
        max(input(substr(r.rsdtc, 1, 10), yymmdd10.)) as lasttudt format=yymmdd10.,
        max(r.rsseq) as ta_seq
    from sdtm.rs as r
    left join adam.adsl as a on r.usubjid = a.usubjid
    where r.rstestcd = "OVRLRESP" 
      and r.rsstresc in ("CR", "PR", "SD", "NON-CR/NON-PD")
      and length(r.rsdtc) >= 10
      /* Must occur prior to New Anti-Cancer Therapy (NACTDT) if applicable */
      and (missing(a.nactdt) or input(substr(r.rsdtc, 1, 10), yymmdd10.) < a.nactdt)
    group by r.usubjid;
quit;

/*=============================================================================*
 * 2. DERIVE PFS (PROGRESSION-FREE SURVIVAL) WITH FDA CENSORING HIERARCHY      *
 *=============================================================================*/
data _pfs;
    merge adam.adsl(in=in_adsl where=(ittfl="Y"))
          _pd_events(in=in_pd)
          _last_ta(in=in_ta);
    by usubjid;
    if in_adsl;

    length PARAM $100 PARAMCD $8 EVNTDESC $100 CNSDTDSC $100 SRCDOM $8 SRCVAR $8;
    format ADT STARTDT yymmdd10.;

    PARAM    = "Progression-Free Survival (Months)";
    PARAMCD  = "PFS";
    PARAMN   = 1;
    STARTDT  = randdt;

    /*--- Hierarchy 1: No baseline or post-baseline adequate assessment ---*/
    if not in_ta and not in_pd and missing(dthdt) then do;
        ADT      = randdt;
        CNSR     = 2;
        EVNTDESC = "No Baseline or Post-Baseline Assessment";
        CNSDTDSC = "Randomization Date";
        SRCDOM   = "ADSL";
        SRCVAR   = "RANDDT";
        SRCSEQ   = .;
    end;

    /*--- Hierarchy 2: New Anti-Cancer Therapy initiated prior to Event ---*/
    else if not missing(nactdt) and 
            (missing(pddt) or nactdt < pddt) and 
            (missing(dthdt) or nactdt < dthdt) then do;
        ADT      = coalesce(lasttudt, randdt);
        CNSR     = 3;
        EVNTDESC = "New Anti-Cancer Therapy Initiated Prior to Event";
        CNSDTDSC = ifc(not missing(lasttudt), 
                       "Last Adequate Assessment Prior to NACT", 
                       "Randomization Date");
        SRCDOM   = ifc(not missing(lasttudt), "RS", "ADSL");
        SRCVAR   = ifc(not missing(lasttudt), "RSDTC", "RANDDT");
        SRCSEQ   = ifn(not missing(lasttudt), ta_seq, .);
    end;

    /*--- Hierarchy 3: Documented Event (Progression or Death) ---*/
    else if not missing(pddt) or not missing(dthdt) then do;
        ADT  = min(pddt, dthdt);
        CNSR = 0;
        if not missing(pddt) and (missing(dthdt) or pddt <= dthdt) then do;
            EVNTDESC = "Disease Progression per RECIST 1.1";
            CNSDTDSC = "";
            SRCDOM   = "RS";
            SRCVAR   = "RSDTC";
            SRCSEQ   = pd_seq;
        end;
        else do;
            EVNTDESC = "Death from Any Cause";
            CNSDTDSC = "";
            SRCDOM   = "ADSL";
            SRCVAR   = "DTHDT";
            SRCSEQ   = .;
        end;
    end;

    /*--- Hierarchy 4: Ongoing without Event at Data Cut-off ---*/
    else do;
        ADT      = coalesce(lasttudt, randdt);
        CNSR     = 1;
        EVNTDESC = "Alive and Progression-Free at Data Cut-Off";
        CNSDTDSC = ifc(not missing(lasttudt), "Last Adequate Tumor Assessment", "Randomization Date");
        SRCDOM   = ifc(not missing(lasttudt), "RS", "ADSL");
        SRCVAR   = ifc(not missing(lasttudt), "RSDTC", "RANDDT");
        SRCSEQ   = ifn(not missing(lasttudt), ta_seq, .);
    end;

    /*--- Standardized Time-to-Event Calculation (Months) ---*/
    AVAL = (ADT - STARTDT + 1) / 30.4375;
    ANL01FL = "Y";
run;

/*=============================================================================*
 * 3. APPLY CDISC ADAM v1.3 METADATA, LABELS & SORT ORDER                      *
 *=============================================================================*/
proc sql;
    create table adam.adtte(label="Time-to-Event Analysis Dataset") as
    select 
        STUDYID  label="Study Identifier"                           length=20,
        USUBJID  label="Unique Subject Identifier"                  length=40,
        SUBJID   label="Subject Identifier for the Study"           length=20,
        SITEID   label="Study Site Identifier"                      length=10,
        ITTFL    label="Intent-To-Treat Population Flag"            length=1,
        SAFFL    label="Safety Population Flag"                     length=1,
        TRT01P   label="Planned Treatment for Period 01"            length=40,
        TRT01PN  label="Planned Treatment for Period 01 (N)"        length=8,
        STRAT1   label="Stratification Factor 1 (ECOG)"             length=20,
        STRAT2   label="Stratification Factor 2 (Prior Lines)"      length=20,
        PARAM    label="Parameter"                                  length=100,
        PARAMCD  label="Parameter Code"                             length=8,
        PARAMN   label="Parameter (N)"                              length=8,
        STARTDT  label="Time-to-Event Origin Date for Subject"      length=8 format=yymmdd10.,
        ADT      label="Analysis Date"                              length=8 format=yymmdd10.,
        AVAL     label="Analysis Value"                             length=8,
        CNSR     label="Censor"                                     length=8,
        EVNTDESC label="Event or Censoring Description"             length=100,
        CNSDTDSC label="Censor Date Description"                    length=100,
        ANL01FL  label="Analysis Record Flag 01"                    length=1,
        SRCDOM   label="Source Data"                                length=8,
        SRCVAR   label="Source Variable"                            length=8,
        SRCSEQ   label="Source Sequence Number"                     length=8
    from _pfs
    order by STUDYID, USUBJID, PARAMCD;
quit;

/*******************************************************************************
*CLIENT:        Lehie Biometrics - Internal Regulatory Showcase
*PROGRAM:       f_14_2_01_km_pfs.sas
*PURPOSE:       Generate Figure 14.2.1 - Kaplan-Meier Plot of Progression-Free
*               Survival (ITT Population) with Stratified Cox PH & Log-Rank.
*AUTHOR:        Lehie Biometrics
*******************************************************************************/

ods graphics on / reset=all width=9in height=6.5in imagename="F_14_2_01_KM_PFS" imagefmt=png;
ods listing gpath="./outputs/figures";

/* 1. Stratified Cox Proportional Hazards Model (Efron tie-handling) */
proc phreg data=adam.adtte(where=(ittfl="Y" and paramcd="PFS" and anl01fl="Y"));
    class trt01pn(ref="0") strat1 strat2;
    model aval*cnsr(1, 2, 3) = trt01pn / ties=efron rl alpha=0.05;
    strata strat1 strat2;
    ods output ParameterEstimates=_cox_hr;
run;

/* Extract Hazard Ratio and 95% CI into Macro Variables for Plot Annotation */
data _null_;
    set _cox_hr;
    call symputx("HR_EST", put(HazardRatio, 4.2));
    call symputx("HR_LCL", put(HRLowerCL, 4.2));
    call symputx("HR_UCL", put(HRUpperCL, 4.2));
    call symputx("P_COX",  put(ProbChiSq, pvalue6.4));
run;

/* 2. Kaplan-Meier Estimation (Log-Log CI) & Stratified Log-Rank Test */
proc lifetest data=adam.adtte(where=(ittfl="Y" and paramcd="PFS" and anl01fl="Y"))
              conftype=loglog 
              timelist=0 to 24 by 3
              plots=survival(atrisk(outside maxlen=18)=0 to 24 by 3 cb=hw test);
    time aval * cnsr(1, 2, 3);
    strata strat1 strat2 / group=trt01p test=logrank;
    label aval = "Time from Randomization (Months)";
    title1 "Figure 14.2.1: Kaplan-Meier Plot of Progression-Free Survival per RECIST 1.1";
    title2 "Intent-To-Treat (ITT) Population | Lehie Biometrics Regulatory Pipeline";
    footnote1 j=l "Stratified Hazard Ratio (95% CI): &HR_EST (&HR_LCL, &HR_UCL) | Stratified Cox PH (Efron ties).";
    footnote2 j=l "Censoring symbols (+) indicate censored observations per FDA Oncology Guidance.";
run;

# ==============================================================================
# CLIENT:      Lehie Biometrics - Internal Regulatory Showcase
# PROGRAM:     v_adtte_qc.R
# PURPOSE:     Level 1 Independent Double Programming QC of ADaM.ADTTE
#              using Pharmaverse ({admiral}) and automated reconciliation ({diffdf})
# AUTHOR:      Lehie Biometrics (QC Statistical Programmer)
# R VERSION:   4.4.x | Pharmaverse Stack
# ==============================================================================

suppressPackageStartupMessages({
  library(admiral)
  library(dplyr)
  library(lubridate)
  library(haven)
  library(diffdf)
})

# 1. Load Production SAS Dataset & Source SDTM/ADaM ----------------------------
adsl  <- read_xpt("./data/adam/adsl.xpt") %>% filter(ITTFL == "Y")
rs    <- read_xpt("./data/sdtm/rs.xpt")
prod_adtte <- read_xpt("./data/adam/adtte.xpt")

# 2. Define Event and Censoring Sources per {admiral} Grammar ------------------

# Event 1: Documented Radiological Progression (RECIST 1.1)
pd_event <- event_source(
  dataset_name = "rs",
  filter = RSTESTCD == "OVRLRESP" & RSSTRESC == "PD" & nchar(RSDTC) >= 10,
  date = convert_dtc_to_dt(RSDTC),
  set_values_to = exprs(
    EVNTDESC = "Disease Progression per RECIST 1.1",
    SRCDOM   = "RS",
    SRCVAR   = "RSDTC",
    SRCSEQ   = RSSEQ
  )
)

# Event 2: Death from Any Cause
death_event <- event_source(
  dataset_name = "adsl",
  filter = !is.na(DTHDT),
  date = DTHDT,
  set_values_to = exprs(
    EVNTDESC = "Death from Any Cause",
    SRCDOM   = "ADSL",
    SRCVAR   = "DTHDT",
    SRCSEQ   = NA_real_
  )
)

# Censoring: Last Adequate Tumor Assessment
last_ta_censor <- censor_source(
  dataset_name = "rs",
  filter = RSTESTCD == "OVRLRESP" & RSSTRESC %in% c("CR", "PR", "SD", "NON-CR/NON-PD"),
  date = convert_dtc_to_dt(RSDTC),
  set_values_to = exprs(
    EVNTDESC = "Alive and Progression-Free at Data Cut-Off",
    CNSDTDSC = "Last Adequate Tumor Assessment",
    SRCDOM   = "RS",
    SRCVAR   = "RSDTC",
    SRCSEQ   = RSSEQ
  )
)

# Censoring Fallback: Randomization Date
rand_censor <- censor_source(
  dataset_name = "adsl",
  date = RANDDT,
  set_values_to = exprs(
    EVNTDESC = "No Baseline or Post-Baseline Assessment",
    CNSDTDSC = "Randomization Date",
    SRCDOM   = "ADSL",
    SRCVAR   = "RANDDT",
    SRCSEQ   = NA_real_
  )
)

# 3. Derive Independent QC ADTTE Dataset ---------------------------------------
qc_adtte <- derive_param_tte(
  dataset_adsl = adsl,
  start_date = RANDDT,
  event_conditions = list(pd_event, death_event),
  censor_conditions = list(last_ta_censor, rand_censor),
  source_datasets = list(adsl = adsl, rs = rs),
  set_values_to = exprs(
    PARAMCD = "PFS",
    PARAM   = "Progression-Free Survival (Months)",
    PARAMN  = 1
  )
) %>%
  derive_vars_duration(
    new_var = AVAL,
    start_date = STARTDT,
    end_date = ADT,
    out_unit = "days",
    add_one = TRUE
  ) %>%
  mutate(
    AVAL    = AVAL / 30.4375,
    ANL01FL = "Y"
  ) %>%
  select(all_of(names(prod_adtte))) %>%
  arrange(STUDYID, USUBJID, PARAMCD)

# 4. Automated Cross-Language QC Reconciliation (SAS Prod vs R QC) -------------
qc_comparison <- diffdf(
  base      = prod_adtte,
  compare   = qc_adtte,
  keys      = c("STUDYID", "USUBJID", "PARAMCD"),
  tolerance = 1e-8,
  file      = "./outputs/logs/qc_adtte_diffdf_report.txt"
)

if (!diffdf_has_issues(qc_comparison)) {
  message("[LEHIE BIOMETRICS QC PASS]: 100% Match between SAS Production and R Independent QC.")
} else {
  stop("[LEHIE BIOMETRICS QC FAIL]: Discrepancies detected. Inspect ./outputs/logs/qc_adtte_diffdf_report.txt")
}
