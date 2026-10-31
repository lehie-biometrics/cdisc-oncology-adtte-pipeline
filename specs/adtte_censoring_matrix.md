### CDISC ADaM Specification — `ADTTE` Parameter & FDA Censoring Hierarchy

**Study:** `ONCO-2026-301` | **Parameter:** `PFS` (Progression-Free Survival) | **Author:** Lehie Biometrics

| Situation / Clinical Scenario | `ADT` (Analysis Date) | `CNSR` | `EVNTDESC` | `CNSDTDSC` | `SRCDOM.SRCVAR` |
| :--- | :--- | :---: | :--- | :--- | :--- |
| **1. No baseline or post-baseline adequate assessment** and no death | `ADSL.RANDDT` | `2` | No Baseline or Post-Baseline Assessment | Randomization Date | `ADSL.RANDDT` |
| **2. New anti-cancer therapy (`NACTDT`)** initiated prior to documented PD or death | `LASTTUDT` prior to `NACTDT` | `3` | New Anti-Cancer Therapy Initiated Prior to Event | Last Adequate Assessment Prior to NACT | `RS.RSDTC` |
| **3a. Documented RECIST 1.1 Progression (`PD`)** prior to death or `NACTDT` | `PDDT` | `0` | Disease Progression per RECIST 1.1 | *Blank* | `RS.RSDTC` |
| **3b. Death from any cause** without prior documented `PD` or `NACTDT` | `ADSL.DTHDT` | `0` | Death from Any Cause | *Blank* | `ADSL.DTHDT` |
| **4. Alive and progression-free** at clinical cut-off date | `LASTTUDT` | `1` | Alive and Progression-Free at Data Cut-Off | Last Adequate Tumor Assessment | `RS.RSDTC` |
