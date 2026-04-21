# ============================================================
# GWAS Configuration File
# ============================================================
#
# Copy this file to:
#   config_study.R
#
# Then update ALL paths below to match your environment.
#
# This config controls:
#   - outcome ↔ phenotype column mapping
#   - ancestry-specific inputs
#   - cohort definitions (explicit arm + sex)
#   - GWAS settings
#
# ============================================================

cfg <- list(

  # ============================================================
  # Outcome definitions
  # ============================================================

  outcomes = list(

    # Continuous outcomes WITHOUT baseline adjustment
    continuous_no_baseline = list(
      "1"  = "pheno1",
      "3"  = "pheno3",
      "5"  = "pheno5",
      "7"  = "pheno7",
      "9"  = "pheno9",
      "11" = "pheno11",
      "13" = "pheno13",
      "15" = "pheno15",
      "17" = "pheno17",
      "19" = "pheno19"
    ),

    # Continuous outcomes WITH baseline adjustment
    continuous_with_baseline = list(
      "2"  = list(pheno = "pheno2",  baseline = "BL_pheno_2"),
      "4"  = list(pheno = "pheno4",  baseline = "BL_pheno_4"),
      "6"  = list(pheno = "pheno6",  baseline = "BL_pheno_6"),
      "8"  = list(pheno = "pheno8",  baseline = "BL_pheno_8"),
      "10" = list(pheno = "pheno10", baseline = "BL_pheno_10"),
      "12" = list(pheno = "pheno12", baseline = "BL_pheno_12"),
      "14" = list(pheno = "pheno14", baseline = "BL_pheno_14"),
      "16" = list(pheno = "pheno16", baseline = "BL_pheno_16"),
      "18" = list(pheno = "pheno18", baseline = "BL_pheno_18"),
      "20" = list(pheno = "pheno20", baseline = "BL_pheno_20")
    ),

    # Binary / event outcomes
    binary = list(
      "21" = "pheno_21",
      "22" = "pheno_22",
      "23" = "pheno_23",
      "26" = "pheno_26",
      "27" = "pheno_27",
      "28" = "pheno_28"
    )
  ),

  # ============================================================
  # Ancestry-specific configuration
  # ============================================================

  ancestry = list(

    # ----------------------------
    # Mixed ancestry (primary)
    # ----------------------------
    mixed = list(

      genotype_prefix = "/path/to/qc/postimputation/mixed_gwas_prefix",

      pcair_outdir = "/path/to/gwas/mixed/pcair",
      grm_rds      = "/path/to/gwas/mixed/pcair/grm.rds",

      # All-participant cohorts
      cohorts_all = list(

        list(
          label = "albi_all",
          file  = "/path/to/cohort/mixed_fullcohort/albi_all.rds",
          arm   = "I",
          sex   = "ALL"
        ),

        list(
          label = "placebo_all",
          file  = "/path/to/cohort/mixed_fullcohort/placebo_all.rds",
          arm   = "C",
          sex   = "ALL"
        )
      ),

      # Sex-stratified cohorts
      cohorts_sex = list(

        list(
          label = "albi_female",
          file  = "/path/to/cohort/mixed_cohort/albi_female.rds",
          arm   = "I",
          sex   = "F"
        ),

        list(
          label = "albi_male",
          file  = "/path/to/cohort/mixed_cohort/albi_male.rds",
          arm   = "I",
          sex   = "M"
        ),

        list(
          label = "placebo_female",
          file  = "/path/to/cohort/mixed_cohort/placebo_female.rds",
          arm   = "C",
          sex   = "F"
        ),

        list(
          label = "placebo_male",
          file  = "/path/to/cohort/mixed_cohort/placebo_male.rds",
          arm   = "C",
          sex   = "M"
        )
      )
    ),

    # ----------------------------
    # European ancestry (sensitivity)
    # ----------------------------
    eur = list(

      genotype_prefix = "/path/to/qc/postimputation/eur_gwas_prefix",

      pcair_outdir = "/path/to/gwas/eur/pcair",
      grm_rds      = "/path/to/gwas/eur/pcair/grm.rds",

      cohorts_all = list(

        list(
          label = "albi_all",
          file  = "/path/to/cohort/eur_fullcohort/albi_all.rds",
          arm   = "I",
          sex   = "ALL"
        ),

        list(
          label = "placebo_all",
          file  = "/path/to/cohort/eur_fullcohort/placebo_all.rds",
          arm   = "C",
          sex   = "ALL"
        )
      ),

      cohorts_sex = list(

        list(
          label = "albi_female",
          file  = "/path/to/cohort/eur_cohort/albi_female.rds",
          arm   = "I",
          sex   = "F"
        ),

        list(
          label = "albi_male",
          file  = "/path/to/cohort/eur_cohort/albi_male.rds",
          arm   = "I",
          sex   = "M"
        ),

        list(
          label = "placebo_female",
          file  = "/path/to/cohort/eur_cohort/placebo_female.rds",
          arm   = "C",
          sex   = "F"
        ),

        list(
          label = "placebo_male",
          file  = "/path/to/cohort/eur_cohort/placebo_male.rds",
          arm   = "C",
          sex   = "M"
        )
      )
    )
  ),

  # ============================================================
  # GWAS settings
  # ============================================================

  gwas = list(
    pc_count       = 10,
    nperbatch      = 200,
    missing_method = "omit",
    maf_range      = c(0, 0.5)
  )
)
