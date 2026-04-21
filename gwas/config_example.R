# Example configuration for GWAS scripts
#
# Copy to config_study.R and edit paths / phenotype mappings as needed.

cfg <- list(

  # Outcome-to-column mappings
  # Keys are outcome IDs used in output filenames.
  # Values are the phenotype column names present in cohort RDS files.
  outcomes = list(

    # Continuous outcomes without baseline adjustment
    continuous_no_baseline = list(
      "1" = "pheno1",
      "3" = "pheno3",
      "5" = "pheno5",
      "7" = "pheno7",
      "9" = "pheno9",
      "11" = "pheno11",
      "13" = "pheno13",
      "15" = "pheno15",
      "17" = "pheno17",
      "19" = "pheno19"
    ),

    # Continuous outcomes with baseline adjustment
    continuous_with_baseline = list(
      "2" = list(pheno = "pheno2", baseline = "BL_pheno_2"),
      "4" = list(pheno = "pheno4", baseline = "BL_pheno_4"),
      "6" = list(pheno = "pheno6", baseline = "BL_pheno_6"),
      "8" = list(pheno = "pheno8", baseline = "BL_pheno_8"),
      "10" = list(pheno = "pheno10", baseline = "BL_pheno_10"),
      "12" = list(pheno = "pheno12", baseline = "BL_pheno_12"),
      "14" = list(pheno = "pheno14", baseline = "BL_pheno_14"),
      "16" = list(pheno = "pheno16", baseline = "BL_pheno_16"),
      "18" = list(pheno = "pheno18", baseline = "BL_pheno_18"),
      "20" = list(pheno = "pheno20", baseline = "BL_pheno_20")
    ),

    # Binary / event outcomes analysed with logistic mixed models
    binary = list(
      "21" = "pheno_21",
      "22" = "pheno_22",
      "23" = "pheno_23",
      "26" = "pheno_26",
      "27" = "pheno_27",
      "28" = "pheno_28"
    )
  ),

  ancestry = list(
    mixed = list(
      genotype_prefix = "/path/to/qc/postimputation/mixed_gwas_prefix",
      pcair_outdir = "/path/to/gwas/mixed/pcair",
      grm_rds = "/path/to/gwas/mixed/pcair/grm.rds",
      full_cohort_dir = "/path/to/cohort/mixed_fullcohort",
      sex_stratified_dir = "/path/to/cohort/mixed_cohort"
    ),

    eur = list(
      genotype_prefix = "/path/to/qc/postimputation/eur_gwas_prefix",
      pcair_outdir = "/path/to/gwas/eur/pcair",
      grm_rds = "/path/to/gwas/eur/pcair/grm.rds",
      full_cohort_dir = "/path/to/cohort/eur_fullcohort",
      sex_stratified_dir = "/path/to/cohort/eur_cohort"
    )
  ),

  gwas = list(
    pc_count = 10,
    nperbatch = 200,
    missing_method = "omit",
    maf_range = c(0, 0.5)
  )
)
