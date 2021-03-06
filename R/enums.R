#' Cohorts participating in the cohort network
#'
#' @noRd
du.enum.cohorts <- function() {
  list(
    DNBC = "dnbc", GECKO = "gecko", ALSPAC = "alspac", GENR = "genr", MOBA = "moba", SWS = "sws", BIB = "bib", CHOP = "chop", ELFE = "elfe",
    EDEN = "eden", NINFEA = "ninfea", HBCS = "hbcs", INMA = "inma", ISGLOBAL = "isglobal", NBFC66 = "nfbc66", NBFC86 = "nfbc86", RAINE = "raine", RHEA = "rhea",
    ABCD = "abcd", BISC = "bisc", ENVIRONAGE = "environage", KANC = "kanc", PELAGIE = "pelagie", SEPAGES = "sepages", TNG = "tng", HGS = "hgs", RECETOX = "recetox", 
    GENXXI = "genxxi", GENESIS="genesis"
  )
}

#' Supported table types
#'
#' @noRd
du.enum.table.types <- function() {
  list(METHYL = "methyl")
}

#' Supported input formats
#'
#' @noRd
du.enum.input.format <- function() {
  list(CSV = "CSV", STATA = "STATA", SPSS = "SPSS", SAS = "SAS", R = "R")
}

#' Actions that can be performed
#'
#' @noRd
du.enum.action <- function() {
  list(ALL = "all", METHYL = "methyl", POPULATE = "populate")
}

#' Dictionary kinds
#'
#' @noRd
du.enum.dict.kind <- function() {
  list(METHYL = "methyl")
}

#' Projects that are containing dictionaries. Repositories containing these dictionaries should be:
#'
#' - ds-dictionaries
#' - ds-beta-dictionaries
#'
#' @noRd
du.enum.projects <- function() {
  list(LIFECYCLE = "lifecycle-project")
}

#' Possible DataSHIELD backends
#'
#' @noRd
du.enum.backends <- function() {
  list(OPAL = "OpalDriver", ARMADILLO = "ArmadilloDriver")
}

#' Run modes in uploading data
#'
#' @noRd
du.enum.run.mode <- function() {
  list(NORMAL = "normal", NON_INTERACTIVE = "non_interactive", TEST = "test")
}

du.enum.dna.source <- function() {
  list(CORD_BLOOD = "cord_blood", PERIPHERLA_BLOOD = "peripheral_blood", PLACENTA="placenta")
}

du.enum.norm.method <- function() {
  list(RAW=0, BMIQ=1, DASEN=2, SWAN=3, SQN=4, RCP=5, NOOB=6, CPACOR=7, FUNNORM=8, OTHER=9)
}


