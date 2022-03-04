# Use environment to store some path variables to use in different functions
ds_upload.globals <- new.env()

#' Uploading methylation data to the DataSHIELD backends
#'
#' @param upload do we need to upload the DataSHIELD backend
#' @param cohort_id cohort name from the dictonary
#' @param action action to be performed, can be 'populate', 'methyl' or 'all'
#' @param methyl_data_input_path path to the methylation data
#' @param covariate_data_input_path path to the covariate data to measure the age
#' @param dict_version version of the dictionary
#' @param data_version is the mean age (this needs to be a number)
#' @param data_format can be CSV or RData
#' @param dna_source can be 'cord_blood', 'peripheral_blood' or 'placenta'
#' @param norm_method can be RAW=0, BMIQ=1, DASEN=2, SWAN=3, SQN=4, RCP=5, NOOB=6, CPACOR=7, FUNNORM=8, OTHER=9 (needs to be a number)
#' @param run_mode can be NORMAL, NON_INTERACTIVE or TEST
#' @param database_name is the name of the data backend of DataSHIELD, default = opal_data
#' @param override_project override the generated project name
#'
#' @importFrom readr read_csv write_csv
#'
#' @examples
#'
#' \dontrun{
#' du.upload.methyl.clocks(
#'   cohort_id = 'gecko',
#'   methyl_data_input_path = "~/path-to-file",
#'   covariate_data_input_path = "~/path-to-file",
#'   dna_source = 'placenta',
#'   norm_method = 0,
#'   dict_version = "2_2",
#'   data_version = "1",
#'   override_project = "custom_project_1_0"
#' )
#' }
#'
#' @export
du.upload.methyl.clocks <- function(upload = TRUE, cohort_id, action = du.enum.action()$ALL, methyl_data_input_path = "", covariate_data_input_path = "", dict_version = '1_1', data_version = "16", data_format = du.enum.input.format()$CSV, dna_source = du.enum.dna.source()$CORD_BLOOD, norm_method = du.enum.norm.method()$RAW, run_mode = du.enum.run.mode()$NORMAL, database_name = "opal_data", override_project = NULL) {
  du.check.package.version()
  du.check.session(upload)

  message("######################################################")
  message("  Start upload methylation data into DataSHIELD backend")
  message("------------------------------------------------------")

  tryCatch(
    {
      workdirs <- du.create.temp.workdir()
      du.check.action(action)
      if(action == du.enum.action()$ALL | action == du.enum.action()$POPULATE) {
        du.dict.download(dict_version = dict_version, dict_kind = du.enum.dict.kind()$METHYL)
      }

      if (data_version == "" || !du.check.version(data_version)) {
        stop("No data version is specified or the data version does not match syntax: 'number*'! Program is terminated.")
      }
      
      if (action == du.enum.action()$ALL | action == du.enum.action()$POPULATE) {
        project <- du.populate(dict_version = dict_version, cohort_id = cohort_id, data_version = data_version, database_name, dict_kind = du.enum.dict.kind()$METHYL, override_project)
      }
      
      if (action == du.enum.action()$ALL | action == du.enum.action()$METHYL) {
        if (missing(methyl_data_input_path)) {
          input_path <- readline("- Specify input path (for your methylation data): ")
        } else if (missing(methyl_data_input_path)) {
          stop("No source file for methylation data specified, please specify your source for methylation data file")
        }
        if (missing(covariate_data_input_path)) {
          input_path <- readline("- Specify input path (for your covoriate data): ")
        } else if (missing(covariate_data_input_path)) {
          stop("No source file for covariate data specified, please specify your source for covariate data file")
        }
        if (missing(cohort_id)) {
          cohort_id <- readline("- Specify cohort identifier (e.g. dnbc): ")
        }
        if (cohort_id == "") {
          stop("No cohort identifier is specified! Program is terminated.")
        } else {
          if (!(cohort_id %in% du.enum.cohorts()) & run_mode != du.enum.run.mode()$TEST) {
            stop(
              "Cohort: [ ", cohort_id, " ] is not a known cohort in the netwprk Please choose from: [ ",
              paste(du.enum.cohorts(), collapse = ", "), " ]"
            )
          }
        }
        
        data_input_format <- data_format

        methyl_data <- du.generate.methyl.data(data_format, methyl_data_input_path, covariate_data_input_path, dna_source = dna_source, norm_method = norm_method)
        file_name <- paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "methyl_", dna_source, "_", data_version, "_", dna_source, ".csv")
        write_csv(methyl_data, paste0(getwd(), "/", file_name), na = "")
        if (upload) {
          if (ds_upload.globals$login_data$driver == du.enum.backends()$OPAL) {
            du.login(ds_upload.globals$login_data)
            du.opal.upload(du.enum.dict.kind()$METHYL, paste0(getwd(), "/", file_name))
          } else if (ds_upload.globals$login_data$driver == du.enum.backends()$ARMADILLO) {
            du.armadillo.import(project = project, data = methyl_data, dict_version = dict_version, data_version = data_version, dna_source = dna_source)
          }
        }
      }
    },
    finally = {
      du.clean.temp.workdir(upload, workdirs)
    }
  )
}

#' Generate the actual clocks.
#' 
#' Make sure you have an "Age" columns in the covariate data. It needs to be spelled exactly like that.
#' 
#' @param data_format can be CSV or Rdata
#' @param methyl_data_input_path input path of the raw methylation data
#' @param covariate_data_input_path input path of the covariate data (this can be used to determine the ages of the methylation clocks)
#' @param dna_source can be 'cord_blood', 'placenta' or 'peripheral_blood' 
#' @param norm_method can be RAW=0, BMIQ=1, DASEN=2, SWAN=3, SQN=4, RCP=5, NOOB=6, CPACOR=7, FUNNORM=8, OTHER=9 (needs to be a number)
#'
#' @importFrom readr read_csv
#' @importFrom tibble add_column
#' @importFrom RCurl url.exists
#'
#' @return the generated clocks with converted columns for child_id and the age_measured attached
#'
#' @noRd
du.generate.methyl.data <- function(data_format, methyl_data_input_path, covariate_data_input_path, dna_source = dna_source, norm_method) {
  requireNamespace("methylclock")
  
  if(data_format == du.enum.input.format()$CSV) {
    if(file.exists(methyl_data_input_path) & file.exists(covariate_data_input_path)) {
      methyl_data <- read_csv(methyl_data_input_path)
      covariate_data <- read_csv(covariate_data_input_path)
    } else {
      methyl_data <- source(methyl_data_input_path)
      covariate_data <- source(covariate_data_input_path)
    }
  } else {
    methyl_data <- load(methyl_data_input_path)
    covariate_data <- load(covariate_data_input_path) 
  }
  
  age <- covariate_data$Age

  if(dna_source == du.enum.dna.source()$CORD_BLOOD | dna_source == du.enum.dna.source()$PLACENTA) {
    message(paste0("* Generate: DNA methylation age: [ ", dna_source, " ]"))
    data <- methylclock::DNAmGA(x = methyl_data, age = age)
  } else {
    message(paste0("* Generate: DNA methylation gestational age: [ ", dna_source, " ]"))
    data <- methylclock::DNAmAge(x = methyl_data, age = age)
  }

  colnames(data)[colnames(data) == "id"] <- "child_id"

  row_id = seq.int(nrow(data))
  data <- add_column(data, row_id, .before = 1)
  data <- add_column(data, norm_method)

  return(data)
}
