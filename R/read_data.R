#' set directory for input files
#'
#' @return string specifying absolute path to directory containing input files
get_input_dir <- function() {
  "../data/"
}

#' get ids in lymphocyte study
#' 
#' @param threeDL2_neg logical. if TRUE, only get ids where 3DL2 is unlicensed
#' @return character vector of ids
#' @export
get_ids <- function() {
  lymphocyte_data <- read.csv(paste0(get_input_dir(), "Labelling_Study_Participant_Lymphocyte_Data.csv"),
                              stringsAsFactors = TRUE)
  ids <- levels(lymphocyte_data$id)
  ids
}

#' get cell populations and licensing statuses for a given id in a lymphocyte study
#' 
#' @param id string: participant id
#' @return data frame with five columns: cells, lic.status, cells_lic.status which
#' is the first two pasted together, disease ("HCV", "HIV-1", "HTLV-1" or "Control"), 
#' functional_iKIR_count(number of unique license genes)
#' @import dplyr
#' @export
get_id_data <- function(id) {
  kir_data <- read_kir_data(id)
  read_lymphocyte_data(id) %>%
    filter(time == 0) %>%
    select(id, cells, lic.status) %>%
    mutate(cells_lic.status = paste(cells, lic.status))
}

#' wrapper to only get cells_lic.status column for id data
#' @inheritParams get_id_data
#' @param CD8_only logical.  if TRUE, include CD8+ T cells only
#' @param NK_only logical. if TRUE, include NK cells only
#' @return character vector
#' @import dplyr
#' @export
get_cells_lic.status <- function(id, CD8_only = FALSE, NK_only = FALSE) {
  get_id_data(id) %>%
    filter_cell_data(CD8_only = CD8_only, NK_only = NK_only) %>%
    pull(cells_lic.status)
}

#' read in granulocyte or monocyte data
#' 
#' @param ids character vector. partiipant ids
#' @param cells "Granulocytes", "Monocytes" or "all".  Read data for these cell types
#' @return data frame with four columns: id, time, Lb, cells
#' @export
read_gm_data <- function(ids, cells) {
  stopifnot(cells %in% c("Granulocytes", "Monocytes", "all"))
  data_filename <- paste0(get_input_dir(), "Labelling_Study_Participant_Mono_and_Granulocyte_Data.csv")
  gm_data <- read.csv(data_filename, header = TRUE, stringsAsFactors = FALSE)
  if(cells != "all") {
    gm_data <- gm_data[gm_data$cells == cells,]
  }
  gm_data <- subset(gm_data, select = -c(weeks, data.fac))
  colnames(gm_data) <- c("id", "cells", "time", "Lb")
  stopifnot(all(ids %in% gm_data$id))
  gm_data <- gm_data[gm_data$id %in% ids,]
  gm_data
}


#' read in lymphocyte data
#' 
#' @param ids character vector. participant ids
#' @param cells_lic_status (optional) either data frame with two columns/list wih two elements: cells and lic.status,
#' or vector which is paste(cells, lic.status)
#' If specified, read those cells and license statuses only.  If not specified, read all.
#' @param drop_cells logical.  if TRUE, output three columns: id, time, L.
#' if FALSE, output five columns: id, time, L, cells, lic.status
#' @return data frame
#' @import dplyr
#' @export
read_lymphocyte_data <- function(ids, cells_lic.status, drop_cells = FALSE) {
  
  data_filename <- paste0(get_input_dir(), "Labelling_Study_Participant_Lymphocyte_Data.csv")
  lymphocyte_data <- read.csv(data_filename, header = TRUE, stringsAsFactors = FALSE) %>%
    filter(id %in% ids) %>%
    rename(time = days, L = frac.label) %>%
    mutate(cells_lic.status = paste(cells, lic.status))
  
  if(!missing(cells_lic.status)) {
    if(is.vector(cells_lic.status, mode = "character")) {
      cells_lic.status_paste <- cells_lic.status
    } else {
      cells_lic.status_paste <- paste(cells_lic.status$cells, cells_lic.status$lic.status)
    }
    if(drop_cells && length(cells_lic.status_paste) > 1) {
      warning("more than one cell population selected but dropping cell population column")
    }
    lymphocyte_data <- lymphocyte_data %>%
      filter(cells_lic.status %in% cells_lic.status_paste)
    missing_cells <- !(cells_lic.status_paste %in% lymphocyte_data$cells_lic.status)
    if(any(missing_cells))
      warning(paste("Cell types missing:", paste(cells_lic.status_paste[missing_cells], collapse = " ")))
  }
  
  keep_cols <- c("id", "time", "L")
  if(!drop_cells) {
    keep_cols <- c(keep_cols, "cells", "lic.status", "cells_lic.status")
  }
  
  lymphocyte_data <- lymphocyte_data %>%
    select(keep_cols)
  stopifnot(nrow(lymphocyte_data) > 0)
  lymphocyte_data
}

#' read in saliva data
#' 
#' @param ids (optional) character vector. if given, subset data for these ids only
#' @return data frame of saliva data
#' @export
read_saliva_data <- function(ids) {
  saliva_data_filename <- paste0(get_input_dir(), "Labelling_Study_Participant_Saliva_Data.csv")
  saliva_data <- read.csv(saliva_data_filename, header = TRUE, stringsAsFactors = FALSE)
  saliva_data <- subset(saliva_data, select = -c(weeks, sd))
  colnames(saliva_data) <- c("id", "time", "U")
  if(!missing(ids)) {
    stopifnot(all(ids %in% saliva_data$id))
    saliva_data <- saliva_data[saliva_data$id %in% ids,]
  }
  saliva_data
}

#' filter cell data frame based on conditions
#' 
#' to do: be able to have more than one filtering condition be true
#' @param cell_data data frame with cells column
#' @param CD8_only logical.  if TRUE, include CD8+ T cells only
#' @param NK_only logical. if TRUE, include NK cells only
#' @return data frame
#' @import dplyr
#' @export
filter_cell_data <- function(cell_data, 
                                   CD8_only = FALSE,
                                   NK_only = FALSE) {
  stopifnot(!(CD8_only && NK_only))
  if(CD8_only) {
    return(cell_data %>%
      filter(grepl("CD8", cells)))
  }
  if(NK_only) {
    return(cell_data %>%
       filter(grepl("NK", cells)))
  }
  warning("no filters applied")
  return(cell_data)
}

#' read in kir data
#' 
#' @param ids character vector. participant ids
#' @return data frame if length(ids) > 1, named vector otherwise
#' @import dplyr
#' @export
read_kir_data <- function(ids) {
  
  data_filename <- paste0(get_input_dir(), "Labelling_Study_Participant_KIR_Data.csv")
  kir_data <- read.csv(data_filename, header = TRUE, stringsAsFactors = FALSE) %>%
    filter(id %in% ids)
  if(nrow(kir_data) == 1) {
    kir_data <- unlist(kir_data)
  }
  kir_data
}

#' read data on the fraction of each cell population in each individual
#' 
#' @inheritParams read_kir_data
#' @return tibble with columns id, cells, lic.status, frac.cells
#' @importFrom dplyr %>% filter
#' @importFrom readr read_csv
#' @export
read_frac_data <- function(ids) {
  data_filename <- paste0(get_input_dir(), "Labelling_Study_KIR_Expression_Data.csv")
  frac_data <- read_csv(data_filename) %>%
    filter(id %in% ids)
  frac_data
}