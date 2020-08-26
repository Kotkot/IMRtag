#' Functions for downloading the data using code developed by Einar Hjorleifsson
#' (https://github.com/fishvice/taggart)
#'
#'
#'

ir2d <-
  function(ir, useI = FALSE)
  {
    lat <- substring(ir, 1, 2)
    lat <- as.numeric(lat)
    lat <- (lat +71)/2 + 0.25
    lon1 <- substring(ir, 3, 3)
    lon1 <- toupper(lon1)
    lon1 <- match(lon1, LETTERS)
    if(!useI) lon1 <- ifelse(lon1 > 8, lon1 - 1, lon1)
    lon1 <- lon1 - 2
    lon2 <- substring(ir, 4)
    lon2 <- as.numeric(lon2)
    lon <- ifelse(lon1 < 0,
                  -44 + lon2 + 0.5,
                  -40 + 10*lon1 + lon2 + 0.5)
    data.frame(lat = lat, lon = lon)
  }

#' Catches
#'
#' @param species "mackerel" or "herring"
#' @param cn.standardized Boolean, if FALSE (default) retains variable names as
#' delivered by the webserver otherwise mri-standaridzed variable names are used
#'
#' @return A dataframe
#' @export
#'
tg_catches <- function(species = "mackerel", cn.standardized = FALSE) {

  if(length(species) > 1) {
    stop(message("Only one species can be specified, 'mackerel' or 'herring'"))
  }

  if(!any(species %in% c("herring", "mackerel"))) {
    stop(message("species has to be either 'mackerel' or 'herring'"))
  }

  d <-
    jsonlite::fromJSON(paste0("http://smartfishsvc.hi.no/api/data/Catches/", species[1])) %>%
    dplyr::as_tibble() %>%
    dplyr::select_all(tolower) %>%
    dplyr::mutate(processing_date = lubridate::ymd_hms(processing_date ),
                  cLon = ir2d(ices_rectangle)$lon,
                  cLat = ir2d(ices_rectangle)$lat,
                  pLon = ir2d(plant_ices_rectangle)$lon,
                  pLat = ir2d(plant_ices_rectangle)$lat,
                  nation = dplyr::case_when(nation == "Eire" ~ "IE",
                                            nation == "Norway" ~ "NO",
                                            nation == "Sweden" ~ "SE",
                                            nation == "GB" ~ "UK",
                                            TRUE ~ nation),
                  species = species[1])

  if(!cn.standardized) {
    return(d)
  } else {
    d %>%
      dplyr::rename(ices = ices_rectangle,
                    year = catchdate,
                    pid = reference_plant,
                    pname = plant_name,
                    pices = plant_ices_rectangle,
                    area = recatch_ices_rectangle,
                    pdate = processing_date) %>%
      return()
  }

}

#' Catches biology
#'
#' @param species "mackerel" or "herring"
#' @param cn.standardized Boolean, if FALSE (default) retains variable names as
#' delivered by the webserver otherwise mri-standaridzed variable names are used
#'
#' @return A dataframe
#' @export
#'
tg_catches_bio <- function(species = "mackerel", cn.standardized = FALSE) {

  if(length(species) > 1) {
    stop(message("Only one species can be specified, 'mackerel' or 'herring'"))
  }

  if(!any(species %in% c("herring", "mackerel"))) {
    stop(message("species has to be either 'mackerel' or 'herring'"))
  }

  d <-
    jsonlite::fromJSON(paste0("http://smartfishsvc.hi.no/api/data/BioRawdataCatches/", species[1])) %>%
    dplyr::as_tibble() %>%
    dplyr::select_all(tolower) %>%
    dplyr::mutate(catch_date = lubridate::ymd_hms(catch_date),
                  lon = ir2d(ices_rectangle)$lon,
                  lat = ir2d(ices_rectangle)$lat,
                  species = species[1]) %>%
    dplyr::rename(maturity = mauturity)

  if(!cn.standardized) {
    return(d)
  } else {
    return(d)
  }
}

#' Expeditions
#'
#' @param species "mackerel" or "herring"
#' @param cn.standardized Boolean, if FALSE (default) retains variable names as
#' delivered by the webserver otherwise mri-standaridzed variable names are used
#' @return A dataframe
#' @export
#'
tg_expeditions <- function(species = "mackerel", cn.standardized = FALSE) {

  if(length(species) > 1) {
    stop(message("Only one species can be specified, 'mackerel' or 'herring'"))
  }

  if(!any(species %in% c("herring", "mackerel"))) {
    stop(message("species has to be either 'mackerel' or 'herring'"))
  }

  d <-
    jsonlite::fromJSON(paste0("http://smartfishsvc.hi.no/api/data/expeditions/", species[1])) %>%
    dplyr::as_tibble() %>%
    dplyr::select_all(tolower) %>%
    dplyr::mutate(when = lubridate::ymd_hms(when),
                  relesedate = lubridate::ymd_hms(relesedate),
                  recapturedate = lubridate::ymd_hms(recapturedate),
                  lo = stringr::str_replace(lo, ",", ".") %>% as.numeric(),
                  la = stringr::str_replace(la, ",", ".") %>% as.numeric(),
                  length = length * 100,
                  species = species[1])

  if(!cn.standardized) {
    return(d)
  } else {
    d %>%
      dplyr::rename(tDate = relesedate,
                    rDate = recapturedate,
                    species = fish,
                    tLon = lo,
                    tLat = la,
                    tLength = length) %>%
      return()
  }
}

#' Expeditions biology
#'
#' @param species "mackerel" or "herring"
#' @return A dataframe
#' @export
#'
tg_expeditions_bio <- function(species = "mackerel") {

  if(length(species) > 1) {
    stop(message("Only one species can be specified, 'mackerel' or 'herring'"))
  }

  if(!any(species %in% c("herring", "mackerel"))) {
    stop(message("species has to be either 'mackerel' or 'herring'"))
  }

  jsonlite::fromJSON(paste0("http://smartfishsvc.hi.no/api/data/BioRawdataExpeditions/", species[1])) %>%
    dplyr::as_tibble() %>%
    dplyr::select_all(tolower) %>%
    dplyr::mutate(catch_date = lubridate::ymd_hms(catch_date),
                  species = species[1]) %>%
    dplyr::rename(lon = lo,
                  lat = la,
                  maturity = mauturity) %>%
    return()

}
