# setup-and-data-prep.R

###### Setup ######

library(tidyverse)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_light())

library(readr)
library(readxl)
library(janitor)
library(lubridate)
library(broom)
library(lwgeom) # needed by sf
library(sf)
library(sfheaders)
#Linking to GEOS 3.10.2, GDAL 3.4.2, PROJ 8.2.1; sf_use_s2() is TRUE
sf_use_s2(FALSE) # since some city simple features are not valid on a spherical globe

library(tigris)
options(tigris_use_cache = TRUE)

library(units)
library(scales)
library(geosphere) # for bearing()
library(gt)
library(glue)
library(patchwork)
library(hexbin) # for stat_binhex()

my_caption <- "Data: NOAA Storm Events Database; plot: Daniel Moul"
carolina_blue <- "#4B9CD3" #"#62C6F2"

set.seed(2023)


###### Data prep ######

###### tornado data

d_raw_tmp <- bind_rows(
  read_csv("./data/raw/tornados_1970-01-01_to_1973-12-31.csv"),
  read_csv("./data/raw/tornados_1974-01-01_to_1993-12-31.csv"),
  x2 <- read_csv("./data/raw/tornados_1994-01-01_to_2003-12-31.csv") |>
    mutate(TOR_LENGTH = as.numeric(TOR_LENGTH),
           TOR_WIDTH = as.numeric(TOR_WIDTH),
           BEGIN_LAT = as.numeric(BEGIN_LAT),
           BEGIN_LON = as.numeric(BEGIN_LON)
    ) |>
    select(-ABSOLUTE_ROWNUMBER),
  read_csv("./data/raw/tornados_2004-01-01_to_2013-12-31.csv"),
  read_csv("./data/raw/tornados_2014-01-01_to_2023-04-30.csv")
) |>
  clean_names()  |>
  relocate(episode_id, .after = event_id) |>
  select(-starts_with("magnitude"), # no data
         -starts_with("absolute"), # not unique since combined two files
         -flood_cause, # no data
         -c(damage_crops_num, state_abbr, cz_timezone, cz_type, source, 
            begin_range, end_range, begin_azimuth, end_azimuth) # not used
  ) |>
  replace_na(list(tor_f_scale = "F0")) |> # assume they were small
  mutate(begin_datetime = mdy_hm(paste(begin_date, str_pad(begin_time, 
                                                           width = 4,
                                                           side = "left",
                                                           pad = "0"))
                                 ),
         end_datetime = mdy_hm(paste(end_date, str_pad(end_time, 
                                                width = 4,
                                                side = "left",
                                                pad = "0")
                                     )
                               ),
         year = year(begin_datetime),
         day_of_year = yday(begin_datetime),
         decade = as.factor(paste0(floor(year / 10) * 10, "s"))
  ) |>
  arrange(begin_datetime) |>
  mutate(idx = row_number(),
         event_id = as.character(event_id),
         episode_id = as.character(episode_id),
         event_narrative = str_replace_all(event_narrative, "[|]", " "),
         episode_narrative = str_replace_all(episode_narrative, "[|]", " "),
         tor_f_scale = if_else(tor_f_scale == "EFU", "EF0", tor_f_scale), # assume unknowns were small
         tor_f_scale = parse_number(tor_f_scale),
         countyfp = str_pad(cz_fips,
                            width = 3,
                            side = "left",
                            pad = "0"),
         assoc_w_hurricane = str_detect(event_narrative, "(?i)hurricane|remnant") |
                             str_detect(episode_narrative, "(?i)hurricane|remnant")
  ) |>
  # fix known data error(s)
  mutate(damage_property_num = ifelse(event_id == "505635", 1.5e8, damage_property_num)
  ) |>
  select(-c(begin_date, begin_time, end_date, end_time, cz_fips))

d_raw <- d_raw_tmp |>
  filter(!is.na(begin_lon) & !is.na(end_lon)) |> # drop records that are missing start or end coordinates
  filter(event_id != "10090840")

year_start <- min(d_raw$year)

###### adjust for inflation (constant 2022 dollars)

d_inflation_ref <- suppressWarnings(
  read_xlsx("./data/inflation/us-inflation.xlsx",
            skip = 13) |>
    mutate(inflation_rate = as.numeric(inflation_rate)),
  classes = "warning")

d_inflation_adj <- 
  bind_rows(
    d_inflation_ref |>
      filter(year >= year_start,
             year <= 2022) |>
      mutate(adj_2022 = cpi[year == 2022] / cpi),
    # add 2023 using 222
    tibble(year = 2023)
  ) %>%
  fill(cpi:adj_2022, .direction = "downup") # add placeholders in first and last row


###### Get city, county and state boundary shapes
nc_cities <- st_read("./data/nc-onemap/NCDOT_City_Boundaries.geojson",
                     quiet = TRUE) |>
  st_transform(crs = "NAD83") |>
  clean_names() |>
  mutate(across(starts_with("county"), str_to_title)) |>
  select(-c(ends_with(c("date", "upper_case", "participate")), intersect_feature_id, changed, change_year)) |>
  filter(population_estimate >= 500)

nc_counties <- counties(state = "NC", cb = TRUE) |>
  clean_names() |>
  st_make_valid()

nc_state = st_union(nc_counties)

# How to convert [(begin_lat, begin_lon) - (end_lat, end_log)] to a linestring
# https://stackoverflow.com/questions/51918536/r-create-linestring-from-two-points-in-same-row-in-dataframe
# (but use tibble instead of data.table)

dta1_tmp <- d_raw[, c("begin_lon", "begin_lat", "idx")] |>
  rename(lon = begin_lon,
         lat = begin_lat)
dta2_tmp <- d_raw[, c("end_lon", "end_lat", "idx")]|>
  rename(lon = end_lon,
         lat = end_lat)

dta_tmp <- bind_rows(dta1_tmp, dta2_tmp) |>
  arrange(idx)

d_linestrings <- sfheaders::sf_linestring(
  obj = dta_tmp
  , x = "lon"
  , y = "lat"
  , linestring_id = "idx"
) |>
  st_set_crs("NAD83") # or "epsg:4326"

dta_beginpoints <- d_raw |>
  filter(!is.na(begin_lon)) |>
  st_as_sf(coords = c(x = "begin_lon", y = "begin_lat")) |>
  st_set_crs("NAD83") |>
  mutate(point_type = "begin")

dta_endpoints <- d_raw |>
  filter(!is.na(end_lon)) |>
  st_as_sf(coords = c(x = "end_lon", y = "end_lat")) |>
  st_set_crs("NAD83") |>
  mutate(point_type = "end")

dta_begin_end_points <- 
  bind_rows(
    dta_beginpoints,
    dta_endpoints
  ) |> 
  arrange(event_id, episode_id, begin_datetime) # TODO: confirm data supports the assumptions in this ordering

dta <- inner_join(d_linestrings, d_raw,
                  by = "idx") |>
  inner_join(d_inflation_adj,
             by = "year") |>
  select(-c(cpi, inflation_rate),
         -c(begin_lon, begin_lat, end_lon, end_lat)
         ) |>
  mutate(length_m = st_length(geometry),
         length_mi = set_units(length_m, "mi"),
         damage_property_mil = damage_property_num / 1e7,
         damage_property_mil_2022 = damage_property_mil * adj_2022) |>
  mutate(bear = if_else(as.numeric(length_mi) > 0,
                          map_dbl(geometry, ~first(bearing(.x, .y))),
                          NA),
         bear = if_else(bear < 0, bear + 360, bear)
  )

dta_non_mapping <- d_raw_tmp |>
  inner_join(d_inflation_adj,
             by = "year") |>
  mutate(damage_property_mil = damage_property_num / 1e7,
         damage_property_mil_2022 = damage_property_mil * adj_2022)

# rm(dta_tmp)
# rm(dta1_tmp)
# rm(dta2_tmp)
# rm(d_linestrings)

# test
# plot(dta |>
#        select(begin_datetime, tor_f_scale, injuries_direct, injuries_indirect, deaths_direct, deaths_indirect))

dta_for_datacheck <- dta |>
  st_drop_geometry() |>
  select(event_id, begin_datetime, begin_location, cz_name_str, tor_f_scale, tor_length, length_mi, year, day_of_year, decade) |>
  mutate(length_mi = drop_units(length_mi),
         length_diff = abs(tor_length - length_mi),
         pct_diff = abs(1 - length_mi / tor_length))

######

### The following was an attempt to disambiguate separate tornadoes from tornado segments.
### I wasn't successful in the time I was willing to spend
### I'm keeping the code commented out, since there are some useful techniques here

# get_endpoint <- function(g) {
#   # following example at https://gis.stackexchange.com/questions/251101/return-end-vertices-of-an-sf-linestring
#   if(length(g) > 0) {
#     st_line_sample(g, sample = 1) # first is 0, last is 1
#   } else
#   {
#     NA
#   }
# }
# 
# get_beginpoint <- function(g) {
#   if(length(g) > 0) {
#     st_line_sample(g, sample = 0) # first is 0, last is 1
#   } else
#   {
#     st_point() #st_is_empty()
#   }
# }
# 
# make_point <- function(g) {
#   if(!is.na(g)) {
#     st_point(unlist(g))
#   } else {
#     st_point() #st_is_empty()
#   }
# }
# 
# d_gaps <- dta |>
#   filter(!is.na(episode_id),
#          !is.na(event_id)) |>
#   select(event_id, episode_id) |>
#   arrange(episode_id, event_id) |> # TODO: confirm assumption that this order works
#   mutate(last_end = map(lag(geometry, default = NA), get_endpoint),
#          current_begin = map(geometry, get_beginpoint),
#          .by = episode_id) |>
#   filter(!is.na(last_end)) |>
#   st_cast(c(last_end, current_begin), to = "LINESTRING")  |>
#   mutate(prior_event_distance_m = st_length(geometry),
#          prior_event_distance_mi = drop_units(set_units(prior_event_distance_m, "mi"))
#          ) |>
#   st_drop_geometry() |>
#   select(-c(last_end, current_begin))
# 
# dta_with_gap_data <- dta |>
#   left_join(d_gaps,
#             by = c("episode_id", "event_id"))

