# DISCLOSURE: Due to rancher privacy agreements, precise geographic coordinates used in the manuscript analyses are not included in this repository. The version of the code shared here uses county-level centroids as placeholders. 
# As a result, the results and figure from this analysis differ slightly from those in the published manuscript. However, all primary findings and interpreations remain consistent.

library(httr)
library(jsonlite)
library(tidyverse)


metadata <- read_csv("metadata.csv") 

# Ensure Latitude and Longitude are numeric
metadata <- metadata %>%
  mutate(
    Latitude = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  )

# Deduplicate latitude and longitude by Location
unique_locations <- metadata %>%
  group_by(Location) %>%
  summarize(
    Latitude = mean(Latitude, na.rm = TRUE),
    Longitude = mean(Longitude, na.rm = TRUE),
    .groups = "drop"
  )


fetch_power_data <- function(lat, lon, start_date, end_date) {
  base_url <- "https://power.larc.nasa.gov/api/temporal/daily/point"
  
  # Make API request
  response <- GET(
    url = base_url,
    query = list(
      parameters = "T2M,T2M_MAX,T2M_MIN,PRECTOTCORR,ALLSKY_SFC_SW_DWN",  # Adjusted parameters
      community = "AG",
      longitude = lon,
      latitude = lat,
      start = gsub("-", "", start_date),  # Format dates as YYYYMMDD
      end = gsub("-", "", end_date),
      format = "JSON"
    )
  )
  
  # Handle API response
  if (http_status(response)$category == "Success") {
    raw_content <- content(response, as = "text", encoding = "UTF-8")
    data <- fromJSON(raw_content)$properties$parameter
    
    # Convert to data frame
    df <- data.frame(
      Date = as.Date(names(data$T2M), format = "%Y%m%d"),
      Temperature = unlist(data$T2M),
      Max_Temperature = unlist(data$T2M_MAX),
      Min_Temperature = unlist(data$T2M_MIN),
      Precipitation = unlist(data$PRECTOTCORR),
      Radiation = unlist(data$ALLSKY_SFC_SW_DWN),
      Latitude = lat,
      Longitude = lon
    )
    
    return(df)
  } else {
    warning("Failed to fetch data for site at (", lat, ", ", lon, "): ", http_status(response)$message)
    return(NULL)
  }
}

# Fetch daily data for all sites
daily_climate_data <- unique_locations %>%
  rowwise() %>%
  mutate(
    Climate_Data = list(fetch_power_data(Latitude, Longitude, "2023-01-01", "2023-12-31"))
  ) %>%
  unnest(cols = c(Climate_Data), names_sep = "_")  # Add a separator to disambiguate columns



# Define the Hargreaves function
calculate_pet <- function(t_mean, t_max, t_min, solar_radiation) {
  PET <- 0.0023 * solar_radiation * (t_mean + 17.8) * sqrt(t_max - t_min)
  return(PET)
}

# Calculate daily and annual climate metrics by Location
mean_annuals <- daily_climate_data %>%
  group_by(Location) %>%
  summarize(
    mean_annual_temp = mean(Climate_Data_Temperature, na.rm = TRUE),      # °C
    annual_high_temp = mean(Climate_Data_Max_Temperature, na.rm = TRUE),  # °C
    annual_low_temp = mean(Climate_Data_Min_Temperature, na.rm = TRUE),   # °C
    mean_daily_precip = mean(Climate_Data_Precipitation, na.rm = TRUE),   # mm/day
    annual_precip = sum(Climate_Data_Precipitation, na.rm = TRUE),        # mm/year
    radiation = mean(Climate_Data_Radiation, na.rm = TRUE)                # MJ/m2/day
  ) %>%
  mutate(
    PET_mm_day = calculate_pet(
      t_mean = mean_annual_temp,
      t_max = annual_high_temp,
      t_min = annual_low_temp,
      solar_radiation = radiation
    ),                                  # mm/day
    PET_mm_year = PET_mm_day * 365,     # mm/year
    aridity_index = annual_precip / PET_mm_year
  ) %>%
  ungroup()



# View the result
print(mean_annuals)


