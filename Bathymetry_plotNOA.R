
library(tidyverse)
library(marmap)
library(ggplot2)
library(scales)
library(grid)
library(maps)

# -----------------------------
# 1. Load sampling data
# -----------------------------
depth_coordinates <- read.csv(
  "/home/lvilleg/Downloads/depth_coordinates.csv",
  stringsAsFactors = FALSE
)

samples <- depth_coordinates %>%
  rename(
    depth_m = depth,
    lon     = long
  ) %>%
  mutate(
    depth_m = as.numeric(depth_m),
    lat     = as.numeric(lat),
    lon     = as.numeric(lon)
  ) %>%
  filter(!is.na(lon), !is.na(lat))

# -----------------------------
# 2. Map extent
# -----------------------------
lon_min <- 34.20
lon_max <- 35.20
lat_min <- 27.75
lat_max <- 29.75

# -----------------------------
# 3. Download bathymetry
# -----------------------------
bathy <- getNOAA.bathy(
  lon1 = lon_min,
  lon2 = lon_max,
  lat1 = lat_min,
  lat2 = lat_max,
  resolution = 0.25,
  keep = TRUE
)

# -----------------------------
# 4. Convert bathymetry to tidy table
# -----------------------------
bathy_df <- marmap::as.xyz(bathy) %>%
  as_tibble() %>%
  rename(
    lon = V1,
    lat = V2,
    depth_bathy_m = V3
  )

# -----------------------------
# 5. Split sea and land from bathymetry
# -----------------------------
sea_df <- bathy_df %>%
  filter(depth_bathy_m <= 0)

land_bathy_df <- bathy_df %>%
  filter(depth_bathy_m > 0)

sea_min <- min(sea_df$depth_bathy_m, na.rm = TRUE)

# -----------------------------
# 6. Extract seafloor depth at sample points
# -----------------------------
sample_bathy <- get.depth(
  bathy,
  x = samples$lon,
  y = samples$lat,
  locator = FALSE
)

samples$seafloor_depth_m <- sample_bathy$depth

write.csv(
  samples,
  "samples_with_bathymetry_wide.csv",
  row.names = FALSE
)

# -----------------------------
# 7. Country outlines from maps package
# -----------------------------
world_map <- map_data("world")

land_map <- world_map %>%
  filter(region %in% c("Egypt", "Saudi Arabia", "Jordan", "Israel")) %>%
  filter(
    long >= lon_min - 1,
    long <= lon_max + 1,
    lat  >= lat_min - 1,
    lat  <= lat_max + 1
  )

# -----------------------------
# 8. Country label positions
# -----------------------------
country_labels <- tibble(
  name = c("Egypt", "Saudi Arabia", "Jordan", "Israel"),
  lon  = c(34.40, 35.03, 35.03, 34.95),
  lat  = c(28.88, 28.98, 29.45, 29.58)
)

# -----------------------------
# 9. Contour levels
# -----------------------------
minor_breaks <- c(-50, -100, -200, -400, -600, -800,
                  -1000, -1200, -1400, -1600)
major_breaks <- c(-200, -600, -1000, -1400)

minor_breaks <- minor_breaks[
  minor_breaks >= sea_min & minor_breaks <= 0
]

major_breaks <- major_breaks[
  major_breaks >= sea_min & major_breaks <= 0
]

# -----------------------------
# 10. Colors
# -----------------------------
water_cols <- rev(c(
  "#f2f7fb",
  "#d6ecf5",
  "#a9dff0",
  "#76cde5",
  "#41b6c4",
  "#1d91c0",
  "#0868ac",
  "#084081"
))

land_fill <- "antiquewhite2"
border_col <- "#9a8a78"
coast_col  <- "#6f5c49"
label_col  <- "#5a4635"

# -----------------------------
# 11. Build plot
# -----------------------------
p <- ggplot() +
  # Land fill from bathymetry
  geom_raster(
    data = land_bathy_df,
    aes(x = lon, y = lat),
    fill = land_fill
  ) +
  # Water bathymetry
  geom_raster(
    data = sea_df,
    aes(x = lon, y = lat, fill = depth_bathy_m)
  ) +
  # Minor bathymetric contours
  geom_contour(
    data = sea_df,
    aes(x = lon, y = lat, z = depth_bathy_m),
    breaks = minor_breaks,
    colour = "black",
    linewidth = 0.10,
    alpha = 0.30
  ) +
  # Major bathymetric contours
  geom_contour(
    data = sea_df,
    aes(x = lon, y = lat, z = depth_bathy_m),
    breaks = major_breaks,
    colour = "black",
    linewidth = 0.18,
    alpha = 0.55
  ) +
  # Subtle country borders
  geom_path(
    data = land_map,
    aes(x = long, y = lat, group = group),
    colour = border_col,
    linewidth = 0.20,
    alpha = 0.80
  ) +
  # Coastline on top
  geom_contour(
    data = bathy_df,
    aes(x = lon, y = lat, z = depth_bathy_m),
    breaks = 0,
    colour = coast_col,
    linewidth = 0.35,
    alpha = 0.95
  ) +
  # Sample points
  geom_point(
    data = samples,
    aes(x = lon, y = lat),
    shape = 21,
    size = 2.8,
    stroke = 0.35,
    fill = "white",
    colour = "black"
  ) +
  # Country labels
  geom_text(
    data = country_labels,
    aes(x = lon, y = lat, label = name),
    size = 3.0,
    fontface = "bold",
    colour = label_col
  ) +
  coord_fixed(
    xlim = c(lon_min, lon_max),
    ylim = c(lat_min, lat_max),
    expand = FALSE
  ) +
  scale_fill_gradientn(
    colours = water_cols,
    limits = c(sea_min, 0),
    oob = squish,
    name = "Depth (m)",
    breaks = pretty(c(sea_min, 0), n = 6)
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(colour = "black", size = 9),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(2.0, "cm"),
    legend.position = "right",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(p)
