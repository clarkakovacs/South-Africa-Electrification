rm(list = ls())
# ============================================================
# Build settlement-day panel from daily VIIRS GeoTIFF rasters
#
# Inputs
# - Settlements polygons (GeoPackage, full geometry): settlements.gpkg
# - Daily GeoTIFFs exported from GEE at 500 m, with 2 bands:
#     band1 = lit   (1 where valid AND above threshold, NA elsewhere)
#     band2 = valid (1 where valid, 0 elsewhere)
#
# Outputs
# - settlement_day.parquet: settlement_id x date panel with p_lit_sett, n_valid, coverage, etc.
# - (optional) settlement_day.gpkg: settlement geometries joined to panel (not recommended for long panels)
#
# Alternative options / tradeoffs
# - exact=TRUE vs exact=FALSE in terra::extract():
#     exact=TRUE weights boundary pixels by overlap fraction (more accurate, slower).
#     exact=FALSE counts pixels by centroid inclusion (faster, more bias for small polygons).
# - Coverage definition:
#     coverage = n_valid / n_max (pixel-count proxy; robust in raster workflow).
#     Alternative: compute area_obs / area_sett using exact weighted pixel areas (more complex).
# - Output storage:
#     Parquet is preferred for panel data; GPKG is for geometry snapshots or monthly aggregates.
# ============================================================

library(sf)      # read/write settlement polygons
library(dplyr)   # aggregation, joins
library(arrow)   # write_parquet
library(terra)   # raster reading + zonal extraction
library(stringr) # date parsing from filenames

sf::sf_use_s2(FALSE) # optional: avoids s2 issues with some polygons

# -----------------------------
# PATHS (EDIT)
# -----------------------------
BASE_PATH <- "/Users/nylan/Library/Mobile Documents/com~apple~CloudDocs/BSE/4. Term 2/23DM017 Geospatial Data Science/Project/1. Main"

# Settlements: FULL geometry for computation (do not use simplified file here)
SETT_GPKG <- file.path(BASE_PATH, "Map Data", "Settlements", "GPKG", "south_africa_dre_atlas_settlements_final.gpkg")

# Folder containing daily GeoTIFFs downloaded from Drive
# Each file corresponds to one day; filenames should include YYYY-MM-DD
RAST_DIR <- file.path(BASE_PATH, "gee_exports-500m_raster", "01")

# Output folder for panel + optional mapping file
OUT_DIR <- file.path(BASE_PATH, "Map Data", "settlement_day_outputs_rasters")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PARQ <- file.path(OUT_DIR, "settlement_day_final.parquet")
OUT_GPKG <- file.path(OUT_DIR, "settlement_month_final.gpkg")  # recommended to write monthly, not daily
OUT_LAYER <- "settlement_month"

# -----------------------------
# PARAMETERS (EDIT)
# -----------------------------

# Pixel resolution in meters (must match your GEE export scale)
PIXEL_M <- 500

# Minimum fraction of the settlement that must be "valid-observed" on a day
# (prevents interpreting low p_lit when most of the settlement was cloud-masked / missing)
MIN_COVERAGE <- 0.50
# Alternatives:
# - 0.25 for lenient (more data, more noise)
# - 0.60–0.80 for strict (less data, higher confidence)

# -----------------------------
# 1) READ SETTLEMENTS
# -----------------------------
sett <- st_read(SETT_GPKG, quiet = TRUE)

# Ensure a stable ID exists (your earlier preprocessing script should already create this)
if (!"settlement_id" %in% names(sett)) {
  stop("settlement_id not found in settlements.gpkg. Create a stable settlement_id at ingestion.")
}

# Clean invalid geometries (important for extraction and area calculations)
sett <- st_make_valid(sett)

# Convert to terra vector for fast extraction
# NOTE: We will reproject to raster CRS after loading a sample raster (below).
sett_v <- terra::vect(sett)

# -----------------------------
# 2) LIST DAILY RASTERS
# -----------------------------
tifs <- list.files(RAST_DIR, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
stopifnot(length(tifs) > 0)

# Helper: parse date from filename
# Expected: ..._YYYY-MM-DD.tif (adapt pattern if your naming differs)
get_date <- function(path) {
  m <- str_match(basename(path), "(\\d{4}-\\d{2}-\\d{2})")
  if (is.na(m[,2])) stop("Could not parse date from filename: ", basename(path))
  as.Date(m[,2])
}

# -----------------------------
# 3) ALIGN CRS (SETTLEMENTS -> RASTER CRS)
# -----------------------------
# Read one raster to get the CRS used in exports.
r0 <- terra::rast(tifs[1])

# If your exported GeoTIFF includes CRS, this will work. If not, you must set it manually.
if (is.na(terra::crs(r0))) {
  stop("Raster CRS is missing. Fix your GEE export to embed CRS or set it manually in R.")
}

# Reproject settlements to raster CRS once (avoid reprojecting inside the loop)
sett_v <- terra::project(sett_v, terra::crs(r0))

# -----------------------------
# 4) PRECOMPUTE N_MAX PER SETTLEMENT (MAX POSSIBLE PIXELS)
# -----------------------------
# n_max is approximated from settlement area / pixel area.
# Used to convert n_valid into a coverage proxy.
pixel_area_m2 <- PIXEL_M * PIXEL_M

# Compute area in a metric CRS.
# If raster CRS is geographic (degrees), compute areas in EPSG:3857 or a local UTM.
sett_area_m2 <- as.numeric(st_area(st_transform(sett, 3857)))

sett_area_tbl <- tibble(
  settlement_id = sett$settlement_id,
  area_m2 = sett_area_m2,
  n_max = pmax(1, ceiling(area_m2 / pixel_area_m2)) # at least 1 pixel
)

# -----------------------------
# 5) LOOP OVER DAYS: ZONAL STATS
# -----------------------------
results <- vector("list", length(tifs))

for (i in seq_along(tifs)) {
  f <- tifs[i]
  d <- get_date(f)
  cat("Processing", d, ":", basename(f), "\n")
  
  # Read raster for that day (2 bands expected)
  r <- terra::rast(f)
  
  # Enforce band naming assumptions:
  # - If your GEE export order differs, swap accordingly.
  # Alternatives:
  # - Use names(r) already embedded in the file (if set by exporter)
  if (terra::nlyr(r) < 2) stop("Expected 2 bands (lit, valid) in: ", basename(f))
  names(r) <- c("lit", "valid")
  
  # Extract zonal sums per settlement
  # - sum(valid) gives n_valid (valid pixels inside polygon; exact=TRUE yields fractional weights)
  # - sum(lit) gives n_lit (lit pixels inside polygon; lit is NA outside valid mask)
  ex_valid <- terra::extract(r[["valid"]], sett_v, fun = sum, na.rm = TRUE, exact = TRUE)
  ex_lit   <- terra::extract(r[["lit"]],   sett_v, fun = sum, na.rm = TRUE, exact = TRUE)
  
  # extract() returns a data.frame with an ID column referring to polygon order
  df_day <- tibble(
    settlement_id = sett$settlement_id,
    date = d,
    n_valid = ex_valid[[2]],
    n_lit   = ex_lit[[2]]
  ) %>%
    mutate(
      # p_lit_sett = lit share among valid pixels for that settlement-day
      # If n_valid is 0, the settlement was not observed (clouds/masks), so p_lit_sett is NA
      p_lit_sett = ifelse(n_valid > 0, n_lit / n_valid, NA_real_)
    )
  
  results[[i]] <- df_day
}

sett_day <- bind_rows(results)

# -----------------------------
# 6) ADD COVERAGE + QC FILTER
# -----------------------------
# Coverage proxy:
# - n_valid / n_max approximates "fraction of settlement observed by valid pixels"
# Alternative:
# - compute exact observed area using weighted pixel areas (more work, marginal gain).
sett_day <- sett_day %>%
  left_join(sett_area_tbl, by = "settlement_id") %>%
  mutate(
    coverage = pmin(1, n_valid / n_max)  # clamp to [0,1] for safety
  ) %>%
  filter(coverage >= MIN_COVERAGE)

# -----------------------------
# 7) SAVE SETTLEMENT-DAY PANEL (ATTRIBUTES ONLY)
# -----------------------------
arrow::write_parquet(sett_day, OUT_PARQ)
cat("Wrote:", OUT_PARQ, "\n")
cat("Rows:", nrow(sett_day), "\n")

sum(is.na(sett_month$mean_p_lit))

# -----------------------------
# 8) OPTIONAL: WRITE MONTHLY GEOMETRY FILE FOR MAPPING (RECOMMENDED)
# -----------------------------
# Writing geometry for every day creates an enormous GPKG and is usually unnecessary.
# Better: compute monthly aggregates and map those.

sett_month <- sett_day %>%
  group_by(settlement_id) %>%
  summarise(
    mean_p_lit = mean(p_lit_sett, na.rm = TRUE),
    sd_p_lit   = ifelse(sum(!is.na(p_lit_sett)) >= 2, sd(p_lit_sett, na.rm = TRUE), NA_real_),
    drop_frac  = mean(p_lit_sett < 0.2, na.rm = TRUE),  # alternative threshold
    n_days_obs = sum(!is.na(p_lit_sett)),
    mean_cov   = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  )

sett_month_sf <- sett %>%
  select(settlement_id) %>%
  left_join(sett_month, by = "settlement_id")

st_write(
  sett_month_sf,
  OUT_GPKG,
  layer = OUT_LAYER,
  delete_layer = TRUE,
  quiet = TRUE
)



cat("Wrote monthly mapping layer:", OUT_GPKG, "layer =", OUT_LAYER, "\n")
message("Done.")
