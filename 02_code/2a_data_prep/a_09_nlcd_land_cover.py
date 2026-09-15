"""
a_09_nlcd_land_cover.py
Modal NLCD 2022 land cover for every cell of the 8-state fishnet.

Table 1 wants land cover for the spillover cells, and the old lookup only ever covered the
Minnesota subset - 73 of the 153 analysis cells, so every non-MN case cell came out NA. This
rebuilds it over the whole flyway grid.

Two traps worth knowing about:
  - fishnet_8state.geojson DECLARES EPSG:4326 but its coordinates are the projected EPSG:5070
    grid, so the CRS has to be overridden rather than transformed. Same gotcha as d_02.
  - land cover is constant within a stratum, so this is descriptive only. It cannot enter a
    case-crossover model as a covariate - it's perfectly collinear with the stratum intercept.

Reads windows per cell rather than the whole CONUS raster; 13,920 cells x ~111k pixels is
manageable that way and 3 billion pixels in memory is not.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds

base = Path(__file__).resolve().parents[2]
env_data = base / "01_data" / "1a_exposure_data" / "meteorological_data"
lc_dir = (base / "01_data" / "1c_supportive_datasets" / "land_cover_data"
          / "Annual_NLCD_LndChg_2022_CU_C1V2")
out_dir = base / "03_output" / "3b_model_output" / "dataframes"

NLCD = {11: "Open Water", 12: "Perennial Ice/Snow", 21: "Developed, Open Space",
        22: "Developed, Low Intensity", 23: "Developed, Medium Intensity",
        24: "Developed, High Intensity", 31: "Barren Land", 41: "Deciduous Forest",
        42: "Evergreen Forest", 43: "Mixed Forest", 52: "Shrub/Scrub",
        71: "Herbaceous", 81: "Pasture/Hay", 82: "Cultivated Crops",
        90: "Woody Wetlands", 95: "Emergent Herbaceous Wetlands"}

grid = gpd.read_file(env_data / "hpai_multistate_era5" / "fishnet_8state.geojson")
grid = grid.set_crs(5070, allow_override=True)      # declared 4326, actually 5070
print(f"fishnet: {len(grid):,} cells across {grid.state.nunique()} states")

tif = lc_dir / "Annual_NLCD_LndCov_2022_CU_C1V2.tif"
rows = []
with rasterio.open(tif) as src:
    print(f"raster crs {src.crs} | {src.width:,} x {src.height:,}")
    if str(src.crs).upper() not in ("EPSG:5070", "ESRI:102039"):
        grid = grid.to_crs(src.crs)
    for i, (zid, geom) in enumerate(zip(grid.zone_id, grid.geometry), 1):
        try:
            w = from_bounds(*geom.bounds, transform=src.transform)
            a = src.read(1, window=w)
        except Exception:
            rows.append({"zone_id": zid}); continue
        # drop nodata and the classes that aren't real cover
        a = a[(a > 0) & (a != src.nodata if src.nodata is not None else a > 0)]
        if a.size == 0:
            rows.append({"zone_id": zid}); continue
        vals, cnt = np.unique(a, return_counts=True)
        tot = cnt.sum()
        k = int(vals[cnt.argmax()])
        # the modal class covers a median of only ~56% of a cell, so the mode alone throws away
        # most of the information. keep the full area fraction per class.
        frac = {NLCD.get(int(v), f"class_{int(v)}"): float(c / tot) for v, c in zip(vals, cnt)}
        rows.append({"zone_id": zid, "land_cover_mode": k, "land_cover_label": NLCD.get(k),
                     "land_cover_frac": float(cnt.max() / tot), **frac})
        if i % 2000 == 0:
            print(f"  {i:,}/{len(grid):,}")

lc = pd.DataFrame(rows)
for c in NLCD.values():
    if c in lc.columns:
        lc[c] = lc[c].fillna(0.0)
lc.to_csv(out_dir / "flyway8_land_cover_lookup.csv", index=False)
print(f"\nwrote flyway8_land_cover_lookup.csv | {lc.land_cover_label.notna().sum():,} cells labelled")
print(lc.land_cover_label.value_counts().head(10))
