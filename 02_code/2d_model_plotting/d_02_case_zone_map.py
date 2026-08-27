"""Figure 2 - HPAI case zones and feedlot negative-control zones, Minnesota.

Same map as the submitted manuscript (originally d_02_case_zone_map.ipynb), rebuilt on the
8-state 10 km grid with the current case definition. Two changes from the notebook: the plain
"Study Cell" outline layer is dropped, so only case cells and feedlot cells are drawn, and
there is no title (captions live in the manuscript text).
"""

from pathlib import Path

import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib_map_utils import scale_bar

base_folder = Path(__file__).resolve().parents[2]
objects_folder = base_folder / "03_output" / "3b_model_output" / "dataframes"
figures_folder = base_folder / "05_figures" / "main"
env_data = base_folder / "01_data" / "1a_exposure_data" / "meteorological_data"

# --- Read the grid. The geojson says WGS84 but the coordinates are the projected EPSG:5070
#     grid, so set the CRS instead of reprojecting; converting it as lon/lat gives empties.
grid = gpd.read_file(env_data / "hpai_multistate_era5" / "fishnet_8state.geojson")
grid = grid.set_crs(5070, allow_override=True)
grid = grid[grid["state"] == "Minnesota"]

state = gpd.read_file(base_folder / "01_data" / "1c_supportive_datasets" / "spatial_data" / "cb_2022_us_state_5m")
state = state[state["NAME"] == "Minnesota"].to_crs(epsg=3857)

layers = pd.read_csv(objects_folder / "map_layers_mn8.csv", dtype={"zone_id": str})
case_ids = set(layers.loc[layers.layer == "case", "zone_id"])
neg_ids = set(layers.loc[layers.layer == "negative_control", "zone_id"])

# --- Project to Web Mercator for the basemap ---
case_zones = grid[grid.zone_id.isin(case_ids)].to_crs(epsg=3857)
negative_control_zones = grid[grid.zone_id.isin(neg_ids)].to_crs(epsg=3857)

# --- Plot ---
fig, ax = plt.subplots(figsize=(10, 10))

negative_control_zones.plot(ax=ax, facecolor="#1f78b4", edgecolor="#1270ae",
                            linewidth=0.4, alpha=0.5)
case_zones.plot(ax=ax, facecolor="#d73027", edgecolor="red", linewidth=0.4, alpha=0.9)
state.boundary.plot(ax=ax, color="black", linewidth=1.6)

# CARTO started stamping "API KEY REQUIRED" across its free Positron tiles some time after
# the original submission. Esri's grey canvas is free and looks near-identical.
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldGrayCanvas, crs=case_zones.crs)

scale_bar(ax=ax, location="lower right", style="boxes",
          bar={"projection": 3857, "length": 1.5, "unit": "km"},
          labels={"style": "major", "loc": "below"})

legend_elements = [
    Line2D([0], [0], color="red", marker="s", markersize=12, linestyle="None",
           markerfacecolor="#d73027", label=f"Positive Case (N = {len(case_zones):,})"),
    Line2D([0], [0], color="#1270ae", marker="s", markersize=12, linestyle="None",
           markerfacecolor="#1f78b4", alpha=0.5,
           label=f"Negative Control (N = {len(negative_control_zones):,})"),
]
ax.legend(handles=legend_elements, loc="upper right", frameon=True,
          bbox_to_anchor=(0.98, 0.98), fontsize=12)

ax.axis("off")
plt.tight_layout()
out = figures_folder / "figure2_case_zone_map.png"
plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print(f"case cells {len(case_zones)} | negative controls {len(negative_control_zones)}")
print(f"wrote {out}")
