"""Figure 2 - HPAI spillover case cells across the northern Mississippi Flyway.

Rebuilt from the Minnesota-only version. Three changes:
  * all seven northern Mississippi Flyway states, not just Minnesota
  * feedlot negative-control cells dropped - they never enter the case-crossover risk set,
    so drawing them invited a comparison the design does not make
  * Minnesota gets a bold outline, the other six a lighter one

Flyway membership follows the USFWS administrative flyways
(https://www.fws.gov/partner/migratory-bird-program-administrative-flyways).
"""

from pathlib import Path

import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patheffects as pe
from shapely.geometry import box
from matplotlib.lines import Line2D
from matplotlib_map_utils import scale_bar

base_folder = Path(__file__).resolve().parents[2]
objects_folder = base_folder / "03_output" / "3b_model_output" / "dataframes"
figures_folder = base_folder / "05_figures" / "main"
env_data = base_folder / "01_data" / "1a_exposure_data" / "meteorological_data"
spatial_data = base_folder / "01_data" / "1c_supportive_datasets" / "spatial_data"

FLYWAY = ["Minnesota", "Wisconsin", "Iowa", "Michigan", "Illinois", "Indiana", "Ohio"]
FOCAL = "Minnesota"

# --- Grid. The geojson declares WGS84 but its coordinates are the projected EPSG:5070 grid,
#     so set the CRS rather than reproject; transforming it as lon/lat returns empties.
grid = gpd.read_file(env_data / "hpai_multistate_era5" / "fishnet_8state.geojson")
grid = grid.set_crs(5070, allow_override=True)
grid = grid[grid["state"].isin(FLYWAY)]

CRS_MAP = 5070   # Albers Equal Area Conic, per the lab convention - NOT Web Mercator

states = gpd.read_file(spatial_data / "cb_2022_us_state_5m")
states = states[states["NAME"].isin(FLYWAY)].to_crs(epsg=CRS_MAP)
focal = states[states["NAME"] == FOCAL]
others = states[states["NAME"] != FOCAL]

layers = pd.read_csv(objects_folder / "map_layers_flyway.csv", dtype={"zone_id": str})
case_ids = set(layers.loc[layers.layer == "case", "zone_id"])
case_zones = grid[grid.zone_id.isin(case_ids)].to_crs(epsg=CRS_MAP)

n_mn = int((layers.state == FOCAL).sum())
n_other = int((layers.state != FOCAL).sum())

# --- Analysis grid, clipped to the seven states. Drawing the whole 8-state fishnet is both
#     slow and muddy at this extent, so dissolve the states to one polygon and keep only the
#     cells that intersect it.
flyway_union = states.geometry.union_all()
grid_web = grid.to_crs(epsg=CRS_MAP)
grid_shown = grid_web[grid_web.intersects(flyway_union)]

# --- Plot ---
fig, ax = plt.subplots(figsize=(12, 11))

# faded analysis grid underneath everything
grid_shown.boundary.plot(ax=ax, color="#7a7a7a", linewidth=0.15, alpha=0.25, zorder=2)
# other states, then case cells, then Minnesota's heavier outline on top
others.boundary.plot(ax=ax, color="#2b2b2b", linewidth=1.2, zorder=3)

# A 10 km cell is under 1% of the width of a seven-state map, so the polygons themselves are
# invisible at this extent. Draw them as centroid markers instead - the cell footprint is not
# the point, the location is.
case_pts = case_zones.copy()
case_pts["geometry"] = case_pts.geometry.centroid
case_pts.plot(ax=ax, color="#d73027", markersize=52, marker="s",
              edgecolor="#5c0000", linewidth=0.7, alpha=0.95, zorder=4)
focal.boundary.plot(ax=ax, color="black", linewidth=3.0, zorder=5)

# Everything outside the seven states - the Great Lakes, Ontario, the neighbouring states -
# is context, not study area, and at this extent it dominates. Wash it out with a translucent
# mask so the eye lands on the flyway.
xmin0, ymin0, xmax0, ymax0 = states.total_bounds
pad = max(xmax0 - xmin0, ymax0 - ymin0) * 0.05
frame = box(xmin0 - pad, ymin0 - pad, xmax0 + pad, ymax0 + pad)
mask = gpd.GeoDataFrame(geometry=[frame.difference(flyway_union)], crs=states.crs)
mask.plot(ax=ax, facecolor="white", edgecolor="none", alpha=0.62, zorder=6)

# Tighten the extent onto the seven states before the basemap call - contextily fetches
# tiles for whatever the axis limits are at the time, so the crop has to happen first.
xmin, ymin, xmax, ymax = states.total_bounds
padx, pady = (xmax - xmin) * 0.01, (ymax - ymin) * 0.01
ax.set_xlim(xmin - padx, xmax + padx)
ax.set_ylim(ymin - pady, ymax + pady)

# Basemap LAST. add_basemap sets the axis limits to whatever tile extent it fetches, so
# calling it before any data is plotted pins the axes to the default 0-1 box and everything
# else lands off-screen - which is exactly how this map came out blank.
# CARTO started stamping "API KEY REQUIRED" across its free Positron tiles after the original
# submission; Esri's topo basemap is free and renders lakes and rivers clearly.
ctx.add_basemap(ax, source=ctx.providers.Esri.WorldTopoMap, crs=case_zones.crs,
                zoom=8, alpha=0.55, zorder=1, attribution=False)

# --- State labels at each polygon's representative point
for _, row in states.iterrows():
    pt = row.geometry.representative_point()
    ax.annotate(row["NAME"], xy=(pt.x, pt.y), ha="center", va="center",
                fontsize=13 if row["NAME"] == FOCAL else 11,
                fontweight="bold" if row["NAME"] == FOCAL else "normal",
                color="#111111", zorder=8,
                path_effects=[pe.withStroke(linewidth=3.2, foreground="white")])

# explicit length: the auto-calculated bar is sized for a city-scale map and warns at this extent
# low, wide bar: "ticks" reads as a horizontal rule rather than a stack of boxes
scale_bar(ax=ax, location="lower right", style="ticks",
          bar={"projection": CRS_MAP, "max": 300, "major_div": 3, "minor_div": 0,
               "height": 0.12, "unit": "km"},
          labels={"style": "major", "loc": "below", "fontsize": 10})

legend_elements = [
    Line2D([0], [0], color="#5c0000", marker="s", markersize=10, linestyle="None",
           markerfacecolor="#d73027",
           label=f"Spillover case cell (N = {len(case_zones):,})"),
    Line2D([0], [0], color="black", linewidth=2.6, label=f"{FOCAL} (N = {n_mn})"),
    Line2D([0], [0], color="#3a3a3a", linewidth=1.1,
           label=f"Other northern Mississippi Flyway states (N = {n_other})"),
    Line2D([0], [0], color="#7a7a7a", linewidth=0.6, alpha=0.6,
           label="10 km analysis grid"),
]
leg = ax.legend(handles=legend_elements, loc="upper right", frameon=True,
                facecolor="white", framealpha=0.94, edgecolor="#666666",
                fontsize=13, labelspacing=0.7, handlelength=2.1,
                borderpad=0.85, handletextpad=0.85)
leg.set_zorder(10)

ax.set_axis_off()
fig.tight_layout()
fig.savefig(figures_folder / "figure2_case_zone_map.png", dpi=300,
            bbox_inches="tight", facecolor="white")
print(f"wrote figure2_case_zone_map.png | {len(case_zones)} case cells "
      f"({n_mn} MN, {n_other} other flyway)")
