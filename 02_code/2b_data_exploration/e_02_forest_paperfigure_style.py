# What the `paper-figure` skill's house style looks like on our forest plot, next to the
# ggplot version in d_23. Kept as a comparison, not as a build step - nothing in the manuscript
# pipeline reads this.
#
# Verdict after rendering both: the one idea worth taking is the per-panel n annotation, which
# the skill's checklist is right to require and d_23 does not have yet. Times at 9-10pt reads
# worse than Avenir at print size, and dropping the panel borders makes the four lag windows
# run together. Also note the log-axis trap below - it cost a render.
#
# Data comes from d_23_figure2_primary_forest.R; run that first.
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
# Figure 2 in the paper-figure skill's house style: serif, 300 dpi, no grid,
# no top/right spines, colourblind-safe palette, no titles inside the figure.
try:
    import pandas as pd, numpy as np
except ModuleNotFoundError as e:
    raise SystemExit(
        f"\n{e}. You're on the wrong interpreter: {sys.executable}\n"
        "pandas/matplotlib only live in the hpai_outbreaks_mn conda env.\n"
        "Run it with:\n"
        "  /opt/homebrew/Caskroom/miniforge/base/envs/hpai_outbreaks_mn/bin/python "
        + __file__ + "\n"
        "or in VS Code: Cmd+Shift+P -> Python: Select Interpreter -> "
        "hpai_outbreaks_mn\n")
from matplotlib.ticker import FixedLocator, NullLocator, NullFormatter
from paper_plot_style import *

d = pd.read_csv(ROOT / "03_output/3b_model_output/dataframes/figure2_primary_forest.csv")
MODELS = ["a. Absolute conditions", "b. Anomalies (departure from the 90-day baseline)"]
PANELS = ["Minnesota", "Other flyway states", "Pooled"]
EXPO   = ["Temperature", "Precipitation", "Soil moisture", "Wind speed", "Anseriformes"]
WINS   = list(dict.fromkeys(d["window"]))
OFF    = {p: o for p, o in zip(PANELS, [-0.26, 0.0, 0.26])}

fig, axes = plt.subplots(2, len(WINS), figsize=(11, 5.8), sharex=True, sharey=True)
for r, mod in enumerate(MODELS):
    for c, win in enumerate(WINS):
        ax = axes[r, c]
        ax.axvline(1, color="0.45", lw=0.6, zorder=0)
        sub = d[(d["model"] == mod) & (d["window"] == win)]
        for k, pan in enumerate(PANELS):
            s = sub[sub["panel"] == pan]
            y = [EXPO.index(v) + OFF[pan] for v in s["variable"]]
            ax.hlines(y, s["low"], s["high"], color=COLORS[k], lw=1.1, zorder=2)
            ax.plot(s["RR"], y, MARKERS[k], color=COLORS[k], ms=4.2, lw=0,
                    label=pan if (r == 0 and c == 0) else None, zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(0.45, 4.3)
        # a log axis keeps its own minor locator/formatter and will overwrite set_xticklabels,
        # which is how the first render came out reading "0.5 x 10^-1"
        ax.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 3]))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xticklabels(["0.5", "1", "2", "3"])
        ax.set_ylim(-0.7, len(EXPO) - 0.3)
        ax.set_yticks(range(len(EXPO)))
        ax.set_yticklabels(EXPO)
        ax.invert_yaxis()
        ax.tick_params(length=3)
        if r == 0:
            ax.set_title(win, fontsize=FONT_SIZE, pad=6)
    # panel label sits outside the axes, so it is a label not an in-figure title
    n_by = d[d["model"] == mod].groupby("panel")["n"].first()
    axes[r, 0].text(0.0, 1.22 if r == 0 else 1.06, mod, transform=axes[r, 0].transAxes,
                    ha="left", va="bottom", fontsize=FONT_SIZE + 1, fontweight="bold")
    # skill checklist: "state n and the unit of replication"
    axes[r, -1].text(1.0, 1.22 if r == 0 else 1.06,
                     "n = " + ", ".join(f"{p} {n_by[p]}" for p in PANELS),
                     transform=axes[r, -1].transAxes, ha="right", va="bottom",
                     fontsize=FONT_SIZE - 2, color="0.35")

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.125, 1.02),
           ncol=3, frameon=False, handletextpad=0.4, columnspacing=1.6)
fig.supxlabel("Rate ratio per +0.5 standard deviation increase in exposure",
              fontsize=FONT_SIZE, y=0.02)
fig.subplots_adjust(hspace=0.30, wspace=0.12)
save_fig(fig, "figure2_forest_paperfigure_style", outdir=str(ROOT / "05_figures/supplement"))
