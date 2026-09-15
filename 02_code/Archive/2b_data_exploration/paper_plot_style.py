# paper_plot_style.py - shared style, per the paper-figure skill's Step 2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FONT_SIZE = 10
DPI = 300
FORMAT = "png"          # skill defaults to pdf; PNG requested here

matplotlib.rcParams.update({
    'font.size': FONT_SIZE,
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
    'axes.labelsize': FONT_SIZE,
    'axes.titlesize': FONT_SIZE + 1,
    'xtick.labelsize': FONT_SIZE - 1,
    'ytick.labelsize': FONT_SIZE - 1,
    'legend.fontsize': FONT_SIZE - 1,
    'figure.dpi': DPI, 'savefig.dpi': DPI,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.05,
    'axes.grid': False,
    'axes.spines.top': False, 'axes.spines.right': False,
    'text.usetex': False, 'mathtext.fontset': 'stix',
})

# skill offers tab10 / Set2 / colorblind; taking the deuteranopia-safe option
COLORS = ["#0072B2", "#D55E00", "#009E73"]
MARKERS = ["o", "^", "s"]

def save_fig(fig, name, outdir=".", fmt=FORMAT):
    p = f"{outdir}/{name}.{fmt}"
    fig.savefig(p)
    print(f"Saved: {p}")
    return p
