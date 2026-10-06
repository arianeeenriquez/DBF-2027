import aerosandbox as asb
import aerosandbox.numpy as np
import json
from aerosandbox.tools import units
import re
from pathlib import Path
import pandas as pd
import numpy as onp   # plain numpy (np above is aerosandbox.numpy, which is for symbolic math)
import casadi as ca
import time
import matplotlib.pyplot as plt
from scipy.optimize import root

COLUMNS = ["alpha", "CL", "CD", "CDp", "CM", "Top_Xtr", "Bot_Xtr", "Top_Itr", "Bot_Itr"]
NAME_RE = re.compile(r"^(?P<airfoil>.+)_Re(?P<Re>\d+(?:\.\d+)?)$")  # greedy, so underscores in the airfoil name are OK

def parse_polar_rows(path):
    """Read the numeric rows (after the dashed line) of one XFOIL polar file."""
    lines = Path(path).read_text(errors="ignore").splitlines()
    start = next((i for i, l in enumerate(lines) if l.strip().startswith("------")), None)
    rows = []
    if start is not None:
        for l in lines[start + 1:]:
            parts = l.split()
            if len(parts) == len(COLUMNS):
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    pass
    return pd.DataFrame(rows, columns=COLUMNS)

def load_airfoils(folder, pattern="*.txt"):
    """Returns {airfoil_name: DataFrame}, with all Reynolds numbers stacked in one table."""
    frames = {}
    for p in sorted(Path(folder).glob(pattern)):
        m = NAME_RE.match(p.stem)
        if not m:
            print(f"Skipping (name doesn't match pattern): {p.name}")
            continue
        df = parse_polar_rows(p)
        df.insert(0, "Re", float(m.group("Re")))
        frames.setdefault(m.group("airfoil"), []).append(df)

    return {
        name: pd.concat(dfs, ignore_index=True).sort_values(["Re", "alpha"]).reset_index(drop=True)
        for name, dfs in frames.items()
    }

class AirfoilPolar:
    """
    Smooth, differentiable CL(alpha) and CD(alpha) for ONE airfoil at ONE Reynolds number,
    built from the DataFrame that load_airfoils() returns (columns: Re, alpha, CL, CD, ...).

    - If Re falls between two Reynolds numbers in the data, the polars are blended in log(Re).
    - Only the pre-stall branch is kept (alpha up to CL_max); post-stall XFOIL data is unreliable.
    - Evaluated as a cubic B-spline (casadi), so it works inside asb.Opti and IPOPT gets real
      derivatives. (Plain linear interpolation has kinks and made IPOPT fail in testing.)
    """
    def __init__(self, df, Re, n_grid=60):
        self.Re_requested = float(Re)
        res = onp.sort(df["Re"].unique())
        if len(res) == 0:
            raise ValueError("polar file(s) contain no converged data points")
        self.Re = float(onp.clip(Re, res.min(), res.max()))
        self.Re_clamped = bool(self.Re != Re)  # True if the data didn't cover the requested Re

        lo, hi = res[res <= self.Re].max(), res[res >= self.Re].min()
        w = 0.0 if lo == hi else (onp.log(self.Re) - onp.log(lo)) / (onp.log(hi) - onp.log(lo))

        def clean(r):
            d = df[df.Re == r].sort_values("alpha").drop_duplicates("alpha")
            return d[onp.isfinite(d.CL) & onp.isfinite(d.CD) & (d.CD > 0)]
        d_lo, d_hi = clean(lo), clean(hi)
        if len(d_lo) < 8 or len(d_hi) < 8:
            raise ValueError("too few converged polar points")

        # blend the two bracketing polars on their common alpha range
        amin = max(d_lo.alpha.min(), d_hi.alpha.min())
        amax = min(d_lo.alpha.max(), d_hi.alpha.max())
        a = onp.linspace(amin, amax, 300)
        cl = (1 - w) * onp.interp(a, d_lo.alpha, d_lo.CL) + w * onp.interp(a, d_hi.alpha, d_hi.CL)
        cd = (1 - w) * onp.interp(a, d_lo.alpha, d_lo.CD) + w * onp.interp(a, d_hi.alpha, d_hi.CD)

        # keep the pre-stall branch only, then resample to a uniform grid for the spline
        i_stall = int(onp.argmax(cl))
        self.alpha_min, self.alpha_stall, self.CL_max = float(a[0]), float(a[i_stall]), float(cl[i_stall])
        if i_stall < 10:
            raise ValueError("polar has no usable pre-stall range")
        grid = onp.linspace(self.alpha_min, self.alpha_stall, n_grid)
        self.alpha_grid = grid
        self.CL_grid = onp.interp(grid, a[: i_stall + 1], cl[: i_stall + 1])
        self.CD_grid = onp.interp(grid, a[: i_stall + 1], cd[: i_stall + 1])
        self.CD_min = float(self.CD_grid.min())

        self._cl = ca.interpolant("cl", "bspline", [grid], self.CL_grid)
        self._cd = ca.interpolant("cd", "bspline", [grid], self.CD_grid)
        self._alpha = ca.interpolant("alpha_from_cl", "bspline", [self.CL_grid], self.alpha_grid)

    def cl(self, alpha):  # alpha in degrees; works with opti variables
        return self._cl(alpha)

    def cd(self, alpha):
        return self._cd(alpha)

    def alpha(self, cl):
        return self._alpha(cl)

x=load_airfoils("airfoil/polars")

y = AirfoilPolar(x["sd7032"], 500000)

cds = {}
for item in x:
    try:
        foil = AirfoilPolar(x[item], 500000)
        if asb.Airfoil(item.lower()).max_thickness() > 0.095 and foil.CL_max>=1.5:   
            a = foil.alpha(0.6)
            cds[item] = foil.cd(a)
            # print(f"{item}, {foil.cd(a)}")
        
    except Exception as e:
        # print(f"Skipping {x[item]}: {e}")
        continue    


top_10 = sorted(cds.items(), key=lambda x: x[1], reverse=False)[:10]
print(top_10)
