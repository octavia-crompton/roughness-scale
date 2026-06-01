#!/usr/bin/env python3
"""Fix OLS dtype error in aniso sensitivity cell."""
import json, pathlib

NB = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "3. roughness_scale-decomp.ipynb"
nb = json.loads(NB.read_text("utf-8"))

TARGET = "Sensitivity of n_e/<n> to anisotropy"
cell = None
for c in nb["cells"]:
    src = "".join(c["source"])
    if TARGET in src:
        cell = c
        break
if cell is None:
    raise SystemExit(f"Could not find cell containing '{TARGET}'")

src = "".join(cell["source"])

src = src.replace(
    "_X_aniso = sm.add_constant(\n"
    "    pd.concat([_sub_sens[['aniso']], _fv_dum, _sig_dum], axis=1))",
    "_X_aniso = sm.add_constant(\n"
    "    pd.concat([_sub_sens[['aniso']], _fv_dum, _sig_dum], axis=1).astype(float))",
)

cell["source"] = src.splitlines(keepends=True)
NB.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n", "utf-8")
print("OK – added .astype(float) to fix OLS dtype")
