"""Patch cell 11 in pattern notebook: change LV → sigma and update labels."""
import json, pathlib

nb_path = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "2. roughness_scale-pattern.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

# Target the FIRST cell that starts with this comment (cell 11, not cell 12 which is the duplicate)
TARGET = "# LV vs effect_ratio \u2014 one panel per fV value, 2 rows \u00d7 3 cols"

found = False
for cell in nb["cells"]:
    src = "".join(cell["source"])
    if TARGET in src:
        new_src = src
        # Comment
        new_src = new_src.replace(
            "# LV vs effect_ratio",
            "# sigma vs effect_ratio",
            1,  # only first occurrence
        )
        # horizontal offset column
        new_src = new_src.replace(
            "_df['LV'] = _df['LV'].astype(float)",
            "_df['sigma'] = _df['sigma'].astype(float)",
        )
        new_src = new_src.replace(
            "_df.loc[_df['_aniso_label'] == _lbl, 'LV'] += _off",
            "_df.loc[_df['_aniso_label'] == _lbl, 'sigma'] += _off",
        )
        # scatterplot x
        new_src = new_src.replace(
            "x='LV', y='effect_ratio'",
            "x='sigma', y='effect_ratio'",
        )
        # xlabel
        new_src = new_src.replace(
            "renameit('LV', rename)",
            "renameit('sigma', rename)",
        )
        # suptitle
        new_src = new_src.replace(
            r"$n_e/\\langle n \\rangle$ vs $L_V$",
            r"$n_e/\\langle n \\rangle$ vs $\\sigma$",
        )
        cell["source"] = new_src.splitlines(keepends=True)
        found = True
        break

if found:
    with open(nb_path, "w", encoding="utf-8") as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    print("OK – patched LV → sigma in cell 11")
else:
    print("WARNING – target cell not found")
