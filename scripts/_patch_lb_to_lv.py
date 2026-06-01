"""Patch cell 12 in pattern notebook: change LB → LV and update labels."""
import json, pathlib

nb_path = pathlib.Path(__file__).resolve().parent.parent / "notebooks" / "2. roughness_scale-pattern.ipynb"

with open(nb_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

# Find the cell by its unique comment line
TARGET = "# LB vs effect_ratio \u2014 one panel per fV value, 2 rows \u00d7 3 cols"

found = False
for cell in nb["cells"]:
    src = "".join(cell["source"])
    if TARGET in src:
        new_src = src
        # Comment
        new_src = new_src.replace(
            "# LB vs effect_ratio",
            "# LV vs effect_ratio",
        )
        # horizontal offset column
        new_src = new_src.replace(
            "_df['LB'] = _df['LB'].astype(float)",
            "_df['LV'] = _df['LV'].astype(float)",
        )
        new_src = new_src.replace(
            "_df.loc[_df['_aniso_label'] == _lbl, 'LB'] += _off",
            "_df.loc[_df['_aniso_label'] == _lbl, 'LV'] += _off",
        )
        # scatterplot x
        new_src = new_src.replace(
            "x='LB', y='effect_ratio'",
            "x='LV', y='effect_ratio'",
        )
        # xlabel
        new_src = new_src.replace(
            "renameit('LB', rename)",
            "renameit('LV', rename)",
        )
        # suptitle
        new_src = new_src.replace(
            r"$n_e/\\langle n \\rangle$ vs $L_B$",
            r"$n_e/\\langle n \\rangle$ vs $L_V$",
        )
        cell["source"] = new_src.splitlines(keepends=True)
        found = True
        break

if found:
    with open(nb_path, "w", encoding="utf-8") as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)
    print("OK – patched LB → LV in cell 12")
else:
    print("WARNING – target cell not found")
