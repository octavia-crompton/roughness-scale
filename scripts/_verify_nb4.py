import json
nb4 = json.load(open("notebooks/4. roughness_scale-cf.ipynb"))
print(f"{len(nb4['cells'])} cells")
for i, c in enumerate(nb4["cells"]):
    src = "".join(c.get("source", []))[:65].replace("\n", " | ")
    print(f"  [{i}] {c['cell_type']:8s} {src}")
