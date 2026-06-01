import json
nb = json.load(open("notebooks/2. roughness_scale-pattern.ipynb"))
for i, c in enumerate(nb["cells"]):
    cid = c.get("id", "NO-ID")
    src = "".join(c.get("source", []))[:80].replace("\n", " ")
    print(f"{i:3d}  id={repr(cid):30s}  {src}")
