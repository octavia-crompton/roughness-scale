"""Patch decomp notebook: remove y-tick labels from center & right columns."""
import json

NB = "/Users/octaviacrompton/Projects/roughness-scale/notebooks/3. roughness_scale-decomp.ipynb"
with open(NB) as f:
    nb = json.load(f)

for ci, cell in enumerate(nb['cells']):
    src = ''.join(cell['source'])
    if 'axes[0, col].sharey(axes[0, 0])' in src and 'Sf decomposition' in src:
        print(f"Found target cell at index {ci}")
        lines = cell['source']
        new_lines = []
        for line in lines:
            new_lines.append(line)
            if 'axes[0, col].sharey(axes[0, 0])' in line:
                new_lines.append('        axes[0, col].tick_params(labelleft=False)\n')
            elif 'axes[row, col].sharey(axes[1, 0])' in line:
                new_lines.append('        if col > 0:\n')
                new_lines.append('            axes[row, col].tick_params(labelleft=False)\n')
        cell['source'] = new_lines
        print("Patched cell")
        break
else:
    print("Cell not found!")

with open(NB, 'w') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write('\n')
print("Done")
