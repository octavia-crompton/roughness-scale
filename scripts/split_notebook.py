"""Split roughness_scale-analysis.ipynb into two focused notebooks.

Output:
  notebooks/roughness_scale-pattern.ipynb  -- spatial pattern + storm char.
  notebooks/roughness_scale-decomp.ipynb   -- decomposition + prediction methods
"""
import json, copy, os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', 'notebooks', 'roughness_scale-analysis.ipynb')

with open(SRC) as f:
    nb = json.load(f)
cells = nb['cells']

id_to_idx = {c.get('id'): i for i, c in enumerate(cells)}

def get_cell(cid):
    if cid not in id_to_idx:
        raise KeyError(f"Cell id not found: {cid!r}")
    return cells[id_to_idx[cid]]

# -- Shared setup cells (indices 0-6) -----------------------------------------
SETUP_IDS = [
    'bf3f668d',  # markdown intro / out_dir config
    'f96a48ee',  # imports
    'e30afedb',  # figure registry configure
    '13164de1',  # swof modules
    '35e4e46a',  # load summary
    '29eb7207',  # USE_HYDRO toggle
    'f7d0d8f9',  # markdown "Visualise"
]

# -- Pattern-only content cells -----------------------------------------------
PATTERN_EXTRA_IDS = [
    '083d8318',  # ## Spatial pattern (markdown)
    '2782d737',  # blob images
    'e3d979d0',  # band mask setup
    '35ee167b',  # sigma grid setup
    '4496bef1',  # sigma grid plot
    '6283f78f',  # sigma grid save (fig3)
    'a156b859',  # ## Sensitivity to storm characteristics (markdown)
    '0dfc9b0b',  # storm scatter plot
]

# -- Decomp content: everything from CF section onward, excluding pattern cells
PATTERN_ONLY_SET = set(PATTERN_EXTRA_IDS)

DECOMP_START_ID = 'dbb03d73'  # ## Correction factor
decomp_start_idx = id_to_idx[DECOMP_START_ID]
decomp_extra_cells = [
    c for c in cells[decomp_start_idx:]
    if c.get('id') not in PATTERN_ONLY_SET
]

# -- Build notebooks ----------------------------------------------------------
def build_nb(selected_cells):
    new = copy.deepcopy(nb)
    new['cells'] = copy.deepcopy(selected_cells)
    return new

pattern_cells = [get_cell(cid) for cid in SETUP_IDS + PATTERN_EXTRA_IDS]
decomp_cells  = [get_cell(cid) for cid in SETUP_IDS] + decomp_extra_cells

pattern_nb = build_nb(pattern_cells)
decomp_nb  = build_nb(decomp_cells)

out_dir_nb = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'notebooks')
pattern_path = os.path.join(out_dir_nb, 'roughness_scale-pattern.ipynb')
decomp_path  = os.path.join(out_dir_nb, 'roughness_scale-decomp.ipynb')

with open(pattern_path, 'w') as f:
    json.dump(pattern_nb, f, indent=1)

with open(decomp_path, 'w') as f:
    json.dump(decomp_nb, f, indent=1)

print(f"Created {pattern_path}")
print(f"  -> {len(pattern_cells)} cells ({len(SETUP_IDS)} setup + {len(PATTERN_EXTRA_IDS)} pattern/storm)")
print(f"Created {decomp_path}")
print(f"  -> {len(decomp_cells)} cells ({len(SETUP_IDS)} setup + {len(decomp_extra_cells)} decomp)")
