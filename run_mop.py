import sys
from process_data import produce_mop_multicell

args = sys.argv
sample_id = args[1] # e.g., "mouse1_sample1"
method = args[2] # e.g., TACCO, cell2location, ...

ref_adata, srt_adata = produce_mop_multicell(sample_id)
ref_annotation = "cluster_label"

if method.lower() == "tacco":
    from methods import run_tacco
    weights = run_tacco(ref_adata, srt_adata, ref_annotation)
elif method.lower() in ["cell2location", "c2l"]:
    from methods import run_cell2location
    weights = run_cell2location(ref_adata, srt_adata, ref_annotation, epochs = 20000)
elif method.lower() == "tangram":
    from methods import run_tangram
    weights = run_tangram(ref_adata, srt_adata, ref_annotation)
elif method.lower() == "novosparc":
    from methods import run_novosparc
    coords = srt_adata.obs[['x', 'y']]
    weights = run_novosparc(ref_adata, srt_adata , ref_annotation, coords, alpha = 0.5, ncells = 500)