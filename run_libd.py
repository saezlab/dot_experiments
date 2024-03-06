import sys
from process_data import load_libd_experiment

args = sys.argv
ref_type = args[1] # e.g., "brain"
sample_index = args[2] # e.g., "1"
method = args[3] # e.g., TACCO, cell2location, ...

ref_adata, srt_adata = load_libd_experiment(ref_type, sample_index)
ref_annotation = "layer"

if method.lower() == "tacco":
    from methods import run_tacco
    weights = run_tacco(ref_adata, srt_adata, ref_annotation)
elif method.lower() in ["cell2location", "c2l"]:
    from methods import run_cell2location
    weights = run_cell2location(ref_adata, srt_adata, ref_annotation, epochs = 2000)
elif method.lower() == "tangram":
    from methods import run_tangram
    weights = run_tangram(ref_adata, srt_adata, ref_annotation)
elif method.lower() == "novosparc":
    from methods import run_novosparc
    coords = srt_adata.obs[['Row', 'Col']]
    weights = run_novosparc(ref_adata, srt_adata , ref_annotation, coords, alpha = 0.5, ncells = 500)