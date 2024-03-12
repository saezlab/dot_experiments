from process_data import load_sc_wu, load_vis_tnbc
import sys

args = sys.argv
vis_id = args[1] # e.g., "1142243F"
method = args[2] # e.g., TACCO, cell2location, ...

ref_adata = load_sc_wu("TNBC")
ref_annotation = "celltype_major"
srt_adata = load_vis_tnbc(vis_id)

if method.lower() == "tacco":
    from methods import run_tacco
    weights = run_tacco(ref_adata, srt_adata, ref_annotation)
elif method.lower() in ["cell2location", "c2l"]:
    from methods import run_cell2location
    weights = run_cell2location(ref_adata, srt_adata, ref_annotation, 
                                epochs = 20000, normalize = False)