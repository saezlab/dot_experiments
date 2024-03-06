import pandas as pd
import scanpy as sc
import numpy as np
from scipy import sparse

def load_libd_sample(sample_name):
    libd_dir = "../LIBD"
    
    layers = pd.read_csv(f"{libd_dir}/barcode_level_layer_map.tsv", sep="\t", header=None, dtype='str')
    layers.columns = ["spot_name", "sample_name", "layer"]
    
    layers = layers[layers["sample_name"] == sample_name].copy()
    layers.set_index("spot_name", inplace=True)
    
    loc = pd.read_csv(f"{libd_dir}/{sample_name}/tissue_positions_list.txt", sep=",", header=None, index_col=0)
    loc = loc[loc[1] == 1].iloc[:, 1:3].copy()
    loc.columns = ["Row", "Col"]
    
    spots = list(set(loc.index) & set(layers.index))
    
    loc = loc.loc[spots].copy()
    loc["layer"] = layers.loc[spots, "layer"]
    loc["sample"] = sample_name
    loc.index.name = "spot"
    
    adata = sc.read_10x_h5(f"{libd_dir}/{sample_name}/filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    adata = adata[spots, :].copy()
    
    rs = np.asarray(adata.X.sum(axis=0)).flatten()
    non_zero_indices = np.where(rs > 0)[0]
    
    adata = adata[:, non_zero_indices].copy()
    adata.obs = loc
    adata.obs_names = [f"{sample_name}_{spot}" for spot in adata.obs_names]

    return adata

def load_libd_batch(sample_names):
    return sc.concat([load_libd_sample(s) for s in sample_names], merge="same")

def load_libd_experiment(ref, srt):
    libd_dir = "../LIBD"
    
    mtdata = pd.read_csv(f"{libd_dir}/libd_clinicaldata.csv", dtype={"SlideID": str})
    mtdata.set_index("SlideID", inplace=True)
    
    sample_names = mtdata.index
    
    srt_samples = []
    if srt in sample_names:
        srt_samples.append(srt)
    else:
        try:
            srt_samples.append(sample_names[int(srt)])
        except (IndexError, ValueError):
            pass
    
    ref_samples = []
    if ref == 'aggregated':
        ref_samples = list(sample_names)
    elif ref == 'brain':
        brain_number = mtdata.loc[srt_samples[0], "Brain.Number"]
        ref_samples = sample_names[mtdata["Brain.Number"] == brain_number].tolist()
    elif ref == 'adjacent':
        brain_number = mtdata.loc[srt_samples[0], "Brain.Number"]
        position = mtdata.loc[srt_samples[0], "Position"]
        ref_samples = sample_names[(mtdata["Brain.Number"] == brain_number) & (mtdata["Position"] == position)].tolist()
    else:
        raise ValueError("Invalid ref.")
    
    ref_samples = set(ref_samples) - set(srt_samples)
    
    return load_libd_batch(ref_samples), load_libd_batch(srt_samples)

def sum_anndata_by(adata, col):
    adata.strings_to_categoricals()
    assert pd.api.types.is_categorical_dtype(adata.obs[col])

    cat = adata.obs[col].values
    indicator = sparse.coo_matrix(
        (
            np.broadcast_to(True, adata.n_obs),
            (cat.codes, np.arange(adata.n_obs))
        ),
        shape=(len(cat.categories), adata.n_obs),
    )

    return sc.AnnData(
        indicator @ adata.X,
        var=adata.var,
        obs=pd.DataFrame(index=cat.categories)
    )

def read_mop_sc():
    sc_dir = "../MOp data/scRNA_10X_v2_A"
    sc_data = sc.read_mtx(f"{sc_dir}/processed_data.mtx")
    sc_data = sc_data.T
    genes = pd.read_csv(f"{sc_dir}/processed_data_genes.csv")
    meta = pd.read_csv(f"{sc_dir}/processed_data_meta.csv", index_col = 0)
    sc_data.var_names = genes['x']
    sc_data.obs_names = meta.index
    sc_data.obs = meta

    return sc_data

def produce_mop_multicell(sample_id = "mouse1_sample1", sc_data = None,
                          seed=1, diameter=100, ann_column="cluster_label"):
    mop_dir = "../MOp data"
    multicell_meta = pd.read_csv(f"{mop_dir}/multicell_meta.csv", dtype={"X": str}, index_col=0)
    
    sc_col = f"sc_seed_{seed}"
    if sc_col not in multicell_meta.columns:
        raise ValueError("Seed not found!")
    
    d_col = f"spot_{diameter}"
    if d_col not in multicell_meta.columns:
        raise ValueError("Diameter not found!")

    if sc_data is None:
        sc_data = read_mop_sc()
    
    mf_cells = multicell_meta.index[multicell_meta['sample_id'] == sample_id].values
    spots = multicell_meta.loc[mf_cells, d_col].values
    sc_cells = multicell_meta.loc[mf_cells, sc_col].values
    
    srt_adata = sc_data[sc_cells, :].copy()
    srt_adata.obs['spot'] = spots
    srt_adata = sum_anndata_by(srt_adata, 'spot')
    
    srt_composition = pd.crosstab(index=spots, columns=sc_data.obs.loc[sc_cells, ann_column].values)
    srt_meta = pd.DataFrame(np.array(list(srt_composition.index.str.split("X", expand=True))))
    srt_meta.columns = ['x', 'y']
    srt_ids = pd.crosstab(index = spots, columns=multicell_meta.loc[mf_cells, 'slice_id'])
    srt_meta['slice_id'] = srt_ids.idxmax(axis=1)

    srt_meta.index = srt_adata.obs_names = srt_composition.index
    srt_adata.obs = srt_meta
    srt_adata.obsm['composition'] = srt_composition
    
    remaining_sc = list(set(sc_data.obs_names) -set(sc_cells))

    return sc_data[remaining_sc].copy(), srt_adata