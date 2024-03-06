# from process_data import *
import scanpy as sc
import scipy
import pandas as pd
import numpy as np
import tempfile
import os

def polish_weights(weights, columns = None, rows = None):
    if rows is None:
        rows = weights.index
    else:
        missing_rows = list(set(rows) - set(weights.index))
        if len(missing_rows) > 0:
            mr = pd.DataFrame(1/weights.shape[1], index=missing_rows, weights=df.columns)
            weights = weights.append(mr, ignore_index=False)
    if columns is None:
        columns = weights.columns
    else:
        missing_cols = list(set(columns) - set(weights.columns))
        if len(missing_cols) > 0:
            weights[missing_cols] = 0
            
    weights = weights.loc[rows, columns].copy()
    weights = weights.div(weights.sum(axis = 1), axis = 0)
    return weights

def sample_by(adata, ann, ncells = None):
    if ncells is None:
        return adata
    if ncells is not None:
        df = pd.DataFrame({"ann": adata.obs[ann], "cell": adata.obs_names})
        selected_cells = df.groupby('ann')['cell'].apply(lambda s: s.sample(min(len(s), ncells)))
        return adata[selected_cells, :].copy()

def run_tacco(ref_adata, srt_adata, ref_annotation):
    import tacco as tc
    tc.tl.annotate(srt_adata, ref_adata, annotation_key=ref_annotation, result_key='tacco_ct_pred')
    
    weights = srt_adata.obsm['tacco_ct_pred']
    return polish_weights(weights, sorted(set(ref_adata.obs[ref_annotation])), srt_adata.obs_names )

def run_tangram(ref_adata, srt_adata, ref_annotation, ncells = None):
    import tangram as tg

    ref_adata = sample_by(ref_adata, ref_annotation, ncells)
    
    tg.pp_adatas(ref_adata, srt_adata, genes=None)
    ad_map = tg.map_cells_to_space(ref_adata, srt_adata)
    tg.project_cell_annotations(ad_map, srt_adata, annotation=ref_annotation)

    weights = srt_adata.obsm['tangram_ct_pred']
    return polish_weights(weights, sorted(set(ref_adata.obs[ref_annotation])), srt_adata.obs_names )

def run_novosparc(ref_adata, srt_adata, ref_annotation, coords, alpha = 0.5, ncells = None):
    import novosparc as nv

    ref_adata = sample_by(ref_adata, ref_annotation, ncells)
        
    common_genes = list(set(ref_adata.var_names) & set(srt_adata.var_names))
    ref_adata = ref_adata[:, common_genes]
    srt_adata = srt_adata[:, common_genes]

    srt_X = srt_adata.X
    if scipy.sparse.issparse(ref_adata.X):
        ref_adata = sc.AnnData(ref_adata.X.A, obs=ref_adata.obs, var=ref_adata.var)
    if scipy.sparse.issparse(srt_X):
        srt_X = srt_X.A
    
    nv_tissue = nv.cm.Tissue(ref_adata, coords, srt_X)
    
    if alpha == 1.0:
        nv_tissue.setup_linear_cost(markers_metric = 'minkowski')
    else:
        nv_tissue.setup_linear_cost.__func__.__defaults__ = (None, None, 'minkowski', 2)
        nv_tissue.setup_reconstruction()
    
    ref_cell_type = pd.get_dummies(ref_adata.obs[ref_annotation])

    epsilon = 5e-4
    while True:
        nv_tissue.reconstruct(alpha, epsilon = epsilon)
        weights = (nv_tissue.gw.T @ ref_cell_type)
        if weights.sum(axis = 1).min() > 0:
            break
        epsilon *= 10

    weights.index = srt_adata.obs_names
    return polish_weights(weights, sorted(set(ref_adata.obs[ref_annotation])), srt_adata.obs_names )

"""
Add the following function to "\lib\site-packages\scvi\model\_utils.py" for scvi-tools version > 1.04 

def parse_use_gpu_arg(
    use_gpu: Optional[Union[str, int, bool]] = None,
    return_device=True,
):
    # Support Apple silicon
    cuda_available = torch.cuda.is_available()
    # If using an older version of torch.
    try:
        mps_available = torch.backends.mps.is_available()
    except AttributeError:
        mps_available = False
    gpu_available = cuda_available
    lightning_devices = None
    if (use_gpu is None and not gpu_available) or (use_gpu is False):
        accelerator = "cpu"
        device = torch.device("cpu")
        lightning_devices = "auto"
    elif (use_gpu is None and gpu_available) or (use_gpu is True):
        current = torch.cuda.current_device() if cuda_available else "mps"
        if current != "mps":
            lightning_devices = [current]
            accelerator = "gpu"
        else:
            accelerator = "mps"
            lightning_devices = 1
        device = torch.device(current)
    # Also captures bool case
    elif isinstance(use_gpu, int):
        device = torch.device(use_gpu) if not mps_available else torch.device("mps")
        accelerator = "gpu" if not mps_available else "mps"
        lightning_devices = [use_gpu] if not mps_available else 1
    elif isinstance(use_gpu, str):
        device = torch.device(use_gpu)
        accelerator = "gpu"
        # changes "cuda:0" to "0,"
        lightning_devices = [int(use_gpu.split(":")[-1])]
    else:
        raise ValueError("use_gpu argument not understood.")

    if return_device:
        return accelerator, lightning_devices, device
    else:
        return accelerator, lightning_devices
"""

def run_cell2location(ref_adata, srt_adata, ref_annotation, epochs = 20000, N_cells_per_location=30, detection_alpha=20, ncells = None):
    os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,force_device=True'
    import cell2location as c2l

    ref_adata = sample_by(ref_adata, ref_annotation, ncells)
    
    with tempfile.TemporaryDirectory(prefix='temp_c2l_', dir='.') as results_folder:
        results_folder = results_folder + '/'
        ref_run_name = f'{results_folder}/reference_signatures'
        run_name = f'{results_folder}/cell2location_map'
    
        srt_adata.var['SYMBOL'] = srt_adata.var_names
        srt_adata.var['MT_gene'] = [gene.startswith('MT-') for gene in srt_adata.var['SYMBOL']]
        srt_adata.var.set_index('SYMBOL', drop=True, inplace=True)
        srt_adata.obs['Sample'] = 'sample'
        
        # remove MT genes for spatial mapping (keeping their counts in the object)
        srt_adata.obsm['MT'] = srt_adata[:, srt_adata.var['MT_gene'].values].X.toarray()
        srt_adata = srt_adata[:, ~srt_adata.var['MT_gene'].values]
    
        ref_adata.var['SYMBOL'] = ref_adata.var_names
        # rename 'GeneID-2' as necessary for your data
        ref_adata.var.set_index('SYMBOL', drop=True, inplace=True)
        ref_adata.obs['Sample'] = 'sample'
        ref_adata.obs['Method'] = 'method'
        
        selected = c2l.utils.filtering.filter_genes(ref_adata, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
    
        # filter the object
        ref_adata = ref_adata[:, selected].copy()
    
        c2l.models.RegressionModel.setup_anndata(adata=ref_adata, 
                        # 10X reaction / sample / batch
                        batch_key='Sample', 
                        # cell type, covariate used for constructing signatures
                        labels_key=ref_annotation, 
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['Method'])
    
        mod = c2l.models.RegressionModel(ref_adata) 
        mod.train(max_epochs=epochs)
    
        ref_adata = mod.export_posterior(
            ref_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
        
        # Save model
        mod.save(f"{ref_run_name}", overwrite=True)
    
        if 'means_per_cluster_mu_fg' in ref_adata.varm.keys():
            inf_aver = ref_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                            for i in ref_adata.uns['mod']['factor_names']]].copy()
        else:
            inf_aver = ref_adata.var[[f'means_per_cluster_mu_fg_{i}' 
                                            for i in ref_adata.uns['mod']['factor_names']]].copy()
    
        
        common_genes = list(set(inf_aver.index) & set(srt_adata.var_names))
        inf_aver = inf_aver.loc[common_genes, :].copy()
        srt_adata = srt_adata[:, common_genes].copy()
    
        c2l.models.Cell2location.setup_anndata(adata=srt_adata, batch_key='Sample')
        
        mod = c2l.models.Cell2location(
            srt_adata, cell_state_df=inf_aver, 
            # the expected average cell abundance: nv_tissue-dependent 
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location=N_cells_per_location,
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection:
            detection_alpha=detection_alpha
        ) 
    
        mod.train(max_epochs=epochs, 
            # train using full data (batch_size=None)
            batch_size=None, 
            # use all data points in training because 
            # we need to estimate cell abundance at all locations
            train_size=1)
        
        srt_adata = mod.export_posterior(
            srt_adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs})
       
        weights = srt_adata.obsm['q05_cell_abundance_w_sf']
        ct_key = 'q05cell_abundance_w_sf_means_per_cluster_mu_fg_'
        weights.columns = [c[len(ct_key):] for c in weights.columns]
        return polish_weights(weights, sorted(set(ref_adata.obs[ref_annotation])), srt_adata.obs_names )

        
