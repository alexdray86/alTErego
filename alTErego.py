import numpy as np
import pandas as pd
import os, sys, glob
import subprocess as sub
import scipy.stats 
import argparse
import csv, gzip
import inspect 
import tqdm
import scipy.integrate as spi
from sh import gunzip
from multiprocessing import Pool, cpu_count
from dask.diagnostics import ProgressBar
from arboreto.utils import load_tf_names
from arboreto.algo import genie3, grnboost2, _prepare_input
from arboreto.core import SGBM_KWARGS, RF_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import to_tf_matrix, target_gene_indices, infer_partial_network

### Parse script parameters ###
parser = argparse.ArgumentParser(description='alTErego : compute TE-derived regulons from expression matrix.')
parser.add_argument('-i', '--in_matrix', type=str,
        help='Expression matrix to use. Format needs to be the following : genes [row] x observations [columns]. No column names are allowed. Row names should be gene name in Ensembl format (e.g. ENS000004757). CSV format (comma separated). Each row has therefore the following shape : ENSG00000001084,0,0,0,0,0,0,0,0,0,0,0.699620226417964,0, ... ,0 ')
parser.add_argument('-o', '--out_dir', type=str,
                    help='Output directory to write results')
parser.add_argument('--tf_db', default='db/tfs_lambert2018.csv', type=str,
        help='Text file representing a table of Transcription Factors to use in the subsequent analysis. The default given is in db/tfs_lambert2018.csv. The 2 first columns matters : First column should match row name index from expression matrix (e.g. Ensembl ref). Second column is the gene symbol for latter usage.')
parser.add_argument('--te_counts', default='db/count_mat.txt', type=str,
        help='Matrix containing TE counts in +/- 50kb around TSS of all coding genes. Gene are named according to ENSEMBL reference and should agree with input matrix gene names.')
parser.add_argument('--te_enrich', default='db/tfs_pyTEnrich_v05_summary.txt', type=str,
        help='Table with 2 columns : first column represent TF name (gene symbol, ref..) and second column TE subfamily previously shown to be enrich for peaks of the associated TF.')
parser.add_argument('--gene_ref', default='db/2005_hg19_ens_coding_genes_body_symbol_win50kb.bed', type=str,
        help='Reference table with gene traduction (by default, ENSEMBL ref -> gene name symbol)')
parser.add_argument('--n_top_gene_auc', default=2000, type=int,
        help='First N genes considered for TE enrichment in modules, to generate regulons.')
parser.add_argument('--min_regulon_size', default=5, type=int,
        help='Minimum number of genes in final regulon')
parser.add_argument('--quantile_importance', default=0.75, type=float,
        help='Quantile of genes to consider for co-expression measures. e.g. quantile_importance = 0.9 means that we consider only the top 10%% co-expression values.')
parser.add_argument('--n_top_targets', default=50, type=int,
        help='Top N targets for each TF we consider to create modules.')
parser.add_argument('--n_top_tf', default=10, type=int,
        help='Top N TFs for each target we consider to create modules.')
parser.add_argument('--min_module_size', default=20, type=int,
        help='Minimum number of genes in modules.')
parser.add_argument('--n_cpu', default=6, type=int,
        help='Number of CPU to use during multi-processing parts of the scripts')

args = parser.parse_args()

### Define constants ###
IN_MATRIX = args.in_matrix 
OUT_DIR   = args.out_dir 
TF_DB     = args.tf_db 
TE_COUNTS_FILE = args.te_counts # TE counts around TSS
TE_ENRICH_FILE = args.te_enrich # TF->TE enrichment
GENE_REF = pd.read_csv(args.gene_ref, sep = "\t", header = None)
GENE_REF.index = GENE_REF[3]
N_TOP = args.n_top_gene_auc
MIN_REG_SIZE = args.min_regulon_size
N_CPU = args.n_cpu
QUANTILE_IM = args.quantile_importance
N_TOP_TARGETS = args.n_top_targets
N_TOP_TF = args.n_top_tf
MIN_MODULE_SIZE = args.min_module_size


### Define general functions ###
# Return the name of the file without directory path nor extension
def basen_no_ext(my_str):
    return os.path.splitext(os.path.basename(my_str))[0]

# Function to print Class::Name - Comments from any classes

def logger(comment):
    """
    logger::Write class / function name with custom comment

    :param comment: custom comment as string
    """
    class_fun_str = inspect.stack()[1][0].f_locals["self"].__class__.__name__ + "::" + \
                    inspect.stack()[1].function + "()"
    print("{:<40} - {:<50}".format(class_fun_str, comment))

# Function to create directory
def create_dir(d):
    mkdir_cmd = "mkdir -p " + d
    sub.run(mkdir_cmd, check=True, shell=True)

# Detect gzip compressed file
def test_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def get_gene_symbol(g):
    if g in GENE_REF.index:
        g_symbol = GENE_REF.loc[g][6]
        if len(np.array([GENE_REF.loc[g][6], 1])) > 2:
            return g
        else :
            return g_symbol
    else:
        return g

def get_gene_symbols(g):
    gs_len = len(g)
    gsymbols = np.empty(gs_len, dtype=object)
    for i in range(gs_len):
        gsymbols[i] = get_gene_symbol(g[i])
    return gsymbols

def erase_if_exists(file):
    if os.path.exists(file):
        os.remove(file)

### CLASS DEFINITION ###

# This class handles calls to bedtools to make intersections, or use genome_subset if needed
class coExpression_launcher(object):
    def __init__(self, in_matrix, tf_db):
        """
        Allows to load all necessary object and lauch co-expression analysis using TFs expression as predictors of gene expression.

        :param in_matrix: expression matrix, can be in following format :
                                * a pandas DataFrame (rows=observations, columns=genes)
                                * a dense 2D numpy.ndarray
                                * a sparse scipy.sparse.csc_matrix
        :param tf_db: a table containing transcription factors to be considered for next steps. By default, TFs as defined in [...] are used. 
        """
        self.in_matrix = None
        self.tf_db = None
        self.adj = None
        self.suffix_exp = basen_no_ext(in_matrix)
        self.ex_matrix_filt_np = None
        self.gene_names = None
        self.tf_names = None
        self.tf_matrix = None
        self.tf_matrix_gene_names = None
        self.method_params = None
        self.target_gene_indices = None

        ## Build output directory :
        create_dir(OUT_DIR + "/results")

        ## Load files
        self.load_files(in_matrix, tf_db)

    def load_files(self, in_matrix_file, tf_db_file):
        logger("loading matrix and TF table files")
        # Load expression matrix
        ex_matrix = pd.read_csv(in_matrix_file, sep=',', header=None, index_col=0)
        self.in_matrix = ex_matrix.T
        #check_matrix(self.in_matrix)
        # Load TF tables
        self.tf_db = pd.read_csv(tf_db_file, header = None)
        #check_tf_db(self.tf_db)
    
    def _prepare_input(self, expression_data,gene_names,tf_names):
        """
        Wrangle the inputs into the correct formats.

        :param expression_data: one of:
                                * a pandas DataFrame (rows=observations, columns=genes)
                                * a dense 2D numpy.ndarray
                                * a sparse scipy.sparse.csc_matrix
        :param gene_names: optional list of gene names (strings).
                           Required when a (dense or sparse) matrix is passed as 'expression_data' instead of a DataFrame.
        :param tf_names: optional list of transcription factors. If None or 'all', the list of gene_names will be used.
        :return: a triple of:
                 1. a np.ndarray or scipy.sparse.csc_matrix
                 2. a list of gene name strings
                 3. a list of transcription factor name strings.
        """

        if isinstance(expression_data, pd.DataFrame):
            expression_matrix = expression_data.values
            gene_names = list(expression_data.columns)
        else:
            expression_matrix = expression_data
            assert expression_matrix.shape[1] == len(gene_names)

        if tf_names is None:
            tf_names = gene_names
        elif tf_names == 'all':
            tf_names = gene_names
        else:
            if len(tf_names) == 0:
                raise ValueError('Specified tf_names is empty')

            if not set(gene_names).intersection(set(tf_names)):
                raise ValueError('Intersection of gene_names and tf_names is empty.')

        return expression_matrix, gene_names, tf_names

    def prepare_input(self):
        logger("prepare data before co-expression analysis")

        # Prepare needed variables
        tf_names = list(self.tf_db[0])
        gene_names = self.in_matrix.columns
        # prepare matrix using arboreto functions
        self.ex_matrix_filt_np, self.gene_names, self.tf_names = self._prepare_input(self.in_matrix, gene_names, tf_names)
        self.tf_matrix, self.tf_matrix_gene_names = to_tf_matrix(self.ex_matrix_filt_np, self.gene_names, self.tf_names)

        # Paramters used by GRNBoost 2 :
        self.method_params = ['GBM',        # regressor_type
                              SGBM_KWARGS ] # regressor_kwargs

        # Index of target genes
        self.target_gene_indices = target_gene_indices(gene_names, target_genes='all')

    def run_infer_partial_network(self, target_gene_index):
        target_gene_name = self.gene_names[target_gene_index]
        target_gene_expression = self.ex_matrix_filt_np[:, target_gene_index]

        n = infer_partial_network(
            regressor_type = self.method_params[0],
            regressor_kwargs = self.method_params[1],
            tf_matrix = self.tf_matrix,
            tf_matrix_gene_names = self.tf_matrix_gene_names,
            target_gene_name = target_gene_name,
            target_gene_expression = target_gene_expression,
            include_meta = False,
            early_stop_window_length = EARLY_STOP_WINDOW_LENGTH,
            seed = None )
        return( n )


class Modules(object):
    def __init__(self, raw_adj, quantile_IM=0.9, n_top_targets=50, n_top_tf=5, n_min_genes=20):
        self.raw_adj  = raw_adj
        self.adj_filt = None
        self.modules  = None
        self.quantile_IM = quantile_IM
        self.threshold_IM = self.raw_adj['importance'].quantile(quantile_IM)
        self.n_top_targets = n_top_targets
        self.n_top_tf = n_top_tf
        self.n_min_genes = n_min_genes

    def get_modules(self):
        logger("Filter co-expression TF->target to form modules")
        # 1. Filter by Importance Mesure
        adj = self.raw_adj.iloc[np.array(self.raw_adj['importance'] > self.threshold_IM)]
        adj.sort_values(by=['TF','importance'], ascending=False)

        # 2. Filter by top N targets per TF
        df_top_targets = []
        for tf_name, df_grp in adj.groupby(by='TF'):
            module = df_grp.nlargest(self.n_top_targets, 'importance')
            df_top_targets.append(module)
        # Update adj data frame
        adj = pd.concat(df_top_targets)

        # 3. Filter by top N TFs per targets
        # Needs to re-sort table first
        adj = adj.sort_values(by=['target','importance'], ascending=False)
        df_top_tf = []
        for tf_name, df_grp in adj.groupby(by='target'):
            module = df_grp.nlargest(self.n_top_tf, 'importance')
            df_top_tf.append(module)
        # Update adj data frame
        adj = pd.concat(df_top_tf)

        # 4. Keep only modules with at least N genes
        adj = adj.sort_values(by=['TF','importance'], ascending=False)
        df_top_size = []
        for tf_name, df_grp in adj.groupby(by='TF'):
            if df_grp.shape[0] > self.n_min_genes:
                df_top_size.append(df_grp)
        # Make final module table
        self.modules = pd.concat(df_top_size)


def compute_areas(arr, mean_arr, n, min_g):
    areas = [0]
    for n_max in range(1,n):
        arr_max = arr[n_max]
        y_arr = np.copy(arr[0:n_max])
        if arr_max >= min_g:
            mu_frac = mean_arr[n_max-1] / y_arr[n_max-1]
            y_arr /= arr_max
            mu_h = y_arr[n_max-1] * mu_frac
            mu_area = mu_h / 2
            area = spi.trapz(y_arr, dx = 1/n_max)
            areas.append(area) # - mu_area)
        else:
            areas.append(0)
    return areas

class Te_Regulon(object):
    def __init__(self, module, tf_gref):
        self.module = module # Record full module info
        self.tf_gref = tf_gref # Ensembl ref of the module tf
        self.tf_gsymb = get_gene_symbol(tf_gref) # Gene symbol for the module tf
        self.target_gref = self.module['target'] # targets' Ensembl references
        self.target_gsymb = get_gene_symbols(np.array(self.target_gref)) # targets' gene names
        self.cumul_mat = None
        self.mean_arr = None
        self.aucs = None

    def cumul_recovery(self, te_counts):
    # Get cumulative recovery for each TE subfam w.r.t a specific gene set
        cumul_mat = np.zeros(shape=(len(te_counts.index), N_TOP))
        i=0
        for te_sfam in te_counts.index:
            # Get the N top genes associated with TEs
            counts_top = te_counts.loc[te_sfam].sort_values(ascending = False).head(N_TOP)
            top_te = counts_top.index
            # Check if top ranked genes are included in module
            arr_bool = pd.Series(top_te).isin(pd.Series(self.module['target']))
            # Transform array of boolean into cumulative distribution using numpy cumsum
            cumul_mat[i,:] = np.cumsum(arr_bool)
            i+=1
        self.cumul_mat = cumul_mat

    def mean_curves(self):
        # Get average curve over all TE subfam
        self.mean_arr = self.cumul_mat.mean(axis=0)

    def get_aucs(self, te_counts_idx):
    # Iterate over TE subfams to compute AUCs values and get max AUC
        # Define empty lists to be used
        tsfam_max_areas, tsfam_idx_max, tsfam_val_max, reg_n_genes = [], [], [], []
        for x in range(np.shape(self.cumul_mat)[0]):
            te_sfam_curve = np.copy(self.cumul_mat[x,:])
            # Compute the actual AUC value for a given te subfam
            areas = compute_areas(te_sfam_curve, self.mean_arr, N_TOP, MIN_REG_SIZE)
            # Get the maximum AUC value and index
            val_max = np.nanmax(areas)
            idx_max = np.nanargmax(areas)
            n_genes_regulon = len(np.unique(te_sfam_curve[0:idx_max+1])) - 1
            # Append results into their specific lists
            tsfam_max_areas.append(areas)
            tsfam_idx_max.append(idx_max)
            tsfam_val_max.append(val_max)
            reg_n_genes.append(n_genes_regulon)
        # Assemble results into one dataframe
        df_res = pd.DataFrame([tsfam_val_max, tsfam_idx_max, reg_n_genes],
                          columns = te_counts_idx, index = ['Max_AUC', 'idx', 'reg_n_genes']).transpose()
        # Get Normalized AUC for subset
        df_res['AUC_norm'] = ( df_res['Max_AUC'] - df_res['Max_AUC'].mean() ) / df_res['Max_AUC'].std()
        # Sort data by AUC values
        df_res = df_res.sort_values(by = 'Max_AUC', ascending = False)
        # Incorporate results in class variable
        self.aucs = df_res


class Te_regulons_maker(object):
    def __init__(self, modules, counts_file, te_enrich_file, out_dir):
        self.modules  = modules
        self.counts   = pd.read_csv(counts_file, index_col = 0)
        self.tfs      = self.modules['TF'].unique()
        self.out_dir  = out_dir
        self.te_enrich= pd.read_csv(te_enrich_file, sep = "\t", header = None, index_col = 0)
        self.regulons = {}
        self.create_dirs()

    def create_dirs(self):
        create_dir(self.out_dir)
        create_dir(self.out_dir + "/figs")

    def create_module(self, tf_ens):
    # Subset module for TF of interest / genes present in module
        single_module = self.modules.iloc[np.array(self.modules['TF'] == tf_ens)]
        single_module = single_module.iloc[np.array(single_module['target'].isin(self.counts.columns))]
        self.regulons[tf_ens] = Te_Regulon(single_module, tf_ens)

    def analyze_single_regulon(self, tf_ens):
        print(get_gene_symbol(tf_ens), end = ',')
        self.create_module(tf_ens)
        self.regulons[tf_ens].cumul_recovery(self.counts)
        self.regulons[tf_ens].mean_curves()
        self.regulons[tf_ens].get_aucs(self.counts.index)
        self.write_detailed_regulon(tf_ens)

        debug_file = self.out_dir + "/debug.txt"
        file = open(debug_file, "a")
        file.write("{}\n".format(tf_ens))
        file.close()

    def write_detailed_regulon(self, tf_ens):
        # Get TF symbol
        tf_name = get_gene_symbol(tf_ens)
        # Erase output file if exists
        out_file = '{0}/{1}_regulons_chip.csv'.format(self.out_dir, tf_name)
        erase_if_exists(out_file)
        # Write Regulons with known TF->TE relation from chipseqs
        if tf_name in self.te_enrich.index:
            str_te_sfam_enrich = self.te_enrich.loc[tf_name][1]
            sel_te_sfam_enrich = list(str_te_sfam_enrich.split(',')) 
            # Final parsing of top regulons and write results in table
            count=1
            for te_sfam in self.regulons[tf_ens].aucs.index:
                if te_sfam in sel_te_sfam_enrich:
                    idx_max = self.regulons[tf_ens].aucs.loc[te_sfam]['idx'].astype(int)
                    # Get the N top genes associated with TEs
                    counts_top = self.counts.loc[te_sfam].sort_values(ascending = False).head(idx_max)
                    top_te = counts_top.index
                    # Check if top ranked genes are included in module
                    mod_sub = self.regulons[tf_ens].module
                    arr_bool = pd.Series(top_te).isin(pd.Series(mod_sub['target']))
                    te_regulon = mod_sub.iloc[np.array(pd.Series(mod_sub['target']).isin(top_te[arr_bool])),:]
                    te_regulon['AUC_norm'] = np.array(np.repeat(self.regulons[tf_ens].aucs.loc[te_sfam]['AUC_norm'],
                                                                te_regulon.shape[0]))
                    # Get TE counts around genes
                    te_regulon.index = te_regulon['target']
                    counts_serie = counts_top[np.array(counts_top.index.isin(top_te[arr_bool]))].astype(int)
                    te_regulon = te_regulon.loc[counts_serie.index] # resort pandas
                    te_regulon['te_around_tss'] = counts_top[np.array(counts_top.index.isin(top_te[arr_bool]))].astype(int)
                    # Select only putative target with more than 1 TE
                    te_regulon = te_regulon.iloc[np.array(te_regulon['te_around_tss'] > 0)]
                    # TE subfam name involved in regulon
                    te_regulon['TE_subfam'] = np.array(np.repeat(te_sfam, te_regulon.shape[0]))
                    # We define TF symbol
                    te_regulon['TF_symbol'] = np.array(np.repeat(tf_name, te_regulon.shape[0]))
                    # Regulon target symbol
                    te_regulon['target_symbol'] = te_regulon['target']
                    te_regulon.index = te_regulon['target']
                    # We check for which target we find a symbol
                    df_targets = GENE_REF.iloc[GENE_REF.index.isin(te_regulon['target'])]
                    # Whenever available, we replace ensembl ref with symbol
                    te_regulon.loc[df_targets.index, 'target_symbol'] = df_targets[6]
                    # Write results
                    te_regulon = te_regulon[["TF_symbol", "TF", "target_symbol", "target", "TE_subfam",
                                                     "te_around_tss", "AUC_norm", "importance"]]
                    if count > 1:
                        te_regulon.to_csv(out_file, mode='a', header=False, sep = "\t", index=False)
                    else:
                        te_regulon.to_csv(out_file, mode='a', header=True, sep = "\t", index=False)
                        count+=1

                    #if self.regulons[tf_ens].aucs.loc[te_sfam]['AUC_norm'] > 1.0:
                        #self.plot_regulon_auc(tf_ens, te_sfam, te_regulon)


if __name__ == "__main__":
    
    ## 1. Analysis of co-expression with coExpression_launcher object (call arboreto GRNBoost2)
    coExpr_obj = coExpression_launcher(IN_MATRIX, TF_DB)
    coExpr_obj.prepare_input()
    # Main multiprocessing function
    print("Launching multi-processing for co-expression analysis. TFs as predictors, one process per target gene.")
    with Pool(N_CPU) as p:
        adjs = list(tqdm.tqdm(p.imap(coExpr_obj.run_infer_partial_network,
                                     list(target_gene_indices(coExpr_obj.gene_names, target_genes='all'))),
                              total=len(coExpr_obj.gene_names),
                              position = 0, leave = True))
    # Concatenate results
    adj = pd.concat(adjs).sort_values(by='importance', ascending=False)

    ## 2. Generate Modules objects by filtering co-expression values
    mod_obj = Modules(adj, quantile_IM=QUANTILE_IM, n_top_targets=N_TOP_TARGETS, n_top_tf=N_TOP_TF, n_min_genes=MIN_MODULE_SIZE)
    mod_obj.get_modules()

    ## 3. Create Te_regulons_maker object and launch multi-processing to compute regulons
    te_reg_mkr = Te_regulons_maker(mod_obj.modules, TE_COUNTS_FILE, TE_ENRICH_FILE, OUT_DIR)
    print("Launching multi-processing for regulon optimisation and generation")
    with Pool(N_CPU) as p:
        adjs = list(tqdm.tqdm(p.imap(te_reg_mkr.analyze_single_regulon, 
                                     te_reg_mkr.tfs),
                              total = len(te_reg_mkr.tfs),
                              position=0, leave=True))
