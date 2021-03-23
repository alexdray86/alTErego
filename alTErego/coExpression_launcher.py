import numpy

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
