import warnings
warnings.simplefilter(action='ignore', category=DeprecationWarning)

import re
from collections import Counter, OrderedDict
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import scipy.spatial.distance as distance
from functools import partial

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher, ParsimonyTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

# Column name mappings - set these to match your input data format
COLUMN_NAMES = {
    'aa_change': 'HGVSp_Short',  # Standard name = actual name in data
    'alt_depth': 't_alt_count',
    'total_depth': 't_depth',
    'sample_id': 'Tumor_Sample_Barcode',
    'chrom': 'Chromosome',
    'start_pos': 'Start_Position',
    'end_pos': 'End_Position',
    'ref_allele': 'Reference_Allele',
    'alt_allele': 'Tumor_Seq_Allele2',
    'vaf': 'VAF',
    'gene': 'Hugo_Symbol',
    'variant_class': 'Variant_Classification'
}

# Configure matplotlib for retina display
plt.rcParams['figure.dpi'] = 150

# Get the directory containing this file
PACKAGE_DIR = Path(__file__).parent.parent

def data_factory(query):
    """Load mutation data based on query type.
    
    Args:
        query (str): One of "01full", "03full", or "full"
        
    Returns:
        pd.DataFrame: Filtered mutation data
    """
    # Load the single MAF file
    fn = PACKAGE_DIR / "data/v3/ALL_somatic.snv.filter.passed.variantfilter.ann.vep.maf"
    df = pd.read_table(fn, low_memory=False)
    
    # Drop columns that are completely NA
    df = df.dropna(axis=1, how='all')
    
    # Calculate VAF
    df[COLUMN_NAMES['vaf']] = df[COLUMN_NAMES['alt_depth']] / df[COLUMN_NAMES['total_depth']]
    
    # Filter out specific samples
    exclude_samples = ['FAP01_P2_B1', 'FAP01_P6_B1', 'FAP01_P6_R5_G1_L2', 'FAP01_T3_B1', 'FAP03_N2_R4_G6']
    df = df[~df[COLUMN_NAMES['sample_id']].isin(exclude_samples)]
    
    match query:
        case "01full":
            return df[df[COLUMN_NAMES['sample_id']].str.startswith('FAP01')]
        case "03full":
            return df[df[COLUMN_NAMES['sample_id']].str.startswith('FAP03')]
        case "full":
            return df
        case _:
            raise NotImplementedError

class SampleSet:
    """Class for handling and analyzing sample sets."""
    var_cols = [COLUMN_NAMES['chrom'], COLUMN_NAMES['start_pos'], COLUMN_NAMES['end_pos']]
    sample_col = COLUMN_NAMES['sample_id']

    def __init__(self, tbl, prefix):
        self.tbl = tbl
        self.prefix = prefix
        self.set_sample_df()
        
        self._vafdm = None
        self._alignment = None
        self._vaf_alignment = None
        self._alignment_type = 'binary'  # or 'vaf'
        self._freq_tree = None # Not yet cached
        self._seq_tree = None # Not yet cached
        self._pars_tree = None # Not yet cached
        self.colordict = {}
        self.vaf_thresholds = {'present': 0.3, 'absent': 0.1}  # Default thresholds

    def set_sample_df(self):
        """Filter dataframe to include only samples with specified prefix."""
        self.sample_df = self.tbl[self.tbl[self.sample_col].str.startswith(self.prefix)]

    def get_dims(self):
        """Get dimensions of the sample set."""
        num_vars = self.sample_df[self.var_cols].drop_duplicates().shape[0]
        num_samples = self.sample_df[self.sample_col].drop_duplicates().shape[0]
        return num_vars, num_samples

    def get_sample_table(self):
        """Get sample table with unique sample IDs."""
        return self.sample_df[[self.sample_col]].drop_duplicates().reset_index(drop=True).reset_index().set_index(self.sample_col)

    @property
    def vafdm(self):
        """Get VAF distance matrix."""
        if self._vafdm is None:
            self._vafdm = self.get_vafdm()
        return self._vafdm

    @property
    def alignment(self):
        """Get sequence alignment based on current alignment type."""
        if self._alignment_type == 'binary':
            if self._alignment is None:
                self._alignment = self.get_binary_alignment()
            return self._alignment
        else:  # vaf-based
            if self._vaf_alignment is None:
                self._vaf_alignment = self.get_vaf_alignment()
            return self._vaf_alignment

    def set_alignment_type(self, alignment_type, vaf_thresholds=None):
        """Set the type of alignment to use and optionally update VAF thresholds.
        
        Args:
            alignment_type (str): Either 'binary' or 'vaf'
            vaf_thresholds (dict, optional): Dictionary with 'present' and 'absent' thresholds
        """
        if alignment_type not in ['binary', 'vaf']:
            raise ValueError("alignment_type must be either 'binary' or 'vaf'")
        
        self._alignment_type = alignment_type
        if vaf_thresholds is not None:
            self.vaf_thresholds = vaf_thresholds
        
        # Clear cached alignments when switching types
        self._alignment = None
        self._vaf_alignment = None

    def get_binary_alignment(self):
        """Generate binary sequence alignment based on presence/absence of alternative alleles."""
        sample_table = self.get_sample_table()
        num_vars, num_samples = self.get_dims()
        matrix = np.zeros((num_samples+1, num_vars), dtype='object')
        for i, (idx, group) in enumerate(self.sample_df.groupby(self.var_cols)):
            ref = group[COLUMN_NAMES['ref_allele']].drop_duplicates().squeeze()

            try:
                assert len(ref) == 1
            except AssertionError:
                # Handle rare case of 1-base insertion at same site as SNV
                ref = ref[ref != "-"].squeeze()
                assert len(ref) == 1
                group = group.query(f'{COLUMN_NAMES["ref_allele"]} != "-"')
            
            matrix[:, i] = ref
            for _, row in group.iterrows():
                sid = sample_table.loc[row[self.sample_col]].squeeze()
                matrix[sid, i] = row[COLUMN_NAMES['alt_allele']]
    
        seq_records = []
        for i, row in enumerate(matrix):
            sequence = ''.join(row)
            seq = Seq(sequence)
            try:
                seq_record = SeqRecord(seq, id=sample_table.index[i])
            except IndexError:
                seq_record = SeqRecord(seq, id="Reference")
            seq_records.append(seq_record)
        return MultipleSeqAlignment(seq_records)

    def get_vaf_alignment(self):
        """Generate VAF-based sequence alignment.
        
        Uses VAF thresholds to determine variant presence:
        - VAF ≥ present_threshold: Use alternative allele
        - VAF ≤ absent_threshold: Use reference allele
        - Otherwise: Use 'N' to indicate uninformative
        """
        sample_table = self.get_sample_table()
        num_vars, num_samples = self.get_dims()
        matrix = np.zeros((num_samples+1, num_vars), dtype='object')
        
        for i, (idx, group) in enumerate(self.sample_df.groupby(self.var_cols)):
            ref = group[COLUMN_NAMES['ref_allele']].drop_duplicates().squeeze()

            try:
                assert len(ref) == 1
            except AssertionError:
                # Handle rare case of 1-base insertion at same site as SNV
                ref = ref[ref != "-"].squeeze()
                assert len(ref) == 1
                group = group.query(f'{COLUMN_NAMES["ref_allele"]} != "-"')
            
            matrix[:, i] = ref
            for _, row in group.iterrows():
                sid = sample_table.loc[row[self.sample_col]].squeeze()
                vaf = row[COLUMN_NAMES['vaf']]
                
                if vaf >= self.vaf_thresholds['present']:
                    matrix[sid, i] = row[COLUMN_NAMES['alt_allele']]
                elif vaf <= self.vaf_thresholds['absent']:
                    matrix[sid, i] = ref
                else:
                    matrix[sid, i] = 'N'
    
        seq_records = []
        for i, row in enumerate(matrix):
            sequence = ''.join(row)
            seq = Seq(sequence)
            try:
                seq_record = SeqRecord(seq, id=sample_table.index[i])
            except IndexError:
                seq_record = SeqRecord(seq, id="Reference")
            seq_records.append(seq_record)
        return MultipleSeqAlignment(seq_records)

    def set_colordict(self, prefixes, colors):
        """Set color dictionary for visualization."""
        sample_table = self.get_sample_table()
        for row in sample_table.itertuples():
            for prefix, color in zip(prefixes, colors):
                if row.Index.startswith(prefix):
                    sample_table.loc[row.Index, "color"] = color
                    break
        self.colordict = sample_table.fillna("black").color.to_dict()
        
    def get_vafdm(self):
        """Calculate VAF distance matrix."""
        sample_table = self.get_sample_table()
        num_vars, num_samples = self.get_dims()
        matrix = np.zeros((num_samples+1, num_vars))
        for i, (idx, group) in enumerate(self.sample_df.groupby(self.var_cols)):
            for _, row in group.iterrows():
                sid = sample_table.loc[row[self.sample_col]].squeeze()
                matrix[sid, i] = row[COLUMN_NAMES['vaf']]
    
        dist = distance.squareform(distance.pdist(matrix)).tolist()
        lower_triangle = []
        for i in range(0, len(dist)):
            lower_triangle.append(dist[i][:i+1])
        return DistanceMatrix(list(sample_table.index) + ["Reference",], lower_triangle)

    def freq_nj_tree(self):
        """Generate neighbor-joining tree based on VAF distances."""
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(self.vafdm)
        tree.root_with_outgroup('Reference')
        return tree

    def seq_nj_tree(self):
        """Generate neighbor-joining tree based on sequence alignment."""
        calculator = DistanceCalculator('identity')
        constructor = DistanceTreeConstructor(calculator, 'nj')
        tree = constructor.build_tree(self.alignment)
        tree.root_with_outgroup('Reference')
    
        # Adjust branch lengths to be in units of number of mutations
        for clade in tree.find_clades():
            clade.branch_length = clade.branch_length * self.alignment.get_alignment_length()
        return tree

    def seq_pars_tree(self):
        """Generate parsimony tree."""
        starting_tree = self.seq_nj_tree()
        scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher, starting_tree)
        tree = constructor.build_tree(self.alignment)
        return tree

    def get_site_statistics(self):
        """Calculate site statistics for the alignment."""
        ref = [seq for seq in self.alignment if seq.id == "Reference"]
        assert len(ref) == 1
        ref = ref[0]
        site_statistics = Counter()
        for i in range(1, self.get_dims()[-1] + 1):
            site_statistics[i] = 0
        for i in range(self.alignment.get_alignment_length()):
            column = self.alignment[:, i]
            base_counts = Counter(column)
            base_counts.pop(ref[i], None)
            if len(base_counts) < 1:
                continue
            for count in base_counts.values():
                # For multiallelic sites, count each ALT individually
                site_statistics[count] += 1
        return site_statistics

    def reset_alignments(self):
        """Reset all cached alignments and derived trees.
        
        This is useful when:
        - Switching alignment types
        - Modifying the underlying data
        - Wanting to force recalculation of alignments
        """
        self._alignment = None
        self._vaf_alignment = None
        self._freq_tree = None
        self._seq_tree = None
        self._pars_tree = None

def plot_tree(tree, color_dict, ax):
    """Plot phylogenetic tree with specified colors."""
    ax.set_frame_on(False)
    ax.get_yaxis().set_ticks([])
    ax.get_yaxis().get_label().set_visible(False)
    Phylo.draw(
        tree,
        label_func=lambda lbl: None if str(lbl).startswith("Inner") else str(lbl),
        axes=ax, do_show=False,
        label_colors=color_dict,
    )

def plot_cmpr_nj(samples: SampleSet):
    """Compare VAF and SNP-based neighbor-joining trees."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    ax1.set_title(f"{samples.prefix} - VAFs")
    plot_tree(samples.freq_nj_tree(), samples.colordict, ax1)

    ax2.set_title(f"{samples.prefix} - SNPs")
    plot_tree(samples.seq_nj_tree(), samples.colordict, ax2)
    
    return fig

def plot_cmpr_pars(samples: SampleSet):
    """Compare neighbor-joining and parsimony trees."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    ax1.set_title(f"{samples.prefix} - NJ")
    plot_tree(samples.seq_nj_tree(), samples.colordict, ax1)

    ax2.set_title(f"{samples.prefix} - Parsimony")
    plot_tree(samples.seq_pars_tree(), samples.colordict, ax2)
    
    return fig

def size_filter(tbl):
    """Filter for single nucleotide variants."""
    return (
        tbl.assign(ref_len=lambda df: df[COLUMN_NAMES['ref_allele']].str.len())
        .assign(alt_len=lambda df: df[COLUMN_NAMES['alt_allele']].str.len())
        .query('(ref_len == 1) & (alt_len == 1)')
    )

def vaf_filter(tbl, vaf_thres=0.25):
    """Filter variants by VAF threshold."""
    return tbl.query(f'{COLUMN_NAMES["vaf"]} >= {vaf_thres}')

def coverage_filter(tbl, min_depth=10, max_depth=100):
    """Filter variants by coverage depth.
    
    Args:
        tbl (pd.DataFrame): Input dataframe
        min_depth (int, optional): Minimum total depth. Defaults to 10.
        max_depth (int, optional): Maximum total depth. Defaults to 100.
        
    Returns:
        pd.DataFrame: Filtered dataframe
    """
    return tbl.query(f'{COLUMN_NAMES["total_depth"]} >= {min_depth} & {COLUMN_NAMES["total_depth"]} <= {max_depth}')

def apply_filters(tbl, vaf_thres=0.25, min_depth=10, max_depth=100):
    """Apply all filters in sequence: size, coverage, then VAF.
    
    Args:
        tbl (pd.DataFrame): Input dataframe
        vaf_thres (float, optional): VAF threshold. Defaults to 0.25.
        min_depth (int, optional): Minimum total depth. Defaults to 10.
        max_depth (int, optional): Maximum total depth. Defaults to 100.
        
    Returns:
        pd.DataFrame: Filtered dataframe
    """
    return (tbl
            .pipe(size_filter)
            .pipe(coverage_filter, min_depth=min_depth, max_depth=max_depth)
            .pipe(vaf_filter, vaf_thres=vaf_thres)
           )

def make_region_tuple(query):
    """Create tuple of region identifiers."""
    seven_regions = [f"R{num}" for num in range(1, 8)]
    return tuple((f"{query}_{region}" for region in seven_regions))

def add_normal(query):
    """Add normal sample information."""
    if query.startswith("FAP01"):
        return (("FAP01_N",), ("FAP01_N",), ("#377abb",))
    elif query.startswith("FAP03"):
        return (("FAP03_N",), ("FAP03_N",), ("#377abb",))
    else:
        raise ValueError

def polyp_factory(query, with_normal=True):
    """Create polyp analysis configuration."""
    seven_color_tuple = ('#ff8c2a', '#bc6ca4', '#fe919f', '#4bb341', '#ff383a', '#67c5be', '#992b00')
    target = ((f"{query}_R",), make_region_tuple(query), seven_color_tuple)

    if with_normal:
        return tuple(map(lambda x, y: x + y, target, add_normal(query)))
    else:
        return target

# Load metadata and driver gene list
def load_metadata():
    """Load sample metadata."""
#    return pd.read_table(PACKAGE_DIR / "data/sgWGS_meta.txt")
    return pd.read_csv(PACKAGE_DIR / "data/TableS1_sgWGS_mutMeta.csv")

def load_driver_genes():
    """Load list of driver genes."""
    return pd.read_csv(PACKAGE_DIR / "data/PanCanDrivers_COADREAD_Cell2018.csv")['Gene']

def get_all_queries():
    """Get list of all tumor sample queries, excluding blood and normal tissue."""
    metadata = load_metadata()
    return metadata.query('~Region.str.startswith("B")').query('~Tissue.str.startswith("N")').Tumor.unique()

def df_filter_targets(df):
    """Filter dataframe to include only driver gene mutations.
    
    Args:
        df (pd.DataFrame): Input dataframe with mutation data
        
    Returns:
        pd.DataFrame: Filtered dataframe containing only driver gene mutations
    """
    driver_list = load_driver_genes()
    return (df
            .query(f'{COLUMN_NAMES["gene"]}.isin(@driver_list)')
            .query(f'{COLUMN_NAMES["variant_class"]}.str.contains("Mutation") | {COLUMN_NAMES["variant_class"]}.str.contains("Frame_Shift")')
            .sort_values([COLUMN_NAMES['gene'], COLUMN_NAMES['sample_id']])
            .assign(mut_str=lambda df: df[COLUMN_NAMES['gene']] +":"+ df[COLUMN_NAMES['aa_change']])
           )

def annotate_tree_with_targets(tree, ax, sample_df, fontsize=6):
    """Annotate a phylogenetic tree with driver gene mutations.
    
    Args:
        tree (Bio.Phylo.BaseTree.Tree): Phylogenetic tree to annotate
        ax (matplotlib.axes.Axes): Axes containing the tree plot
        sample_df (pd.DataFrame): Sample dataframe containing mutation data
        fontsize (int, optional): Font size for annotations. Defaults to 6.
    """
    annot_df = df_filter_targets(sample_df)
    annot_dict = annot_df.groupby(COLUMN_NAMES['sample_id'])['mut_str'].apply(lambda x: '\n'.join(x)).to_dict()
    
    text_labels = [text for text in ax.get_children() if isinstance(text, plt.Text)]
    for text in text_labels:
        leaf_label = text.get_text().strip()
        if leaf_label in annot_dict:
            x, y = text.get_position()
            x = x - tree.find_any(name=leaf_label).branch_length * 0.95
            y = y - 0.1
            annotation_text = annot_dict[leaf_label]
            ax.text(x, y, annotation_text, fontsize=fontsize, verticalalignment='bottom')
