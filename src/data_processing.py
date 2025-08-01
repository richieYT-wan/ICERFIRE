import copy

import pandas as pd
import numpy as np
import multiprocessing
import math
from src.utils import pkl_load
import os
import warnings
# import peptides

warnings.filterwarnings('ignore')

DATADIR = '../data/' if os.path.exists('../data/') else './data/'
OUTDIR = '../output/' if os.path.exists('../output/') else './output/'


def _init(DATADIR):
    VAL = math.floor(4 + (multiprocessing.cpu_count() / 1.5))
    N_CORES = VAL if VAL <= multiprocessing.cpu_count() else int(multiprocessing.cpu_count() - 2)

    MATRIXDIR = f'{DATADIR}Matrices/'
    ICSDIR = f'{DATADIR}ic_dicts/'
    AA_KEYS = [x for x in 'ARNDCQEGHILKMFPSTWYV']

    CHAR_TO_INT = dict((c, i) for i, c in enumerate(AA_KEYS))
    INT_TO_CHAR = dict((i, c) for i, c in enumerate(AA_KEYS))

    CHAR_TO_INT['-'] = -1
    INT_TO_CHAR[-1] = '-'

    BG = np.loadtxt(f'{MATRIXDIR}bg.freq.fmt', dtype=float)
    BG = dict((k, v) for k, v in zip(AA_KEYS, BG))

    # BLOSUMS 50
    BL50 = {}
    _blosum50 = np.loadtxt(f'{MATRIXDIR}BLOSUM50', dtype=float).T
    for i, letter_1 in enumerate(AA_KEYS):
        BL50[letter_1] = {}
        for j, letter_2 in enumerate(AA_KEYS):
            BL50[letter_1][letter_2] = _blosum50[i, j]
    BL50_VALUES = {k: np.array([v for v in BL50[k].values()]) for k in BL50}
    # BLOSUMS 62
    BL62_DF = pd.read_csv(f'{MATRIXDIR}BLOSUM62', sep='\s+', comment='#', index_col=0)
    BL62 = BL62_DF.to_dict()
    BL62_VALUES = BL62_DF.drop(columns=['B', 'Z', 'X', '*'], index=['B', 'Z', 'X', '*'])
    BL62_VALUES = dict((x, BL62_VALUES.loc[x].values) for x in BL62_VALUES.index)

    # BLOSUMS 62 FREQS
    _blosum62 = np.loadtxt(f'{MATRIXDIR}BLOSUM62.freq_rownorm', dtype=float).T
    BL62FREQ = {}
    BL62FREQ_VALUES = {}
    for i, letter_1 in enumerate(AA_KEYS):
        BL62FREQ[letter_1] = {}
        BL62FREQ_VALUES[letter_1] = _blosum62[i]
        for j, letter_2 in enumerate(AA_KEYS):
            BL62FREQ[letter_1][letter_2] = _blosum62[i, j]
    HLAS = pkl_load(ICSDIR + 'ics_shannon.pkl')[9].keys()

    return VAL, N_CORES, DATADIR, AA_KEYS, CHAR_TO_INT, INT_TO_CHAR, BG, BL62FREQ, BL62FREQ_VALUES, BL50, BL50_VALUES, BL62, BL62_VALUES, HLAS


# Bad practice but hard-coded way to have all the matrices etc. I need
VAL, N_CORES, DATADIR, AA_KEYS, CHAR_TO_INT, INT_TO_CHAR, BG, BL62FREQ, BL62FREQ_VALUES, BL50, BL50_VALUES, BL62, BL62_VALUES, HLAS = _init(DATADIR)


######################################
####      assertion / checks      ####
######################################

def verify_df(df, seq_col, hla_col, target_col):
    df = copy.deepcopy(df)
    # unique_labels = sorted(df[target_col].dropna().unique())
    # # Checks binary label
    # assert ([int(x) for x in sorted(unique_labels)]) in [[0, 1], [0], [1]], f'Labels are not 0, 1! {unique_labels}'
    # Checks if any seq not in alphabet
    try:
        df = df.drop(df.loc[df[seq_col].apply(lambda x: any([z not in AA_KEYS and not z == '-' for z in x]))].index)
    except:
        print(len(df), df.columns, seq_col, AA_KEYS)
        raise ValueError
    # Checks if HLAs have correct format
    if all(df[hla_col].apply(lambda x: not x.startswith('HLA-'))):
        df[hla_col] = df[hla_col].apply(lambda x: 'HLA-' + x)
    df[hla_col] = df[hla_col].apply(lambda x: x.replace('*', '').replace(':', ''))
    # Check HLA only in subset
    try:
        df = df.query(f'{hla_col} in @HLAS')
    except:
        print(type(df), type(HLAS), HLAS, hla_col)
        raise ValueError(f'{type(df)}, {type(HLAS)}, {HLAS}, {hla_col}, {df[hla_col].unique()}')

    return df


def assert_encoding_kwargs(encoding_kwargs, mode_eval=False):
    """
    Assertion / checks for encoding kwargs and verify all the necessary key-values
    are in
    """
    # Making a deep copy since dicts are mutable between fct calls
    encoding_kwargs = copy.deepcopy(encoding_kwargs)
    if encoding_kwargs is None:
        encoding_kwargs = {'max_len': 12,
                           'encoding': 'onehot',
                           'blosum_matrix': None,
                           'standardize': False}
    essential_keys = ['max_len', 'encoding', 'blosum_matrix', 'standardize']
    keys_check = [x in encoding_kwargs.keys() for x in essential_keys]
    keys_check_dict = {k: v for (k, v) in zip(essential_keys, keys_check) if v == False}
    assert all(keys_check), f'Encoding kwargs don\'t contain the essential key-value pairs! ' \
                            f"{list(keys_check_dict.keys())} are missing!"

    if mode_eval:
        if any([(x not in encoding_kwargs.keys()) for x in ['seq_col', 'hla_col', 'target_col', 'rank_col']]):
            if 'seq_col' not in encoding_kwargs.keys():
                encoding_kwargs.update({'seq_col': 'icore_mut'})
            if 'hla_col' not in encoding_kwargs.keys():
                encoding_kwargs.update({'hla_col': 'HLA'})
            if 'target_col' not in encoding_kwargs.keys():
                encoding_kwargs.update({'target_col': 'agg_label'})
            if 'rank_col' not in encoding_kwargs.keys():
                encoding_kwargs.update({'rank_col': 'EL_rank_mut'})

        # This KWARGS not needed in eval mode since I'm using Pipeline and Pipeline
        del encoding_kwargs['standardize']
    return encoding_kwargs


######################################
####      SEQUENCES ENCODING      ####
######################################

#
# def get_aa_properties(df, seq_col='icore_mut', do_vhse=True, prefix=''):
#     """
#     Compute some AA properties that I have selected
#     keep = ['aliphatic_index', 'boman', 'hydrophobicity',
#         'isoelectric_point', 'VHSE1', 'VHSE3', 'VHSE7', 'VHSE8']
#     THIS KEEP IS BASED ON SOME FEATURE DISTRIBUTION AND CORRELATION ANALYSIS
#     Args:
#         df (pandas.DataFrame) : input dataframe, should contain at least the peptide sequences
#         seq_col (str) : column name containing the peptide sequences
#
#     Returns:
#         out (pandas.DataFrame) : The same dataframe but + the computed AA properties
#
#     """
#     out = df.copy()
#
#     out[f'{prefix}aliphatic_index'] = out[seq_col].apply(lambda x: peptides.Peptide(x).aliphatic_index())
#     out[f'{prefix}boman'] = out[seq_col].apply(lambda x: peptides.Peptide(x).boman())
#     out[f'{prefix}hydrophobicity'] = out[seq_col].apply(lambda x: peptides.Peptide(x).hydrophobicity())
#     out[f'{prefix}isoelectric_point'] = out[seq_col].apply(lambda x: peptides.Peptide(x).isoelectric_point())
#     # out['PD2'] = out[seq_col].apply(lambda x: peptides.Peptide(x).physical_descriptors()[1])
#     # out['charge_7_4'] = out[seq_col].apply(lambda x: peptides.Peptide(x).charge(pH=7.4))
#     # out['charge_6_65'] = out[seq_col].apply(lambda x: peptides.Peptide(x).charge(pH=6.65))
#     if do_vhse:
#         vhse = out[seq_col].apply(lambda x: peptides.Peptide(x).vhse_scales())
#         # for i in range(1, 9):
#         #     out[f'VHSE{i}'] = [x[i - 1] for x in vhse]
#         for i in [1, 3, 7, 8]:
#             out[f'VHSE{i}'] = [x[i - 1] for x in vhse]
#
#     # Some hardcoded bs
#     return out, ['aliphatic_index', 'boman', 'hydrophobicity',
#                  'isoelectric_point', 'VHSE1', 'VHSE3', 'VHSE7', 'VHSE8']


def encode(sequence, max_len=None, encoding='onehot', blosum_matrix=None):
    """
    encodes a single peptide into a matrix, using 'onehot' or 'blosum'
    if 'blosum', then need to provide the blosum dictionary as argument
    """
    assert (encoding == 'onehot' or encoding.lower().startswith("blosum")), 'wrong encoding type'
    encoded = np.empty((20,))
    # One hot encode by setting 1 to positions where amino acid is present, 0 elsewhere
    size = len(sequence)
    if encoding == 'onehot':
        int_encoded = [CHAR_TO_INT[char] for char in sequence]
        onehot_encoded = list()
        for value in int_encoded:
            letter = [0 for _ in range(len(AA_KEYS))]
            letter[value] = 1 if value != -1 else 0
            onehot_encoded.append(letter)
        encoded = np.array(onehot_encoded)

    # BLOSUM encode
    if encoding == 'blosum':
        if blosum_matrix is None or not isinstance(blosum_matrix, dict):
            raise Exception('No BLOSUM matrix provided!')

        encoded = np.zeros([size, len(AA_KEYS)], dtype=np.float32)
        for idx in range(size):
            encoded[idx, :] = blosum_matrix[sequence[idx]]

    # Padding if max_len is provided
    if max_len is not None and max_len > size:
        diff = int(max_len) - int(size)
        try:
            encoded = np.concatenate([encoded, np.zeros([diff, len(AA_KEYS)], dtype=np.float32)],
                                 axis=0)
        except:
            print('Here in encode', type(encoded), encoded.shape, len(AA_KEYS), type(diff), type(max_len), type(size), sequence)
            #     return tmp, diff, len(AA_KEYS)
            raise Exception
    return encoded


def encode_batch(sequences, max_len=None, encoding='onehot', blosum_matrix=None):
    """
    Encode multiple sequences at once.
    """
    if max_len is None:
        max_len = max([len(x) for x in sequences])

    return np.stack([encode(seq, max_len, encoding, blosum_matrix) for seq in sequences]).astype(np.float32)


def get_ic_weights(df, ics_dict: dict, max_len=None, seq_col='Peptide', hla_col='HLA', mask=False,
                   invert=False, threshold=0.234):
    """

    Args:
        df:
        ics_dict:
        max_len:
        seq_col:
        hla_col:
        invert: Invert the behaviour; for KL/Shannon, will take IC instead of 1-IC as weight
                For Mask, will amplify MIA positions (by 1.3) instead of setting to 0

    Returns:

    """
    # if 'len' not in df.columns:

    df['len'] = df[seq_col].apply(len)
    if max_len is not None:
        df = df.query('len<=@max_len')
    else:
        max_len = df['len'].max()
    # Weighting the encoding wrt len and HLA
    lens = df['len'].values
    pads = [max_len - x for x in lens]
    hlas = df[hla_col].str.replace('*', '').str.replace(':', '').values
    # If mask is true, then the weight is just a 0-1 mask filter
    # Using the conserved / MIAs positions instead of the ICs
    if mask:
        # Get mask for where the values should be thresholded to 0 and 1
        weights = np.stack([np.pad(ics_dict[l][hla][0.25], pad_width=(0, pad), constant_values=(1, 1)) \
                            for l, hla, pad in zip(lens, hlas, pads)])
        # IC > 0.3 goes to 0 because anchor position
        # IC <= 0.3 goes to 1 because "MIA" position
        idx_min = (weights > threshold)
        idx_max = (weights <= threshold)
        if invert:
            weights[idx_min] = 1
            weights[idx_max] = 0
        else:
            weights[idx_min] = 0
            weights[idx_max] = 1

    else:
        if invert:  # If invert, then uses the actual IC as weight
            weights = np.stack([np.pad(ics_dict[l][hla][0.25], pad_width=(0, pad), constant_values=(0, 1)) \
                                for l, hla, pad in zip(lens, hlas, pads)])
        else:  # Else we get the weight with the 1-IC depending on the IC dict provided
            weights = 1 - np.stack([np.pad(ics_dict[l][hla][0.25], pad_width=(0, pad), constant_values=(1, 1)) \
                                    for l, hla, pad in zip(lens, hlas, pads)])

    weights = np.expand_dims(weights, axis=2).repeat(len(AA_KEYS), axis=2)
    return weights.astype(np.float32)


# Here stuff for extra AA bulging out:

def find_extra_aa(core, icore):
    """
    Finds the bulging out AA between an icore and its corresponding core, returning the extra AA as "frequencies"
    Args:
        core:
        icore:

    Returns:

    """
    assert len(core) == 9, f'Core is not of length 9 somehow: {core}'
    if len(icore) == len(core) or len(icore) == 8:
        return np.zeros((20)), np.array(0)

    elif len(icore) > len(core):
        results = []
        j = 0
        for i, char in enumerate(icore):
            if char != core[j]:
                results.append(char)
            else:
                j += 1
        # Here changed to len icore - len core to get len of bulge
        # return (encode(''.join(results)).sum(axis=0) / (len(icore)-len(core))).astype(np.float32)

        # Here, changed to return the extra + the length so that we can do the weighted division
        return encode(''.join(results)).sum(axis=0), np.array(len(icore) - len(core))


def batch_find_extra_aa(core_seqs, icore_seqs):
    """
    Same as above but by batch
    Args:
        core_seqs:
        icore_seqs:

    Returns:

    """
    mapped = list(map(find_extra_aa, core_seqs, icore_seqs))
    encoded, lens = np.array([x[0] for x in mapped]), np.array([x[1] for x in mapped])
    return encoded, lens


def encode_batch_weighted(df, ics_dict=None, max_len=None, encoding='onehot', blosum_matrix=None, seq_col='Peptide',
                          hla_col='HLA', mask=False, invert=False, threshold=.200):
    """
    Takes as input a df containing sequence, len, HLA;
    Batch onehot-encode all sequences & weights them with (1-IC) depending on the ICs dict given

    Args:
        df (pandas.DataFrame): DF containing pep sequence, HLA, optionally 'len'
        ics_dict (dict): Dictionary containing the ICs
        max_len (int): Maximum length to consider
        encoding (str) : 'onehot' or 'blosum'
        blosum_matrix : The blosum matrix dictionary; Should just use the BL62_VALUES that's initialized by default
        seq_col (str): Name of the column containing the Peptide sequences (default = 'Peptide')
        hla_col (str): Name of the column containing the HLA alleles (default = 'HLA')

    Returns:
        weighted_sequence (numpy.array): Tensor containing the weighted onehot-encoded peptide sequences.
    """

    if seq_col == 'expanded_input':
        df['seq_len'] = df[seq_col].apply(lambda x: len(x) - x.count('-'))
    else:
        df['seq_len'] = df[seq_col].apply(len)
    if max_len is not None:
        df = df.query('seq_len<=@max_len')
    else:
        max_len = df['seq_len'].max()

    # Encoding the sequences
    encoded_sequences = encode_batch(df[seq_col].values, max_len, encoding=encoding, blosum_matrix=blosum_matrix)
    if ics_dict is not None:
        weights = get_ic_weights(df, ics_dict, max_len, seq_col, hla_col, mask, invert, threshold)
    else:
        # Here, if no ics_dict is provided, the normal weight will just be ones everywhere
        # In case we are not doing weighted sequence (either for onehot-input or frequency computation)
        weights = np.ones(encoded_sequences.shape)
    weighted_sequences = weights * encoded_sequences
    true_lens = df['seq_len'].values
    return weighted_sequences.astype(np.float32), true_lens


def get_dataset(df, ics_dict, max_len=12, encoding='onehot', blosum_matrix=None, seq_col='icore_mut', hla_col='HLA',
                target_col='agg_label', rank_col='EL_rank_mut', mut_col=None, mask=False, invert=False, add_rank=False,
                threshold=.200):
    """
    """
    # df = verify_df(df, seq_col, hla_col, target_col)

    encoded_weighted, true_lens = encode_batch_weighted(df, ics_dict, max_len, encoding, blosum_matrix, seq_col,
                                                        hla_col, mask, invert, threshold)
    x = batch_compute_frequency(encoded_weighted, true_lens)
    if add_rank:
        ranks = np.expand_dims(df[rank_col].values, 1)
        try:
            x = np.concatenate([x, ranks], axis=1)
        except:
            print(rank_col, ranks.shape, x.shape)
            print('\n\n\n', x[:5])
            print('\n\n\n', ranks[:5])
            raise ValueError

    y = df[target_col].values

    if mut_col is not None and type(mut_col) == list:
        if len(mut_col) > 0:
            mut_scores = df[mut_col].values
            x = np.concatenate([x, mut_scores], axis=1)

    return x, y


def get_test_dataset(df, ics_dict, max_len=12, encoding='onehot', blosum_matrix=None, seq_col='icore_mut', hla_col='HLA',
                     target_col='agg_label', rank_col='EL_rank_mut', mut_col=None, mask=False, invert=False, add_rank=False,
                     threshold=.200):
    """
    """
    # df = verify_df(df, seq_col, hla_col, target_col)

    encoded_weighted, true_lens = encode_batch_weighted(df, ics_dict, max_len, encoding, blosum_matrix, seq_col,
                                                        hla_col, mask, invert, threshold)
    x = batch_compute_frequency(encoded_weighted, true_lens)
    if add_rank:
        ranks = np.expand_dims(df[rank_col].values, 1)
        try:
            x = np.concatenate([x, ranks], axis=1)
        except:
            print(rank_col, ranks.shape, x.shape)
            print('\n\n\n', x[:5])
            print('\n\n\n', ranks[:5])
            raise ValueError

    if mut_col is not None and type(mut_col) == list:
        if len(mut_col) > 0:
            mut_scores = df[mut_col].values
            x = np.concatenate([x, mut_scores], axis=1)

    return x


def batch_compute_frequency(encoded_sequences, true_lens=None):
    """

    Args:
        encoded_sequences:
    Returns:

    """
    # This is the new way with mask and .all(dim=2) which works with both BLOSUM and OH
    if true_lens is None:
        mask = (encoded_sequences == 0).all(2)  # checking on second dim that every entry == 0
        # true_lens = (mask.shape[1] - torch.bincount(torch.where(mask)[0])).unsqueeze(1) if type(
        #     mask) == torch.Tensor else \
        #     np.expand_dims(mask.shape[1] - np.bincount(np.where(mask)[0]), 1)
        true_lens = np.expand_dims(mask.shape[1] - np.bincount(np.where(mask)[0]), 1)
        frequencies = encoded_sequences.sum(axis=1) / true_lens
    else:
        frequencies = encoded_sequences.sum(axis=1) / np.repeat(true_lens, axis=0, repeats=20).reshape(len(true_lens),
                                                                                                       20)

    return frequencies


#### ==== PSSM, pfm etc ==== ####
def get_weights(onehot_seqs):
    """
    Compute the sequences weights using heuristics (1/rs)
    :param onehot_seqs:
    :param counts:
    :return:
    """
    # Get counts
    counts = onehot_seqs.sum(axis=0)
    # Get absolute counts (i.e. # of diff AA in position K --> r)
    abs_counts = counts.copy()
    abs_counts[abs_counts > 0] = 1
    rs = abs_counts.sum(axis=1)
    # Get total count of each aa per position --> s
    ss = (onehot_seqs * counts).sum(axis=2)
    # weights = 1/sum(r*s)
    weights = (1 / (np.multiply(rs, ss))).sum(axis=1)
    n_eff = rs.sum() / onehot_seqs.shape[1]
    # Reshaping to fit the right shape to allow pointwise mul with onehot
    weights = np.expand_dims(np.tile(weights, (onehot_seqs.shape[1], 1)).T, axis=2).repeat(20, axis=2)
    return weights, n_eff


def compute_pfm(sequences, how='shannon', seq_weighting=False, beta=50):
    """
    Computes the position frequency matrix or pseudofrequency given a list of sequences
    """
    try:
        max_len = max([len(x) for x in sequences])
    except:
        print(sequences)
        raise Exception(sequences)
    N = len(sequences)
    onehot_seqs = encode_batch(sequences, max_len, encoding='onehot', blosum_matrix=None)

    if how == 'shannon':
        freq_matrix = onehot_seqs.sum(axis=0) / N
        return freq_matrix

    elif how == 'kl':
        weights, neff = get_weights(onehot_seqs) if seq_weighting else (1, len(sequences))
        # return weights, neff
        onehot_seqs = weights * onehot_seqs
        alpha = neff - 1
        freq_matrix = onehot_seqs.sum(axis=0) / N
        g_matrix = np.matmul(BL62FREQ, freq_matrix.T).T
        p_matrix = (alpha * freq_matrix + beta * g_matrix) / (alpha + beta)
        return p_matrix


def compute_ic_position(matrix, position):
    """

    Args:
        matrix:
        position:

    Returns:

    """
    row = matrix[position]
    row_log20 = np.nan_to_num(np.log(row) / np.log(20), neginf=0)
    ic = 1 + np.sum(row * row_log20)
    return ic


def compute_ic(sequences, how='shannon', seq_weighting=True, beta=50):
    """
    returns the IC for sequences of a given length based on the frequency matrix
    Args:
        sequences (list) : list of strings (sequences) from which to compute the IC
        how (str): 'shannon' or 'kl' for either shannon or kullback leibler PFM
    Returns:
        ic_array (np.ndarray) : A Numpy array of the information content at each position (of shape max([len(seq) for seq in sequences]), 1)
    """
    # if type(sequences) == np.ndarray:
    #     return np.array([compute_ic_position(sequences, pos) for pos in range(sequences.shape[0])])
    pfm = compute_pfm(sequences, how, seq_weighting, beta)
    ic_array = np.array([compute_ic_position(pfm, pos) for pos in range(pfm.shape[0])])
    return ic_array


def get_mia(ic_array, threshold=.2):
    return np.where(ic_array < threshold)[0]
