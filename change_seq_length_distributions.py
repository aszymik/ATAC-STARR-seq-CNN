import argparse
import numpy as np
import random
from scipy.interpolate import interp1d
from Bio import SeqIO
from pad_fasta_with_Ns import pad_with_Ns


def count_sequences_in_fasta(input_path):
    return sum(1 for _ in SeqIO.parse(input_path, 'fasta'))


def cut_sequence(seq, length):
    assert length <= len(seq)
    n = (len(seq) - length) // 2
    return seq[n:n + length]


def convert_to_specified_length(seq, length=2000):
    if len(seq) > length:
        seq = cut_sequence(seq, length)
    else:
        seq = pad_with_Ns(seq, length)
    return seq


def get_lengths(input_path, min_length=1):
    """Returns sequence lenghts from a fasta file longer than min_length"""
    lengths = []

    for record in SeqIO.parse(input_path, 'fasta'):
        sequence = str(record.seq)
        trimmed_sequence = sequence.strip('N')  # not counting 'N's padding a sequence
        if len(trimmed_sequence) > 2000:
            lengths.append(2000)
        elif len(trimmed_sequence) > min_length:  
            lengths.append(len(trimmed_sequence))
    return lengths


def get_icdf(lengths):
    """Calculates icdf, divided into sequences longer and shorter than 2000"""
    assert len(lengths) > 0
    under_2000 = [length for length in lengths if length < 2000]
    over_2000_ratio = (len(lengths) - len(under_2000)) / len(lengths)
    under_2000.sort()

    # Define CDF
    cdf = np.cumsum(under_2000) / np.sum(under_2000)

    # Interpolate ICDF
    min_length = min(under_2000)
    print(f'{min_length=}')
    max_length = max(under_2000)
    print(f'{max_length=}')

    icdf_interp = interp1d(cdf, under_2000, kind='cubic', fill_value='extrapolate')

    def icdf(prob):
        # Returns length of the sequence for a given probability
        return icdf_interp(prob)

    return icdf, over_2000_ratio


def get_sample_from_distribution(n, lengths, exact=False, seed=0):
    """ Imitates a given length distribution from lengths list 
        If exact: sampling from the list, else: using ICDF
    """
    if not exact:
        icdf, over_2000_ratio = get_icdf(lengths)
        over_2000_count = int(over_2000_ratio * n)
        under_2000_count = n - over_2000_count

        over_2000 = [2000] * over_2000_count
        np.random.seed(seed)
        under_2000 = [int(numb) for numb in np.round(icdf(np.random.uniform(size=under_2000_count)))]
        chosen_lengths = over_2000 + list(under_2000)
        random.shuffle(chosen_lengths)

        return chosen_lengths

    if n <= len(lengths):
        np.random.seed(seed)
        chosen_lengths = np.random.choice(lengths, size=n, replace=False)
    else:
        np.random.seed(seed)
        chosen_lengths = np.random.choice(lengths, size=n, replace=True)
    return chosen_lengths


def create_subsets(input_path, output_path, size=None, lengths=None, seed=0):
    """ Creates a subset of the given size of data from the FASTA file, 
        if a list of sequence lengths is given, changes the sequence lengths accordingly 
    """
    data = list(SeqIO.parse(input_path, "fasta"))

    if not size:
        size = len(data) if lengths is None else len(lengths)

    # Choose random sequence indices
    all_indices = np.arange(len(data))
    np.random.seed(seed)
    indices = np.random.choice(all_indices, size=size, replace=False)
    indices.sort()

    subset = []
    for i in range(size):
        record = data[indices[i]]
        if lengths is None:
            subset.append(record)
        else:
            seq = convert_to_specified_length(str(record.seq), lengths[i])
            seq = convert_to_specified_length(seq, 2000)
            subset.append(record)
            subset[-1].seq = seq

    with open(output_path, 'w') as output_file:
        SeqIO.write(subset, output_file, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process and manipulate FASTA files.")
    parser.add_argument('--reference', type=str, nargs='+', required=True, help='FASTA file paths of which distribution will be mimicked (at least one required)')
    parser.add_argument('--input', type=str, required=True, help='Input FASTA file path')
    parser.add_argument('--output', type=str, required=True, help='Output FASTA file path')
    parser.add_argument('--length', type=int, default=2000, help='Target length for sequences')
    parser.add_argument('--min_length', type=int, default=1, help='Minimum length for sequences to be included')
    parser.add_argument('--subset_size', type=int, help='Size of the subset to create')
    parser.add_argument('--exact', action='store_true', help='Use exact lengths from distribution')
    args = parser.parse_args()

    reference_paths = args.reference
    input_path = args.input
    output_path = args.output
    length = args.length
    min_length = args.min_length
    subset_size = args.subset_size
    exact = args.exact

    combined_lengths = []
 
    # Process each reference file in terms of sequence lengths
    for ref_path in reference_paths:
        lengths = get_lengths(ref_path, min_length)
        combined_lengths.extend(lengths)

    if not subset_size:
        subset_size = count_sequences_in_fasta(input_path)

    subset_lengths = get_sample_from_distribution(subset_size, combined_lengths, exact=exact)
    create_subsets(input_path, output_path, size=subset_size, lengths=subset_lengths)

