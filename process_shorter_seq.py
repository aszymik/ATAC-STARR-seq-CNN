import argparse
import numpy as np
import random
from scipy.interpolate import interp1d
from Bio import SeqIO


def process_fasta_to_list(input_path):
    records = list(SeqIO.parse(input_path, 'fasta'))
    lines = []
    for record in records:
        lines.append(f'>{record.id}\n')
        lines.append(str(record.seq))
    return lines


def cut_sequence(seq, length):
    assert length <= len(seq)
    n = (len(seq) - length) // 2
    return seq[n:n + length]


def elongate_sequence(seq, length):
    assert length >= len(seq)
    n = (2000 - len(seq)) // 2

    if len(seq) % 2 == 1:
        seq = n * 'N' + seq + (n + 1) * 'N'
    else:
        seq = n * 'N' + seq + n * 'N'
    return seq


def convert_to_specified_length(seq, length=2000):
    if len(seq) > length:
        seq = cut_sequence(seq, length)
    else:
        seq = elongate_sequence(seq, length)
    return seq


def process_fasta_to_specified_length(length, input_path, output_path, lengths=[]):
    assert(length % 2 == 0)
    lines = process_fasta_to_list(input_path)

    for i in range(len(lines)):
        if lines[i][0] != '>':
            lengths.append(len(lines[i]))
            # Convert sequences to given length
            lines[i] = convert_to_specified_length(lines[i], length)

    # Elongate each sequence to 2 kb
    for i in range(len(lines)):
        if lines[i][0] != '>':
            lines[i] = convert_to_specified_length(lines[i], 2000) + '\n'

    # Write output
    with open(output_path, 'w') as output_file:
        for line in lines:
            output_file.write(line)


def original_length_to_2kb(input_path, output_path, min_length=1):
    filelines = process_fasta_to_list(input_path)
    lines = []

    for i in range(len(filelines)):
        if filelines[i][0] != '>':

            if len(filelines[i]) >= min_length:
                # Accept only the sequences at least 200 bp long
                lines.append(filelines[i - 1])
                lines.append(convert_to_specified_length(filelines[i], 2000) + '\n')

    with open(output_path, 'w') as output_file:
        for line in lines:
            output_file.write(line)


def filter_length(length, input_path):
    lines = []
    filtered = []

    f = open(input_path, 'r')
    filelines = f.read().split('>')
    f.close()

    print(f'Before filtering: {len(filelines)}')

    lines = process_fasta_to_list(input_path)
    for i in range(len(lines)):
        if len(lines[i][1]) >= length:
            filtered.append(lines[i])

    print(f'After filtering: {len(filtered)}')


def get_lengths(input_path, min_length=1):
    lengths = []
    # Process fasta so that each sequence is in one line
    lines = process_fasta_to_list(input_path)
    beg = 1 if lines[0] == '' else 0
    end = len(lines) - 1 if lines[-1] == '' else len(lines)

    # Write out sequence lengths
    for i in range(beg, end):
        if lines[i][0] != '>':
            if len(lines[i]) >= min_length:
                if len(lines[i]) > 2000:
                    lengths.append(2000)
                else:
                    lengths.append(len(lines[i]))
    return lengths


def get_icdf(lengths):
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


def blocked_gauss(mu, sigma, n):
    numbers = []
    while len(numbers) < n:
        numb = random.gauss(mu, sigma)
        if (numb > 0 and numb < 1):
            numbers.append(numb)
    return numbers


def get_sample_from_distribution(n, lengths, exact=False, seed=0):
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


def create_subsets(input_path, output_path, size=None, lengths=None):
    records = list(SeqIO.parse(input_path, "fasta"))
    data = [record for record in records]

    if not size:
        size = len(data) if lengths is None else len(lengths)

    all_indices = np.arange(len(data))
    np.random.seed(0)
    indices = np.random.choice(all_indices, size=size, replace=False)
    indices.sort()

    testset = []
    for i in range(size):
        if lengths is None:
            testset.append(data[indices[i]])
        else:
            seq = convert_to_specified_length(str(data[indices[i]].seq), lengths[i])
            seq = convert_to_specified_length(seq, 2000)
            testset.append(data[indices[i]])
            testset[-1].seq = seq

    with open(output_path, 'w') as output_file:
        SeqIO.write(testset, output_file, "fasta")


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
    for input_path in reference_paths:
        original_length_to_2kb(input_path, output_path, min_length)
        lengths = get_lengths(input_path, min_length)
        combined_lengths.extend(lengths)

    if subset_size:
        subset_lengths = get_sample_from_distribution(subset_size, combined_lengths, exact=exact)
        create_subsets(input_path, output_path, size=subset_size, lengths=subset_lengths)
    else:
        process_fasta_to_specified_length(length, input_path, output_path)

