import argparse
from Bio import SeqIO
from change_seq_length_distributions import convert_to_specified_length


def process_fasta_to_specified_length(length, input_path, output_path):
    assert(length % 2 == 0)
    records = list(SeqIO.parse(input_path, 'fasta'))

    # Convert sequences to given length
    for i, record in enumerate(records):
        records[i].seq = convert_to_specified_length(record.seq, length)

    # Elongate each sequence to 2 kb, by padding with Ns
    for i, record in enumerate(records):
            records[i].seq = convert_to_specified_length(record.seq, 2000) + '\n'

    # Write output to file
    with open(output_path, 'w') as output_file:
        SeqIO.write(records, output_file, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process a FASTA file to a specified sequence length and pad sequences to 2 kb.')
    parser.add_argument('length', type=int, help='Length to which sequences should be converted (must be even).')
    parser.add_argument('input_path', type=str, help='Path to the input FASTA file.')
    parser.add_argument('output_path', type=str, help='Path to the output FASTA file.')

    args = parser.parse_args()
    process_fasta_to_specified_length(args.length, args.input_path, args.output_path)

