import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def pad_with_Ns(seq, output_length=2000):
    assert output_length >= len(seq)
    n = (output_length - len(seq)) // 2

    if len(seq) % 2 == 1:
        seq = n*'N' + seq + (n+1)*'N'
    else:
        seq = n*'N' + seq + n*'N'
    return seq


def pad_fasta_and_filter_shortest(input_path, output_path, min_length=1, output_length=2000):
    records_filtered = []

    for record in SeqIO.parse(input_path, 'fasta'):
        # Accept only the sequences longer than min_length
        if len(record.seq) >= min_length:
            seq = pad_with_Ns(str(record.seq), output_length=output_length)
            records_filtered.append(record)
            records_filtered[-1].seq = Seq(seq)

    with open(output_path, 'w') as output_file:
        SeqIO.write(records_filtered, output_file, 'fasta')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pad sequences in a FASTA file with 'N's to a specified length and filter by minimum length.")
    parser.add_argument("input_path", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_path", type=str, help="Path to the output FASTA file.")
    parser.add_argument("--min_length", type=int, default=1, help="Minimum length of sequences to be included in the output file. Default is 1.")
    parser.add_argument("--output_length", type=int, default=2000, help="Length to which sequences will be padded. Default is 2000.")
    
    args = parser.parse_args()
    pad_fasta_and_filter_shortest(args.input_path, args.output_path, args.min_length, args.output_length)
