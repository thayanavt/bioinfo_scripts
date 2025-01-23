from Bio import SeqIO
import argparse


def extract_region(start, end, genome_path, scaffold_id=None):
    start, end = int(start), int(end)
    reverse_complement = start > end

    if reverse_complement:
        start, end = end, start

    with open(genome_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if scaffold_id is None or record.id == scaffold_id:
                sequence = record.seq[start - 1:end]
                if reverse_complement:
                    sequence = sequence.reverse_complement()
                print(f'>{record.id}:{start}-{end}{" (reverse complement)" if reverse_complement else ''}\n{sequence}')
                break

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--start', required=True, help='Start position (1-based)')
    parser.add_argument('-e', '--end', required=True, help='End position (1-based)')
    parser.add_argument('-g', '--genome', required=True, help='Path to the genome FASTA file')
    parser.add_argument('--scaffold', help='Scaffold ID for multifasta files. If it')

    args = parser.parse_args()
    extract_region(args.start, args.end, args.genome, args.scaffold)

if __name__ == '__main__':
    main()
