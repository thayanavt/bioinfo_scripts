from Bio.Seq import Seq
from pathlib import Path
from typing import Tuple
import argparse
import pandas as pd

'''Neb file in the context of this script refers to a LAMP primer set that
can be created on the website https://lamp.neb.com/#!/

File content example:
F3	69	85	17	59.97	-5.35	-5.30	65	CCAGACGCCGTTTTCCG
B3	238	256	19	61.43	-7.02	-5.51	63	CGCTGGAGGAGATGACCAC
FIP			41					TCCGTACAGCCCTATGGCACC GCTTCATCCAGTGCCAGTGA
BIP			38					CGCCCCCAGTCCAGCGTTT GCCGTACATTTTCCCGACG
F2	86	105	20	62.05	-5.26	-4.91	55	GCTTCATCCAGTGCCAGTGA
F1c	126	146	21	65.87	-5.27	-6.85	62	TCCGTACAGCCCTATGGCACC
B2	219	237	19	60.58	-6.26	-6.96	58	GCCGTACATTTTCCCGACG
B1c	175	193	19	67.22	-7.97	-5.84	68	CGCCCCCAGTCCAGCGTTT

IMPORTANT: primers B3, F1c, and B2 have their sequences in reverse complement in the final fasta.
'''

def process_primer_row(row: pd.Series) -> Tuple[str, str]:
    sequence = row['seq']
    identifier = row['primer_id']
    if identifier in {'F1c', 'B2', 'B3'}:
        return f"{identifier}_rc", str(Seq(sequence).reverse_complement())
    return identifier, sequence

def read_and_process_neb_file(file_path: Path) -> pd.DataFrame:
    df_neb = pd.read_csv(file_path, sep=r'\s+', header=None).iloc[:, [0, -1]]
    df_neb.columns = ['primer_id', 'seq']
    df_neb[['primer_id', 'seq']] = df_neb.apply(
        lambda x: pd.Series(process_primer_row(x)), axis=1
    )
    return df_neb

def convert_to_fasta(df: pd.DataFrame) -> str:
    df_fasta = df.apply(lambda x: f">{x['primer_id']}\n{x['seq']}", axis=1)
    return '\n'.join(df_fasta)

def save_fasta_file(output_path: Path, fasta_content: str) -> None:
    with open(output_path, 'w') as file_content:
        file_content.write(fasta_content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-neb', '--neb_file', required=True)
    parser.add_argument('-o', '--output_path', required=True)
    args = parser.parse_args()

    neb_file = Path(args.neb_file)
    output_file = Path(args.output_path).joinpath(neb_file.with_suffix('.fna').name)

    df_neb = read_and_process_neb_file(file_path=neb_file)
    multifasta = convert_to_fasta(df=df_neb)
    save_fasta_file(output_path=output_file, fasta_content=multifasta)
