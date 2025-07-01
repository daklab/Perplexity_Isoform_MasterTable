import csv, typer
from pathlib import Path
from typing import Dict, Tuple, List


def extract_last_junction(psl_input: Path) -> Dict[str, int]:
    with open(psl_input, "rt") as psl:
        read_last = {}
        for line in psl:
            read = line.strip().split("\t")
            name = read[9].lower() # to ignore how cpat changes the letter case
            qstarts = read[19].split(",")[:-1]
            strand = read[8]
            if strand == "+":
                last_junc = int(qstarts[-1])
            if strand == "-":
                qsize = read[10]
                last_junc = int(qsize) - int(qstarts[1])
            read_last[name] = last_junc

    return read_last


def extract_orf_end(cpat_tsv: Path) -> Dict[str, int]:
    with open(cpat_tsv, "rt") as orf_tsv:
        next(orf_tsv)
        orf_end_dict = {}
        for line in orf_tsv:
            row = line.strip().split("\t")
            name = row[0].lower() # to ignore how cpat changes the letter case
            orf_end = row[6]
            orf_end_dict[name] = int(orf_end)

    return orf_end_dict


def filter_NMD_id(
    last_junctions: Dict[str, int], orf_end_dict: Dict[str, int]
) -> Tuple[List[str], List[str]]:
    non_NMD = []
    NMD = []
    for k, v in last_junctions.items():
        if k in orf_end_dict:
            if orf_end_dict[k]+50 < v:
                NMD.append(k)
            else:
                non_NMD.append(k)

    return non_NMD, NMD


def write_filtered_NMD_fasta(
    fa_input: Path, fa_output: Path, psl_input: Path, cpat_tsv: Path
) -> None:
    last_junctions = extract_last_junction(psl_input)
    orf_end_dict = extract_orf_end(cpat_tsv)
    non_NMD, NMD = filter_NMD_id(last_junctions, orf_end_dict)

    with open(fa_input, "rt") as orf_fa:
        with open(fa_output, "wt") as outfile:
            tsv_writer = csv.writer(outfile, delimiter="\t")
            for line in orf_fa:
                row = line.strip()
                if row.startswith(">"):
                    isoform = row.split(">")[-1]
                else:
                    continue
                if isoform.lower() in non_NMD:
                    tsv_writer.writerow([row])
                    tsv_writer.writerow([next(orf_fa, "").strip()])
    for i in NMD:
        print(i)


if __name__ == "__main__":
    typer.run(write_filtered_NMD_fasta)
