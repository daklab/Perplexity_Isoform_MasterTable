#!/usr/bin/env python
# coding: utf-8

import csv 
import typer
from pathlib import Path
from typing import List

def read_sqanti(sqanti_output: Path) -> List[str]:
    isoform_set = []
    with open(sqanti_output, "rt") as sqanti:
        for line in sqanti:
            row = line.strip().split()
            isoform = row[0].split(".")[0]
            subcategory = row[14]

            if subcategory != "intron_retention":
                if isoform in isoform_set:
                    raise NameError("Something Wrong")
                isoform_set.append(isoform)
    
    return isoform_set


def filter_psl(input_psl: Path, isoform_set: List[str], output_psl: Path) -> None:
    with open(output_psl, "wt") as output_file:
        tsv_writer = csv.writer(output_file, delimiter="\t")
        with open(input_psl, "rt") as collapse:
            for line in collapse:
                row = line.strip().split("\t")
                isoform = row[9].split("_")[0]
                isoform = isoform.split(".")[0]
                if isoform in isoform_set:
                    tsv_writer.writerow(row)


def main(sqanti_output: Path, input_psl: Path, output_psl: Path) -> None:
    isoform_set = read_sqanti(sqanti_output)
    print(len(isoform_set))
    filter_psl(input_psl, isoform_set, output_psl)


if __name__ == "__main__":
    typer.run(main)
