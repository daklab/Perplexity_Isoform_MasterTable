#!/usr/bin/env python
# coding: utf-8


# add the longest UTR to each isoform (corrected read psl file, _gid psl file)

import argparse, csv, subprocess

parser = argparse.ArgumentParser(prog='add_utr.py', description='"Add the longest first and last block to the collapse output."')
parser.add_argument('map_file', action='store', type=str, metavar='map_file.txt', help='map file from collapse output')
parser.add_argument('read_psl', action='store', type=str, metavar='corrected_read.psl', help='corrected read psl file')
parser.add_argument('iso_psl', action='store', type=str, metavar='collapsed_isoform.psl', help='collapse output psl file')
parser.add_argument('-o', '--outfile', action='store', type=str, metavar='ourfile_name', help='output file name', required=True)
args = parser.parse_args()


def store_map(map_file: str) -> dict:
    with open(map_file, 'r') as map_txt:
        iso_to_map = {}
        for line in map_txt:
            isoform_map = line.strip().split("\t")
            iso_to_map[isoform_map[0]] = isoform_map[1].split(",")
    
    return iso_to_map


def replace_with_longest(n_th: int, iso_to_map: dict, read_psl: str, iso_psl: str) -> tuple:
    if n_th != 0 and n_th != -1:
        print("can only replace first or last block")
        raise ValueError
    
    with open(read_psl, 'r') as read_file:
        read_to_blocksizes = {}
        read_to_psl_line = {}
        for line in read_file:
            read = line.strip().split("\t")
            name = read[9]            
            blocksizes = read[18].split(",")[:-1]
            read_to_psl_line[name] = read
            read_to_blocksizes[name] = blocksizes
    
    iso_to_replace_line = {}
    for k, v in iso_to_map.items():
        mapped_to_nsize = {}
        max_size = 0
        for each in v:
            current = int(read_to_blocksizes[each][n_th])
            if max_size <= current:
                max_size = current
                max_read_line = read_to_psl_line[each]
            
        iso_to_replace_line[k] = max_read_line
    
    with open(iso_psl, 'r') as iso_file:
        plot_blocksizes = {}
        plot_matches = {}
        plot_qsize = {}
        plot_qend = {}
        plot_tstart = {}
        plot_tend = {}
        plot_qstarts = {}
        plot_tstarts = {}

        for line in iso_file:
            old_read = line.strip().split("\t")
            name = old_read[9]
            new_read = iso_to_replace_line[name]
            
            old_blocksizes = old_read[18].split(",")[:-1]
            new_blocksizes = new_read[18].split(",")[:-1]
            blocksizes = old_blocksizes
            n_old_blocksize = old_blocksizes[n_th]
            blocksizes[n_th] = new_blocksizes[n_th]
            plot_blocksizes[name] = [int(i) for i in blocksizes]
            
            length_diff = int(new_blocksizes[n_th]) - int(n_old_blocksize)
            old_matches = old_read[0]
            matches = int(old_matches) + length_diff
            plot_matches[name] = matches
            
            plot_qsize[name] = matches
            
            plot_qend[name] = matches
            
            old_tstart = old_read[15]
            if n_th == 0:
                tstart = int(old_tstart) - length_diff
            else:
                tstart = int(old_tstart)
            plot_tstart[name] = tstart
            
            old_tend = old_read[16]
            if n_th == 0:
                tend = int(old_tend)
            else:
                tend = int(old_tend) + length_diff
            plot_tend[name] = str(tend)
            
            qstarts = []
            point = 0
            for i in blocksizes:
                qstarts.append(point)
                point += int(i)
            plot_qstarts[name] = qstarts
            
            old_tstarts = old_read[20].split(",")[:-1]
            if n_th == 0:
                tstarts = [int(i) for i in old_tstarts]
                tstarts[0] = tstart
            else:
                tstarts = [int(i) for i in old_tstarts]
            plot_tstarts[name] = tstarts
            
    return plot_matches, plot_qsize, plot_qend, plot_tstart, plot_tend, plot_blocksizes, plot_qstarts, plot_tstarts   
            

def update_pslfile(outfile: str, plot_tup: tuple, psl_file: str):
    with open(outfile+"_isoform_full.psl", "wt") as isoform_psl:
        tsv_writer = csv.writer(isoform_psl, delimiter='\t')
        with open(psl_file, "r") as iso_psl:
            for line in iso_psl:
                read = line.split("\t")
                name = read[9] 

                read[0] = str(plot_tup[0][name])
                read[10] = str(plot_tup[1][name])
                read[12] = str(plot_tup[2][name])
                read[15] = str(plot_tup[3][name])
                read[16] = str(plot_tup[4][name])
                read[18] = str(plot_tup[5][name])[1:-1] + ","
                read[19] = str(plot_tup[6][name])[1:-1] + ","
                read[20] = str(plot_tup[7][name])[1:-1] + ","
                read[-1] = read[-1].strip()
                tsv_writer.writerow(read)
    print("The output file is isoform_full.psl file")                


iso_to_map = store_map(args.map_file)

read_psl = args.read_psl
iso_psl = args.iso_psl
outfile = args.outfile
plot_tup = replace_with_longest(0, iso_to_map, read_psl, iso_psl)
update_pslfile(outfile+"_first", plot_tup, iso_psl)



iso_psl_first = outfile+'_first_isoform_full.psl'
plot_tup = replace_with_longest(-1, iso_to_map, read_psl, iso_psl_first)
update_pslfile(outfile, plot_tup, iso_psl_first)
subprocess.call(['rm', iso_psl_first])


