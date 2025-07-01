#!/usr/bin/env python
# coding: utf-8

# SP Collapse Isoforms

# must have bedtools loaded to $PATH
import string, csv, subprocess, sys, argparse

parser = argparse.ArgumentParser(prog='sp_collapse_isoforms.py', description='"Collapse reads and give a representative set of isoforms."')
parser.add_argument('psl_file', action='store', type=str, metavar='corrected_reads.psl', help='corrected reads file to collapse')
parser.add_argument('--quantify', action='store_true', default=False)
parser.add_argument('-tsv', '--tsv_file', action='store', type=str, metavar='sample_info.tsv', help='tsv file for quantification containing sample id and information')
parser.add_argument('-tss', '--tss_bed', action='store', type=str, metavar='tss_region.bed', help='bed file for the region tss should overlap', required=True)
parser.add_argument('-tes', '--tes_bed', action='store', type=str, metavar='tes_region.bed', help='bed file for the region tes should overlap', required=True)
parser.add_argument('-id', '--pick_id', action='store_true', default=False)
parser.add_argument('-gtf', '--annotation_gtf', action='store', type=str, metavar='genome_annotation.gtf', help='annotation file for picking up transcript/gene id')
parser.add_argument('-o', '--outfile', action='store', type=str, metavar='outfile_name', help='outfile name without file extension', required=True)
parser.add_argument('-n', '--n_reads', action='store', type=int, metavar='collapse_reads_minimum', help='minimum number of reads required to collapse into one isoform', default=1)
args = parser.parse_args()



def get_key(my_dict, val):
    dict_key = list()
    for key, value in my_dict.items():
        if val == value:
            dict_key.append(key)
            
    if len(dict_key) > 1:
        return dict_key
    elif len(dict_key) == 1:
        return dict_key[0]
 
    return "none"


# read psl file and filter for promoter_supported
def promoter_filter(psl_file, promoter_bed, outfile):
    with open(psl_file, "r") as input_file:
        with open(outfile+"_tss.bed", "wt") as tss_bed:
            for line in input_file:
                read = line.split("\t")
                name = read[9]
                chromosome = read[13]
                strand = read[8]
                start = int(read[15])
                end = int(read[16])

                tsv_writer = csv.writer(tss_bed, delimiter='\t')
                if strand == "+":
                    tsv_writer.writerow([chromosome, start, start, name])
                elif strand == "-":
                    tsv_writer.writerow([chromosome, end, end, name])
                else:
                    tsv_writer.writerow([chromosome, start, start, name])
                    tsv_writer.writerow([chromosome, end, end, name])

    subprocess.check_call(
        ["bedtools", "intersect", "-a", outfile+"_tss.bed", "-b", promoter_bed],
        stdout=open(outfile+"_promoter_tss.bed", "wt")
    )
    
    with open(outfile+"_promoter_tss.bed", "r") as promoter_tss_bed:
        supported_names = set()
        for line in promoter_tss_bed:
            supported_name = line.split("\t")[3].strip()
            supported_names.add(supported_name)

    with open(psl_file, "r") as input_file:
        with open(outfile+"_promoter_supported.psl", "wt") as promoter_supported_psl:
            tsv_writer = csv.writer(promoter_supported_psl, delimiter='\t')
            for line in input_file:
                read = line.strip().split("\t")
                name = read[9]
                if name in supported_names:
                    tsv_writer.writerow(read)
                    
    print("Done filtering TSS region")


def tes_filter(psl_file, termination_bed, outfile):
    with open(psl_file, "r") as input_file:
        with open(outfile+"_tes.bed", "wt") as tes_bed:
            for line in input_file:
                read = line.split("\t")
                name = read[9]
                chromosome = read[13]
                strand = read[8]
                start = int(read[15])
                end = int(read[16])

                tsv_writer = csv.writer(tes_bed, delimiter='\t')
                if strand == "+":
                    tsv_writer.writerow([chromosome, end, end, name])
                elif strand == "-":
                    tsv_writer.writerow([chromosome, start, start, name])
                else:
                    tsv_writer.writerow([chromosome, start, start, name])
                    tsv_writer.writerow([chromosome, end, end, name])

    subprocess.check_call(["bedtools", "intersect", "-a", outfile+"_tes.bed", "-b", termination_bed], stdout=open(outfile+"_termination_tes.bed", "wt"))
    
    with open(outfile+"_termination_tes.bed", "r") as termination_tes_bed:
        supported_names = set()
        for line in termination_tes_bed:
            supported_name = line.split("\t")[3].strip()
            supported_names.add(supported_name)

    with open(psl_file, "r") as input_file:
        with open(outfile+"_tes_supported.psl", "wt") as tes_supported_psl:
            tsv_writer = csv.writer(tes_supported_psl, delimiter='\t')
            for line in input_file:
                read = line.strip().split("\t")
                name = read[9]
                if name in supported_names:
                    tsv_writer.writerow(read)
                    
    print("Done filtering TES region")                


# read psl file and collapse reads
def collapse_group_reads(pre_collapse_psl, outfile):
    with open(pre_collapse_psl, "r") as file:
        all_names = []
        all_tstarts = {}
        collapse_tstarts = dict()
        all_tends = dict()
        collapse_tends = dict()
        for line in file:
            read = line.split("\t")
            name = read[9]
            blocksizes = read[18][:-1].split(",")
            blocksizes = [int(i) for i in blocksizes]
            tstarts = read[20].strip()[:-1].split(",")
            tstarts = [int(i) for i in tstarts]
            tends = [x + y for (x, y) in zip(tstarts, blocksizes)]

            all_names.append(name)
            all_tstarts[name] = tstarts
            collapse_tstarts[name] = tstarts[1:]
            all_tends[name] = tends
            collapse_tends[name] = tends[:-1]
            
    list_tstarts = list()
    tuple_tstarts = list()
    id_tstarts = list()
    i = 0
    for k, v in collapse_tstarts.items():
        i += 1
        if v not in list_tstarts:
            list_tstarts.append(v)
            tuple_tstarts.append(tuple(v))
            unique_tstarts = 's'+str(hash(tuple(v)))
            id_tstarts.append(unique_tstarts) #hash id
    ref_tstarts = dict(zip(id_tstarts, list_tstarts))
    flipped_ref_tstarts = dict(zip(tuple_tstarts, id_tstarts))

    list_tends = list()
    tuple_tends = list()
    id_tends = list()
    for k, v in collapse_tends.items():
        i += 1
        if v not in list_tends:
            list_tends.append(v)
            tuple_tends.append(tuple(v))
            unique_tends = 'e'+str(hash(tuple(v)))
            id_tends.append(unique_tends) #hash id
    ref_tends = dict(zip(id_tends, list_tends))
    flipped_ref_tends = dict(zip(tuple_tends, id_tends))
    print("Extracted junction chains")

    flipped_tstarts = dict()
    read_start_id = dict()
    for k, v in collapse_tstarts.items():
        ids = flipped_ref_tstarts[tuple(v)]
        read_start_id[k] = ids
        if ids not in flipped_tstarts:
            flipped_tstarts[ids] = [k]
        else:
            flipped_tstarts[ids].append(k)     

    flipped_tends = dict()
    read_end_id = dict()
    for k, v in collapse_tends.items():
        ids = flipped_ref_tends[tuple(v)]
        read_end_id[k] = ids
        if ids not in flipped_tends:
            flipped_tends[ids] = [k]
        else:
            flipped_tends[ids].append(k)
    print("Junction chain grouped")

    read_positions = dict() # [read_starts, read_ends]
    for read in all_names:
        if read in read_start_id:
            read_starts = read_start_id[read]
        else:
            read_starts = "not found"
        
        if read in read_end_id:
            read_ends = read_end_id[read]
        else:
            read_ends = "not found"

        read_positions[read] = [read_starts, read_ends]

    flipped_positions = dict()
    for k, v in read_positions.items():
        position_id = ":".join(v)
        if position_id not in flipped_positions:
            flipped_positions[position_id] = [k]
        else:
            flipped_positions[position_id].append(k) 

    grouped_reads = {}
    unassigned_groups = {}
    for k, v in flipped_positions.items():
        if len(v) >= args.n_reads: # number of supporting reads
            grouped_reads[k] = v
        else:
            unassigned_groups[k] = v
    
    print("Collapsed the reads")
    
    return grouped_reads, unassigned_groups, ref_tstarts, ref_tends, flipped_ref_tstarts, flipped_ref_tends


# exclude single block isoforms
def exclude_single_block(grouped_reads_dict, ref_tstarts, ref_tends):
    multi_block_dict = dict()
    for key, value in grouped_reads_dict.items():
        if "ENST" in key:
            multi_block_dict[key] = value
        else:
            start_key, end_key = key.split(":")
            cds_start = ref_tstarts[start_key][0] if len(ref_tstarts[start_key]) > 0 else "single-block"
            cds_end = ref_tends[end_key][-1] if len(ref_tstarts[start_key]) > 0 else "single-block"
            if type(cds_start) is int and type(cds_end) is int:
                multi_block_dict[key] = value
    print("Excluded single block isoforms")
    return multi_block_dict
        


# remove UTRs and make first and last block length 1
# before picking up transcript id
def remove_utr(grouped_reads_multi_dict, ref_tstarts, ref_tends):
    plot_tstarts = dict()
    plot_tends = dict()
    plot_blocksizes = dict()
    plot_qstarts = dict()
    plot_qsize = dict()
    plot_matches = dict()
    plot_qend = dict()
    plot_tstart = dict()
    plot_tend = dict()
    for key in grouped_reads_multi_dict:
        start_key, end_key = key.split(":")
        plot_tstarts[key] = [ref_tends[end_key][0] - 1] + ref_tstarts[start_key]
        plot_tends[key] = ref_tends[end_key] + [ref_tstarts[start_key][-1] + 1] 
        plot_blocksizes[key] = [y - x for (x, y) in zip(plot_tstarts[key], plot_tends[key])]
        
        i = 0
        plot_qstarts[key] = list()
        for j in plot_blocksizes[key]:
            plot_qstarts[key].append(i)
            i += j
        
        plot_qsize[key] = plot_qstarts[key][-1] + plot_blocksizes[key][-1]
        plot_matches[key] = plot_qsize[key]
        plot_qend[key] = plot_qsize[key]
        plot_tstart[key] = plot_tstarts[key][0]
        plot_tend[key] = plot_tends[key][-1]

    print("Removed UTR and made the first and last block length 1")
    return plot_matches, plot_qsize, plot_qend, plot_tstart, plot_tend, plot_blocksizes, plot_qstarts, plot_tstarts
        


# pick up annotated transcipt id
def pick_transcript_id(gencode_transcripts, grouped_reads_dict, plot_noutr, ref_tstarts, ref_tends, flipped_ref_tstarts, flipped_ref_tends):
    with open(gencode_transcripts, "r") as annotated_file:
        annotated_tstarts = dict()
        annotated_tends = dict()
        for line in annotated_file:
            transcript = line.split("\t")
            name = transcript[9].split("_")[0]
            blocksizes = transcript[18][:-1].split(",")
            blocksizes = [int(i) for i in blocksizes]
            tstarts = transcript[20].strip()[:-1].split(",")
            tstarts = [int(i) for i in tstarts]
            tends = [x + y for (x, y) in zip(tstarts, blocksizes)]

            annotated_tstarts[name] = tstarts[1:]
            annotated_tends[name] = tends[:-1]

    transcript_positions = dict()
    for k, v in annotated_tstarts.items():
        if tuple(v) in flipped_ref_tstarts:
            annotated_tstarts[k] = flipped_ref_tstarts[tuple(v)]
        else: 
            annotated_tstarts[k] = "start_not_found"
        transcript_positions[k] = annotated_tstarts[k]
    for k, v in annotated_tends.items():
        if tuple(v) in flipped_ref_tends:
            annotated_tends[k] = flipped_ref_tends[tuple(v)]
        else: 
            annotated_tends[k] = "end_not_found"
        transcript_positions[k] += ":"+annotated_tends[k]
    
    flipped_transcript_positions = dict()
    for k, v in transcript_positions.items():
        if v not in flipped_transcript_positions:
            flipped_transcript_positions[v] = [k]
        else:
            flipped_transcript_positions[v].append(k)
    
    
    grouped_reads_tid = dict()    
    for k, v in grouped_reads_dict.items():
        if k in flipped_transcript_positions:
            if len(flipped_transcript_positions[k]) > 1:
                transcripts_group = "&".join(flipped_transcript_positions[k])
                grouped_reads_tid[transcripts_group] = v
            elif len(flipped_transcript_positions[k]) == 1:
                grouped_reads_tid[flipped_transcript_positions[k][0]] = v
        else:
            grouped_reads_tid[k] = v
    
    plot_noutr_tid = tuple()
    for feature in plot_noutr:
        feature_tid = dict()
        for k, v in feature.items():
            if k in flipped_transcript_positions:
                if len(flipped_transcript_positions[k]) > 1:
                    transcripts_group = "&".join(flipped_transcript_positions[k])
                    feature_tid[transcripts_group] = v
                elif len(flipped_transcript_positions[k]) == 1:
                    feature_tid[flipped_transcript_positions[k][0]] = v
            else:
                feature_tid[k] = v
        feature_list = list()
        feature_list.append(feature_tid)
        plot_noutr_tid += tuple(feature_list)
        
    print("Picked up annotated transcript ids for those who share the same juntion chain")
    return grouped_reads_tid, plot_noutr_tid


# read count file
def write_countfile(outfile, grouped_reads_dict):
    with open(outfile+"_counts.txt", "wt") as count_file:
        tsv_writer = csv.writer(count_file, delimiter='\t')
        for k, v in grouped_reads_dict.items():
            row = [k, len(v)]
            tsv_writer.writerow(row)
    print("Readcounts are logged inside the count.txt file")


# map file
def write_mapfile(outfile, grouped_reads_dict):
    with open(outfile+"_read_map.txt", "wt") as map_file:
        tsv_writer = csv.writer(map_file, delimiter='\t')
        for k, v in grouped_reads_dict.items():
            mapped_reads = ','.join(v)
            row = [k, mapped_reads]
            tsv_writer.writerow(row)
    print("Grouped reads are logged inside the map.txt file")


# write output psl file
def write_pslfile(outfile, grouped_reads_dict, plot_noutr, psl_file):
    representing_reads = dict()
    for k, v in grouped_reads_dict.items():
       # representing_reads[k] = v[0]
        representing_reads[v[0]] = k
    with open(outfile+"_isoform.psl", "wt") as isoform_psl:
        tsv_writer = csv.writer(isoform_psl, delimiter='\t')
        with open(psl_file, "r") as reads_psl:
            for line in reads_psl:
                read = line.split("\t")
                name = read[9]
                if name in representing_reads:
                    isoform = representing_reads[name] #get_key(representing_reads, name)
                else:
                    isoform = "none"

                if isoform != "none":
                    # just make the representing read to the isoform
                    read[0] = str(plot_noutr[0][isoform])
                    read[9] = isoform
                    read[10] = str(plot_noutr[1][isoform])
                    read[12] = str(plot_noutr[2][isoform])
                    read[15] = str(plot_noutr[3][isoform])
                    read[16] = str(plot_noutr[4][isoform])
                    read[18] = str(plot_noutr[5][isoform])[1:-1] + ","
                    read[19] = str(plot_noutr[6][isoform])[1:-1] + ","
                    read[20] = str(plot_noutr[7][isoform])[1:-1] + ","
                    read[-1] = read[-1].strip()
                    tsv_writer.writerow(read)
    print("The output file is isoform.psl file")                
                    
import gzip
import pandas as pd

def quantify_junctions(outfile, grouped_reads_dict, tsv_file):
    with open(tsv_file, "r") as tsv_in:
        read_sample = dict()
        df = dict()
        for line in tsv_in:
            grepped = list()
            sample = line.strip().split("\t")
            sample_id = ":".join([sample[0], sample[1]])
            sample_path = sample[2]
            if ".fastq" in sample_path or ".fq" in sample_path:
                if ".gz" in sample_path:
                    with gzip.open(sample_path, "rt") as sample_file:
                        for row in sample_file:
                            if "@" in row:
                                grepped.append(row)
                else:
                    sample_file = sample_path
                    for row in sample_file:
                        if "@" in row:
                            grepped.append(row)
            
            
            elif ".fasta" in sample_path or ".fa" in sample_path:
                if ".gz" in sample_path:
                    with gzip.open(sample_path, "rt") as sample_file:
                        for row in sample_file:
                            if ">" in row:
                                grepped.append(row)
                else:
                    sample_file = sample_path
                    for row in sample_file:
                        if ">" in row:
                            grepped.append(row)
            
            for read in grepped:
                if " " in read:
                    read = read.split(" ")[0]
                read = read[1:].strip()
                read_sample[read] = sample_id
            df[sample_id] = list()
    
    idx = list()
    for transcript, mapped_reads in grouped_reads_dict.items():
        idx.append(transcript)
        for sample in df:
            df[sample].append(0)
        for each in mapped_reads:
            each = each.split(";")[0]
            df[read_sample[each]][-1] += 1
            
    df_tsv=pd.DataFrame(data=df, index=idx)
    df_tsv.to_csv(outfile+"_quantify.tsv", sep="\t")

psl_file = args.psl_file
promoter_bed = args.tss_bed
termination_bed = args.tes_bed
gtf = args.annotation_gtf
outfile = args.outfile
tsv_file = args.tsv_file

# filter for promoter region and termination region
promoter_filter(psl_file, promoter_bed, outfile)
pre_collapse = outfile+"_promoter_supported.psl"
tes_filter(pre_collapse, termination_bed, outfile)
pre_collapse = outfile+"_tes_supported.psl"

# collapse reads and store the id info
grouped_reads, unassigned_groups, ref_tstarts, ref_tends, flipped_ref_tstarts, flipped_ref_tends = collapse_group_reads(pre_collapse, outfile)

# exclude single block reads
grouped_reads_multi = exclude_single_block(grouped_reads, ref_tstarts, ref_tends)
grouped_reads_multi_unassigend = exclude_single_block(unassigned_groups, ref_tstarts, ref_tends)

# remove UTRs
plot_noutr = remove_utr(grouped_reads_multi, ref_tstarts, ref_tends)
plot_noutr_unassigned = remove_utr(grouped_reads_multi_unassigend, ref_tstarts, ref_tends)

# pick up transcript id # should be fixed
if args.pick_id:
    subprocess.check_call([sys.executable, "/gpfs/commons/home/spark/sw/flair/src/flair/gtf_to_psl.py", gtf, outfile+"_annotated_transcript.psl", "--include_gene"])
    transcripts_psl = outfile+"_annotated_transcript.psl"

    grouped_reads_tid, plot_noutr_tid = pick_transcript_id(transcripts_psl, grouped_reads_multi, plot_noutr, ref_tstarts, ref_tends, flipped_ref_tstarts, flipped_ref_tends)
    
    grouped_reads = grouped_reads_tid
    plot_noutr = plot_noutr_tid
    
    unassigned_tid, plot_unassigned_tid = pick_transcript_id(transcripts_psl, unassigned_groups, plot_noutr_unassigned, ref_tstarts, ref_tends, flipped_ref_tstarts, flipped_ref_tends)
    
    unassigned_groups = unassigned_tid
    plot_noutr_unassigned = plot_unassigned_tid

# write output files
write_countfile(outfile, grouped_reads)
write_mapfile(outfile, grouped_reads)
#print("grouped_reads:", grouped_reads, file=sys.stdout)
#print("plot_noutr:", plot_noutr, file=sys.stdout)
write_pslfile(outfile, grouped_reads, plot_noutr, psl_file)

write_countfile(outfile+"_unassigned", unassigned_groups)
write_mapfile(outfile+"_unassigned", unassigned_groups)
write_pslfile(outfile+"_unassigned", unassigned_groups, plot_noutr_unassigned, psl_file)


# pick up gene-ID from FLAIR #should be fixed
if args.pick_id:
    subprocess.check_call([sys.executable, "/gpfs/commons/home/spark/sw/flair/src/flair/identify_gene_isoform.py", outfile+"_isoform.psl", gtf, outfile+"_isoform_gid.psl", "--gene_only"])
    subprocess.call(["rm", outfile+"_isoform.psl"])
    
    subprocess.check_call([sys.executable, "/gpfs/commons/home/spark/sw/flair/src/flair/identify_gene_isoform.py", outfile+"_unassigned_isoform.psl", gtf, outfile+"_unassigned_isoform_gid.psl", "--gene_only"])
    subprocess.call(["rm", outfile+"_unassigned_isoform.psl"])
    
    with open (outfile+"_isoform_gid.psl", "r") as psl:
        names_dict = dict()
        for line in psl:
            read = line.split("\t")
            name = read[9]
            if "_" in name:
                names_dict[name.split("_")[0]] = name
            else:
                raise ValueError
        
        grouped_reads_gid = dict()
        for k, v in grouped_reads.items():
            grouped_reads_gid[names_dict[k]] = v
    
    grouped_reads = grouped_reads_gid
    write_countfile(outfile+"_gid", grouped_reads_gid)
    subprocess.call(["rm", outfile+"_counts.txt"])
    write_mapfile(outfile+"_gid", grouped_reads_gid)
    subprocess.call(["rm", outfile+"_read_map.txt"])    
    
    print("Picked up gene id from annotated gtf file")
    
# quantification
if args.quantify:
    quantify_junctions(outfile, grouped_reads, args.tsv_file)
    print("Quantified the isoform abundance")
            
                 



