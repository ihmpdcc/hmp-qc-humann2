#!/usr/bin/env python3

'''
Description: This script is used to execute the HUMAnN2 docker pipeline 

Author: Kemi Ifeonu

Input: 
- An SRR ID or a set of fastq files 
- Mode: run qc only, humann2 only, or both
- S3_path if 
- AWS credentials

Output:
- Statistics file 
- QC'ed files
- HUMAnN2 Output files: _humann2_genefamilies.tsv, _humann2_pathabundance.tsv, _humann2_pathcoverage.tsv

'''

import argparse
import json
import pandas as pd
import os
import shutil
import logging
from subprocess import Popen, PIPE, STDOUT
import subprocess
import re
import sys
import boto3
from botocore.exceptions import ClientError
import hashlib
import bz2
import ntpath
import time

def main():
    parser = argparse.ArgumentParser( description='Execute Dockerized HUMAnN2 pipeline')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--srr', type=str, help='SRR ID')
    group.add_argument('-i', '--input', type=str, help='Prefix for a set of FASTQ files')
    parser.add_argument('-p', '--input_pair', type=str, required=False, help='paired-file')
    parser.add_argument('-m', '--mode', type=str, required=False, help='mode: qc, humann2, both',choices=['qc', 'humann2', 'both'])
    parser.add_argument('-b', '--bucket', type=str, required=False, help='Path to S3 bucket')
    

    global args 
    args = parser.parse_args()

    global f
    global l
    global srr
    global srr_2
    global json_string
    global result

    l = open("log.txt", 'w')
    json_string = []

    os.mkdir("input_seqs")

    # If SRR ID is provided, download data
    if args.srr is None:
        input_file = args.input
        if input_file.startswith('s3'):
            input_file_only = ntpath.basename(input_file)
            srr = os.path.splitext(input_file_only)[0]
            os.system("aws s3 cp " + input_file + " input_seqs/" + srr + ".fastq")
        else:
            srr = os.path.splitext(input_file)[0]
            os.system("cp input/" + input_file + " input_seqs/" + srr + ".fastq")
                        
        #os.system("ls -l input_seqs")
        summ_file = srr + "_summary_stats.txt"
        f = open(summ_file, 'w')

    else:
        srr = args.srr
        if srr.startswith('s3'):
            srr_full = srr
            srr = ntpath.basename(srr_full)
            srr = os.path.splitext(srr)[0]
            srr = srr.split('_')[0]
            os.system("aws s3 cp " + srr_full + " input_seqs/" + srr + "_1.fastq")
            os.system("aws s3 cp " + srr_full + " input_seqs/" + srr + "_2.fastq")
            summ_file = srr + "_summary_stats.txt"
            f = open(summ_file, 'w')
            f.write("Copied input files from S3 Bucket:\n")
        else:
            summ_file = srr + "_summary_stats.txt"
            f = open(summ_file, 'w')
            f.write("SRA Download:\n")
            print("\nSRA Download starting...\n")
            #download_cmd = "fasterq-dump " + srr + " --outdir input_seqs"
            #p = Popen(download_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            #result = p.stdout.read()
            #l.write(result.decode())
            fastq_download()
            if re.match('.*err.*',result.decode()):
                print("SRA download failed. Trying again in 1 minute...\n" + result.decode())
                time.sleep(60)
                fastq_download()
                if re.match('.*err.*',result.decode()):
                    print("SRA download failed. Trying again in 5 minutes...\n" + result.decode())
                    time.sleep(300)
                    fastq_download()
                    if re.match('.*err.*',result.decode()):
                        print("SRA download failed. Trying again in 15 minutes...\n" + result.decode())
                        time.sleep(900)
                        fastq_download()
                        if re.match('.*err.*',result.decode()):
                            print("SRA download failed after 4 tries. Exiting...\n" + result.decode())
                            sys.exit(1)

            if os.path.exists("input_seqs/"+ srr + ".sra_1.fastq"):
                os.system("mv input_seqs/" + srr + ".sra_1.fastq input_seqs/" + srr + "_1.fastq")
                os.system("mv input_seqs/" + srr + ".sra_2.fastq input_seqs/" + srr + "_2.fastq")
            elif os.path.exists("input_seqs/"+ srr + ".sra.fastq"):
                os.system("mv input_seqs/" + srr + ".sra.fastq input_seqs/" + srr + "_1.fastq")

            if os.path.isdir("input_seqs/sra"):
                os.system("rm -r input_seqs/sra")

            print("\nSRA Download complete...\n")
            sra_match = re.match('spots read\s+:\s+(\S+)\nreads read\s+:\s+(\S+)\nreads written\s+:\s+(\S+)',result.decode())
            spots_read = sra_match.group(1).replace(",","")
            reads_read = sra_match.group(2).replace(",","")
            reads_written = sra_match.group(3).replace(",","")
        
            hash_sra1 = hashlib.sha256(open("input_seqs/" + srr + "_1.fastq",'rb').read()).hexdigest()
            #hash_sra2 = hashlib.sha256(open("input_seqs/" + srr + "_2.fastq",'rb').read()).hexdigest()
            size_sra1 = os.path.getsize("input_seqs/" + srr + "_1.fastq")
            #size_sra2 = os.path.getsize("input_seqs/" + srr + "_2.fastq")
            #os.system("ls -l input_seqs")

            bz_sra = bz2.compress(open("input_seqs/" + srr + "_1.fastq", 'rb').read())
            sra1_bz = "input_seqs/" + srr + "_1.fastq.bz2"
            fh = open(sra1_bz, "wb")
            fh.write(bz_sra)
            fh.close()

            if os.path.exists("input_seqs/"+ srr + "_2.fastq"):
                hash_sra2 = hashlib.sha256(open("input_seqs/" + srr + "_2.fastq",'rb').read()).hexdigest()
                size_sra2 = os.path.getsize("input_seqs/" + srr + "_2.fastq")
                bz_sra2 = bz2.compress(open("input_seqs/" + srr + "_2.fastq", 'rb').read())
                sra2_bz = "input_seqs/" + srr + "_2.fastq.bz2"
                fh = open(sra2_bz, "wb")
                fh.write(bz_sra2)
                fh.close()
                size_sra2_bz = os.path.getsize("input_seqs/" + srr + "_2.fastq.bz2")
                file_2_name = srr + "_2.fastq.bz2"
            else:
                hash_sra2 = "n/a"
                size_sra2 = "n/a"
                bz_sra2 = "n/a"
                sra2_bz = "n/a"
                size_sra2_bz = "n/a"
                file_2_name = "n/a"

            size_sra1_bz = os.path.getsize("input_seqs/" + srr + "_1.fastq.bz2")
            #size_sra2_bz = os.path.getsize("input_seqs/" + srr + "_2.fastq.bz2")

            download_cmd = "fasterq-dump --version"
            p = Popen(download_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
            sra_version = p.stdout.read()
            sra_version = sra_version.decode().replace("\n","").replace("\"","")

            f.write("software: " + sra_version + "\nspots read: " + spots_read + "\nreads read: " + reads_read + \
                    "\nreads written: " + reads_written + \
                    "\n" + srr + "_1_file_name: " + srr + "_1.fastq.bz2" + \
                    "\n" + srr + "_1_file_size_uncompressed: " + str(size_sra1) + \
                    "\n" + srr + "_1_file_size_compressed: " + str(size_sra1_bz) + \
                    "\n" + srr + "_1_file_sha256: " + hash_sra1 +  \
                   # "\n" + srr + "_2_file_name: " + srr + "_2.fastq.bz2" + \
                    "\n" +  srr + "_2_file_name: " + file_2_name + \
                    "\n" + srr + "_2_file_size_uncompressed: " + str(size_sra2) + \
                    "\n" + srr + "_2_file_size_compressed: " + str(size_sra2_bz) + \
                    "\n" + srr + "_2_file_sha256: " + hash_sra2 + "\n")

            sra_download = {
                    'sra_download':{
                        'software': sra_version,
                        'spots_read': int(spots_read),
                        'reads_read': int(reads_read),
                        'reads_written': int(reads_written),
                        'file1':{
                            'file_name': srr + "_1.fastq",
                            'file_size': size_sra1,
                            'sha256': hash_sra1
                        },
                        'file2':{
                            'file_name': file_2_name,
                            'file_size': size_sra2,
                            'sha256': hash_sra2
                            }
                        }
                    }
            json_string.append(sra_download)

    if args.input_pair is not None:
        input_pair = args.input_pair
        if input_pair.startswith('s3'):
            input_pair_only = ntpath.basename(input_pair)
            srr_2 = os.path.splitext(input_pair_only)[0]
            os.system("aws s3 cp " + input_pair + " input_seqs/" + srr_2 + ".fastq")
        else:
            srr_2 = os.path.splitext(input_pair)[0]
            os.system("cp input/" + input_pair + " input_seqs/" + srr_2 + ".fastq")
                        
        #os.system("ls -l input_seqs")


    os.mkdir("final_output")

    if args.mode == 'qc':
        f.write("\nRunning QC only...\n")
        qc()
        os.system("rm final_output/" + srr +  "_qc.fastq")
    elif args.mode == 'humann2':
        f.write("\nRunning HUMAnN2 only...\n")
        if args.input_pair is None:
            os.system("cp input/" + args.input + " final_output/" + srr +  "_qc.fastq")
        else:
            os.system("cat input/" + args.input + " input/" + args.input_pair + " >  final_output/" + srr + "_qc.fastq")
        humann2()
        os.system("mm final_output/" + srr +  "_qc.fastq")
    else:
        f.write("\nRunning both QC and HUMAnN2...\n")
        qc()
        humann2()
        os.system("rm final_output/" + srr +  "_qc.fastq")
        
        
    f.close()
    os.system("mv " + summ_file + " final_output")
    json_summ = srr + "_summary_stats.json"
    with open("final_output/" + json_summ, 'w') as json_file:
        json.dump(json_string, json_file)
     
    #Upload output to S3 bucket
    print("\nUpload to S3 bucket " + args.bucket + " starting...\n")
    if args.bucket is not None:
        bucket = args.bucket
        bucket = bucket.replace(r's3://','')
        [s3_bucket,s3_path] = bucket.split('/', 1)
        s3_path = s3_path.rstrip('\/')
        s3_client = boto3.client('s3')

        for subdir, dirs, files in os.walk('final_output'):
            for file in files:
                full_path = os.path.join(subdir, file)
                dest_path = s3_path + "/" + file
                try:
                    response = s3_client.upload_file(full_path, s3_bucket, dest_path)
                except ClientError as e:
                    logging.error(e)
                    return False
                #return True
        print("\nUpload to S3 bucket " + args.bucket + " complete...\n")
    
    #clean up
    os.system("rm log.txt")
    os.system("rm out.txt")
    os.system("rm -r final_output")


def fastq_download():
    #download_cmd = "fasterq-dump " + srr + " --outdir input_seqs"
    #p = Popen(download_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    #global result
    #result = p.stdout.read()
    download_cmd = "prefetch " + srr + " -O input_seqs"
    p = Popen(download_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    #os.system("ls -l input_seqs")
    global result
    result = p.stdout.read()
    if re.match('.*err.*',result.decode()):
        return
    download_cmd = "fasterq-dump input_seqs/" + srr + "/" + srr + ".sra --outdir input_seqs"
    p = Popen(download_cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    result = p.stdout.read()

def qc():
    print("\nQC starting...\n")
    if args.input_pair is not None:
        os.system("kneaddata --input input_seqs/" + srr + ".fastq --input input_seqs/" + srr_2 + ".fastq -db /dbs/kneaddata --output kneaddata_output --cat-final-output --threads 6  2>&2 | tee -a out.txt")
    elif args.input is not None:
        os.system("kneaddata --input input_seqs/" + srr + ".fastq -db /dbs/kneaddata --output kneaddata_output --cat-final-output --threads 6  2>&2 | tee -a out.txt")
    elif os.path.exists("input_seqs/"+ srr + ".sra_2.fastq"):
        os.system("kneaddata --input input_seqs/" + srr + "_1.fastq --input input_seqs/" + srr + "_2.fastq -db /dbs/kneaddata --output kneaddata_output --cat-final-output --threads 6  2>&2 | tee -a out.txt")
    else:
        os.system("kneaddata --input input_seqs/" + srr + "_1.fastq -db /dbs/kneaddata --output kneaddata_output --cat-final-output --threads 6  2>&2 | tee -a out.txt")


    qc_file = srr + "_qc.fastq"

    if not os.path.exists(srr + "_qc.fastq"):
        print("\nQC failed\n")
        sys.exit(1)

    os.system("mv kneaddata_output/*_kneaddata.log kneaddata_output/kneaddata.log")
    initial_read_count = 0
    read_count_after_trimming = 0
    read_count_after_decontamination = 0
    log = open("kneaddata_output/kneaddata.log", "r")
    for line in log:
        if re.match('.*ERROR.*',line):
            print("QC failed\n" + line)
            sys.exit(1)

        m0 = re.match(".*Initial number of reads.*:\s(\w+)", line)
        m1 = re.match(".*Total reads after trimming.*:\s(\w+)", line)
        m2 = re.match(".*Total reads after merging results from multiple databases.*:\s(\w+)", line)
        if m0:
            initial_read_count = initial_read_count + int(m0.group(1))
        if m1:
            read_count_after_trimming = read_count_after_trimming + int(m1.group(1))
        if m2:
            read_count_after_decontamination = read_count_after_decontamination + int(m2.group(1))
    
    #delete orig files
    os.system("rm -r input_seqs")

    os.system("mv kneaddata_output/" + srr + "*_kneaddata.fastq final_output/" + qc_file)
    bz = bz2.compress(open("final_output/" + qc_file, 'rb').read())
    qc_file_bz = "final_output/" + qc_file + ".bz2"
    fh = open(qc_file_bz, "wb")
    fh.write(bz)
    fh.close()
    hash_qc = hashlib.sha256(open("final_output/" + qc_file + ".bz2",'rb').read()).hexdigest()
    size_qc = os.path.getsize("final_output/" + qc_file)
    size_bz_qc = os.path.getsize("final_output/" + qc_file + ".bz2")
    #os.system("rm final_output/" + qc_file)

    cmd = "kneaddata --version |awk '{ print $2 }'"
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    kd_version = p.stdout.read()
    kd_version = kd_version.decode().replace("\n","")#.replace("\"","")

    cmd = "bowtie2 --version |awk '{print $3}'"
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    bt_version = p.stdout.read()
    bt_version = bt_version.decode().replace("\n","")#.replace("\"","")

    f.write("QC:\nsoftware: \nkneaddata: " + kd_version + "\nbowtie: " + bt_version + "\ntrimmomatic: v0.33 \ninitial read count: " + str(initial_read_count) + "\nread count after trimming: " + str(read_count_after_trimming) + "\nread count after decontamination: " + str(read_count_after_decontamination) + "\nQC Output file:\n"+ srr+ "_QC_file: " + qc_file + ".bz2\n" + srr + "_QC_file_size_compressed: " + str(size_bz_qc) + "\n" + srr + "_QC_file_size_uncompressed: " + str(size_qc)+ "\n" + srr + "_QC_file_sha256_compressed: " + hash_qc )
    qc_string = {
            'qc':{
                'software':{
                    'kneaddata': kd_version,
                    'bowtie': bt_version,
                    'trimmomatic': 'v0.33'
                    },
                'initial_read_count': initial_read_count,
                'read_count_after_trimming': read_count_after_trimming,
                'read_count_after_decontamination': read_count_after_decontamination,
                'qc_output_file': {
                    "file_name": qc_file + '.bz2' ,
                    "file_size_compressed": size_bz_qc,
                    "file_size_uncompressed": size_qc,
                    "sha256": hash_qc
                    }
                }
            }
    json_string.append(qc_string)
    
    #delete extra kneaddata files
    os.system("rm -r kneaddata_output")
    print("\nQC complete...\n")


def humann2():
    print("\nHUMAnN2 starting...\n")
    return_code = os.system("humann2 --input final_output/" + srr + "_qc.fastq --output humann2_output --threads 6 --metaphlan-options=\"--mpa_pkl /dbs/humann2/metaphlan/mpa_v20_m200.pkl --bowtie2db /dbs/humann2/metaphlan/mpa_v20_m200\"")
    if return_code != 0:
        sys.exit("HUMAnN2 failed" + return_code)

    os.system("mv humann2_output/*_genefamilies.tsv final_output/" + srr + "_humann2_genefamilies.tsv")
    os.system("mv humann2_output/*_pathcoverage.tsv final_output/" + srr + "_humann2_pathcoverage.tsv")
    os.system("mv humann2_output/*_pathabundance.tsv final_output/" + srr + "_humann2_pathabundance.tsv")

    os.system("mv humann2_output/*_humann2_temp/*.log log.txt")
    
    #remove extra ouput files
    os.system("rm -r humann2_output")

    nuc_genes = 0
    nuc_unalign = 0
    trans_genes = 0
    trans_unalign = 0

    log = open("log.txt", "r")
    for line in log:
        h0 = re.match(".*Total gene families from nucleotide alignment.*:\s(\w+)", line)
        h1 = re.match(".*Total gene families from nucleotide alignment.*:\s(\w+)", line)
        h2 = re.match(".*Total gene families after translated alignment.*:\s(\w+)", line)
        h3 = re.match(".*Unaligned reads after translated alignment.*:\s(\w+)", line)
        if h0:
            nuc_genes = int(h0.group(1))
        if h1:
            nuc_unalign = int(h1.group(1))
        if h2:
            trans_genes = int(h2.group(1))
        if h3:
            trans_unalign = int(h3.group(1))



    
    cmd = "humann2 --version |awk '{ print $2 }'"
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    h2_version = p.stdout.read()
    h2_version = h2_version.decode().replace("\n","")#.replace("\"","")

    cmd = "metaphlan2.py --version |awk '{ print $3 }'"
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    m2_version = p.stdout.read()
    m2_version = m2_version.decode().replace("\n","")#.replace("\"","")

    cmd = "diamond --version |awk '{ print $3 }'"
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    dia_version = p.stdout.read()
    dia_version = dia_version.decode().replace("\n","")#.replace("\"","")



    hash_gf = hashlib.sha256(open("final_output/" + srr + "_humann2_genefamilies.tsv",'rb').read()).hexdigest()
    hash_pc = hashlib.sha256(open("final_output/" + srr + "_humann2_pathcoverage.tsv",'rb').read()).hexdigest()
    hash_pa = hashlib.sha256(open("final_output/" + srr + "_humann2_pathabundance.tsv",'rb').read()).hexdigest()

    bz_gf = bz2.compress(open("final_output/" + srr + "_humann2_genefamilies.tsv", 'rb').read())
    humann_gf_bz = "final_output/" + srr + "_humann2_genefamilies.tsv.bz2"
    fh = open(humann_gf_bz, "wb")
    fh.write(bz_gf)
    fh.close()

    bz_pc = bz2.compress(open("final_output/" + srr + "_humann2_pathcoverage.tsv", 'rb').read())
    humann_pc_bz = "final_output/" + srr + "_humann2_pathcoverage.tsv.bz2"
    fh = open(humann_pc_bz, "wb")
    fh.write(bz_pc)
    fh.close()
    
    bz_pa = bz2.compress(open("final_output/" + srr + "_humann2_pathabundance.tsv", 'rb').read())
    humann_pa_bz = "final_output/" + srr + "_humann2_pathabundance.tsv.bz2"
    fh = open(humann_pa_bz, "wb")
    fh.write(bz_pa)
    fh.close()

    size_gf = os.path.getsize("final_output/" + srr + "_humann2_genefamilies.tsv")
    size_pc = os.path.getsize("final_output/" + srr + "_humann2_pathcoverage.tsv")
    size_pa = os.path.getsize("final_output/" + srr + "_humann2_pathabundance.tsv")

    size_gf_bz = os.path.getsize("final_output/" + srr + "_humann2_genefamilies.tsv.bz2")
    size_pc_bz = os.path.getsize("final_output/" + srr + "_humann2_pathcoverage.tsv.bz2")
    size_pa_bz = os.path.getsize("final_output/" + srr + "_humann2_pathabundance.tsv.bz2")

    os.system("rm -r final_output/*.tsv")

    f.write("\n\nHUMAnN2:\nsoftware: \nhumann2: " + h2_version + "\nmetaphlan2: " + m2_version + "\ndiamond: " + dia_version + \
            "\nTotal gene families after translated alignment " + str(nuc_genes) + \
            "\nUnaligned reads after translated alignment: " + str(nuc_unalign) + \
            "%\nTotal gene families after translated alignment: " + str(trans_genes) + \
            "\nUnaligned reads after translated alignment: " + str(trans_unalign) + \
            "%\n" + srr + "_genefamilies_file_name: " +  srr + "_humann2_genefamilies.tsv.bz2" + \
            "\n" + srr + "_genefamilies_file_size_uncompressed: " + str(size_gf) + \
            "\n" + srr + "_genefamilies_file_size_compressed: " + str(size_gf_bz) + \
            "\n" + srr + "_genefamilies_sha256: " + hash_gf + \
            "\n" + srr + "_pathabundance_file_name: " +  srr + "_humann2_pathabundance.tsv.bz2" + \
            "\n" + srr + "_pathabundance_file_size_uncompressed: " + str(size_pa) + \
            "\n" + srr + "_pathabundance_file_size_compressed: " + str(size_pa_bz) + \
            "\n" + srr + "_pathabundance_sha256: " + hash_pa + \
            "\n" + srr + "_pathcoverage_file_name: " +  srr + "_humann2_pathcoverage.tsv.bz2" + \
            "\n" + srr + "_pathcoverage_file_size_uncompressed: " + str(size_pc) + \
            "\n" + srr + "_pathcoverage_file_size_compressed: " + str(size_pc_bz) + \
            "\n" + srr + "_pathcoverage_sha256: " + hash_pc)
    
    humann_string = {
            'humann2':{
                'software':{
                    'humann2': h2_version,
                    'metaphlan2': m2_version,
                    'diamond': dia_version
                    },
                'gene_families_per_nuc_align': nuc_genes,
                'unaligned_after_nuc_align': str(nuc_unalign) + "%",
                'gene_families_per_translated_align': trans_genes,
                'unaligned_after_translated_align': str(trans_unalign) + "%",
                'humann2_genefamilies':{
                    'file_name': srr + "_humann2_genefamilies.tsv.bz2",
                    'file_size': size_gf_bz,
                    'file_size_uncompressed': size_gf,
                    'sha256': hash_gf
                    },
                'humann2_pathabundance':{
                    'file_name': srr + "_humann2_pathabundance.tsv.bz2",
                    'file_size': size_pa_bz,
                    'file_size_uncompressed': size_pa,
                    'sha256': hash_pa
                    },
                'humann2_pathcoverage':{
                    'file_name': srr + "_humann2_pathcoverage.tsv.bz2",
                    'file_size': size_pc_bz,
                    'file_size_uncompressed': size_pc,
                    'sha256': hash_pc
                    }
                }
            }
    json_string.append(humann_string)
    print("\HUMAnN2 complete...\n")


if __name__ == '__main__':
    main()
    

