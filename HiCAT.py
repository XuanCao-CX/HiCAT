#!/usr/bin/env python
# -*- coding:utf-8 -*-

import gzip
import os
import sys
import getopt
from joblib import Parallel, delayed
from datetime import datetime

def qname(fq_path, output):
    if not fq_path.endswith(".gz"):
        fq = open(fq_path)
    else:
        fq = gzip.open(fq_path, 'rb')
    output = open(output, 'w')
    num = 1
    line = fq.readline().rstrip('\n')
    while line:
        if num % 4 == 1:
            name = "@" + str(num/4)
            output.write(name + '\n')
        elif num % 4 == 2:
            output.write(line + '\n')
        elif num % 4 == 3:
            output.write("+\n")
        elif num % 4 == 0:
            output.write(line + '\n')
        num += 1
        line = fq.readline().rstrip('\n')
    fq.close()
    output.close()

def full_mapping(fq_path,bowtie_index, cpu, mismatch, trim5, trim3, unmapped_fq_path, sam_path, bowtie_log_path):
    bowtie_command = 'bowtie %s -p %d -m 1 -n %d -5 %d -3 %d --un %s -q %s -S %s 2> %s'\
                     % (bowtie_index, cpu, mismatch, trim5, trim3, unmapped_fq_path, fq_path, sam_path, bowtie_log_path)
    os.system(bowtie_command)

def RE_mapping(RE_fq_path, bowtie_index, cpu,mismatch, sam_path, bowtie_log_path):
    bowtie_command = 'bowtie %s -p %d -m 1 -n %d -q %s -S %s 2> %s' \
                         % (bowtie_index, cpu, mismatch, RE_fq_path,sam_path, bowtie_log_path)
    os.system(bowtie_command)

def REsite_make(sites5):
    sites5 = [x.upper() for x in sites5]
    rule = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    sites5_reverse = ["".join(map(lambda x: rule[x], y))[::-1] for y in sites5]
    resites = [x + '+' + y for x in sites5 for y in sites5_reverse]
    resites = list(set(resites))
    resites=sorted(resites)
    resites.reverse()
    return resites

def REsite(un_fq_path, sites5, RE_fq_path, trim5, trim3, REsite_log_path):
    resites = REsite_make(sites5)
    resites_num = {x: [0, 0, 0] for x in resites} # REsite: [rds_has_REsite_num, cutted_left_rds_length_than_20bp,cutted_right_rds_length_than_20bp]
    RE_fq = open(RE_fq_path, 'w')
    fq = open(un_fq_path)
    line = fq.readline().rstrip('\n')
    while line:
        qname = line
        seq = fq.readline().rstrip('\n')
        line3 = fq.readline().rstrip('\n')
        quantily = fq.readline().rstrip('\n')
        if trim3 == 0:
            seq = seq[trim5:]
            quantily = quantily[trim5:]
        else:
            seq = seq[trim5:-trim3]
            quantily = quantily[trim5:-trim3]
        for resite in resites:
            s1, s2 = resite.split('+')
            s = s1 + s2
            l1 = len(s1)
            resite_pos = seq.find(s)
            if resite_pos == -1:
                continue
            else:
                resites_num[resite][0] += 1
                cut_site = resite_pos + l1
                left_seq = seq[:cut_site]
                left_quantily = quantily[:cut_site]
                right_seq = seq[cut_site:]
                right_quantily = quantily[cut_site:]
                right_len = len(right_quantily)
                if cut_site >= 20:
                    resites_num[resite][1] += 1
                    RE_fq.write(qname + '\n')
                    RE_fq.write(left_seq + '\n')
                    RE_fq.write(line3 + '\n')
                    RE_fq.write(left_quantily + '\n')
                if right_len >= 20:
                    resites_num[resite][2] += 1
                    RE_fq.write(qname + '\n')
                    RE_fq.write(right_seq + '\n')
                    RE_fq.write(line3 + '\n')
                    RE_fq.write(right_quantily + '\n')
                break
        line = fq.readline().rstrip('\n')
    fq.close()
    RE_fq.close()
    RE_summary_log = open(REsite_log_path, 'w')
    RE_summary_log.write("REsite\treads_have_RE\t5_reads_length_than_20bp\t3_reads_length_than_20bp\n")
    for resite in resites:
        resite_num = resites_num[resite]
        RE_summary_log.write('%s\t%d\t%d\t%d\n' % (resite,resite_num[0], resite_num[1], resite_num[2]))
    RE_summary_log.close()

def sam_split(sam_file, full_RE, R1_R2 ,pair_dir, rds_num):
    split_command = 'awk \'BEGIN{FS=OFS="\\t"}{if ($2 ==0 || $2 == 16){idx=int($1/%d);print $1,$3,$4,"%s","%s" >> "%s/hits."idx}}\' %s' \
                    % (rds_num, full_RE, R1_R2, pair_dir, sam_file)
    os.system(split_command)

def pairing(hit_file, genomic_dis, pair_dir):
    pair_name = os.path.split(hit_file)[1]
    pair_name = pair_name.replace("hits", "paired")
    intra_file = pair_dir + os.sep + pair_name + '.intra.pets'
    inter_file = pair_dir + os.sep + pair_name + '.inter.pets'
    paired_log = pair_dir + os.sep + pair_name + '.log'
    intra_file = open(intra_file, 'w')
    inter_file = open(inter_file, 'w')
    hit_file = open(hit_file)
    hits = {}
    ### read hits_file
    for line in hit_file:
        line = line.rstrip('\n').split('\t')
        qname, chr, pos, full_RE, R1_R2 = int(line[0]), line[1], int(line[2]), line[3], line[4]
        if qname in hits:
            hits[qname].append([chr, pos, full_RE, R1_R2])
        else:
            hits[qname] = [[chr, pos, full_RE, R1_R2]]
    hit_file.close()
    paired_type_names = ["intra_pets", "inter_pets", "genomic",
                         "R1full_R2full_intra", "R1full_R2RE_intra", "R1RE_R2full_intra",
                          "R1RE_R2RE_intra", 'R1RE_R1RE_intra', 'R2RE_R2RE_intra',
                          "R1full_R2full_inter", "R1full_R2RE_inter", "R1RE_R2full_inter",
                          "R1RE_R2RE_inter", 'R1RE_R1RE_inter', 'R2RE_R2RE_inter',
                          "R1full_unpaired", "R1RE_unpaired", "R2full_unpaired", "R2RE_unpaired"]\
                        + ["reads_has_%d_pets" % x for x in range(8)]
    paired_type_number = {x: 0 for x in paired_type_names}
    for qname in hits:
        hits_one_qname = hits[qname]
        hits_one_qname = sorted(hits_one_qname, key=lambda  x: (x[0], x[1]))
        if len(hits_one_qname) == 1:
            paired_type = hits_one_qname[0][3] + hits_one_qname[0][2] + '_unpaired'
            paired_type_number[paired_type] += 1
        else:
            ### remove genomic DNA
            rm_genomic_hits = [hits_one_qname[0]]
            for j in range(1, len(hits_one_qname)):
                h1 = rm_genomic_hits[-1]
                h2 = hits_one_qname[j]
                if h1[0] == h2[0]:
                    dis = h2[1] - h1[1]
                    if dis <= genomic_dis:
                        rm_genomic_hits[-1] = h2
                        paired_type_number['genomic'] += 1
                    else:
                        rm_genomic_hits.append(h2)
                else:
                    rm_genomic_hits.append(h2)
            ### paired pets
            rm_genomic_hits_num = len(rm_genomic_hits)
            for x in range(rm_genomic_hits_num - 1):
                chr_A, pos_A, full_RE_A ,R1_R2_A = rm_genomic_hits[x]
                for y in range(x + 1, rm_genomic_hits_num):
                    chr_B, pos_B, full_RE_B, R1_R2_B = rm_genomic_hits[y]
                    if chr_A == chr_B:
                        paired_type_number["intra_pets"] += 1
                        intra_file.write('%s\t%d\t%s\t%d\n' %(chr_A, pos_A, chr_B, pos_B))
                        type = "intra"
                    else:
                        paired_type_number["inter_pets"] += 1
                        if chr_A > chr_B:
                            chr_A, chr_B = chr_B, chr_A
                            pos_A, pos_B = pos_B, pos_A
                        inter_file.write('%s\t%d\t%s\t%d\n' % (chr_A, pos_A, chr_B, pos_B))
                        type = "inter"
                    #### pets_paired _type
                    paired_type_A = R1_R2_A + full_RE_A
                    paired_type_B = R1_R2_B + full_RE_B
                    if paired_type_A > paired_type_B:
                        paired_type_A,paired_type_B = paired_type_B, paired_type_A
                    paired_type = paired_type_A + '_' + paired_type_B + '_' + type
                    paired_type_number[paired_type] += 1
            ### calculate one qname has the paired pets number
            one_qname_paired_num = rm_genomic_hits_num * (rm_genomic_hits_num - 1) / 2
            paired_type_number["reads_has_%d_pets" % one_qname_paired_num] += 1
    intra_file.close()
    inter_file.close()
    paired_log = open(paired_log, 'w')
    for key in paired_type_names:
        paired_log.write("%s\t%d\n" % (key, paired_type_number[key]))
    paired_log.close()
    return paired_type_number

def split_by_chr(paired_pets_path, pets_dir):
    split_command = 'awk \'BEGIN{FS=OFS="\\t"}{print $0 >> "%s/"$1".vs."$3".pets.no_rmdup"}\' %s' \
                    % (pets_dir, paired_pets_path)
    os.system(split_command)

def pets_remove_dumplicate(pets_file):
    f = open(pets_file)
    f_out = pets_file.replace(".no_rmdup", "")
    f_out = open(f_out, 'w')
    keep_num = 0
    pets = []
    for line in f:
        line = line.rstrip('\n').split('\t')
        chr1, pos1, chr2, pos2 = line[0], int(line[1]), line[2], int(line[3])
        pets.append([chr1, pos1, chr2, pos2])
    f.close()
    pets = sorted(pets, key=lambda x: (x[1], x[3]))
    pet_A = pets[0]
    for i in xrange(1, len(pets)):
        pet_B = pets[i]
        if pet_A[1] == pet_B[1] and pet_A[3] == pet_B[3]:
            pet_A = pet_B
        else:
            keep_num += 1
            f_out.write('%s\t%d\t%s\t%d\n' % (pet_A[0], pet_A[1], pet_A[2], pet_A[3]))
            pet_A = pet_B
    keep_num += 1
    f_out.write('%s\t%d\t%s\t%d\n' % (pet_A[0], pet_A[1], pet_A[2], pet_A[3]))
    f_out.close()
    return keep_num

def usage():
    print 'Uasge: HiC pipline.'
    print '     -1              The input R1 fastq or fastq.gz files. Files seperated by comma.'
    print '     -2              The input R2 fastq or fastq.gz files. Files seperated by comma.'
    print '                     The order of files can paired with R1 files.'
    print '     -o              The output dir.'
    print '     --label              The lable for each paired R1_R2 fastq files, seperated by comma.'
    print '     -g              The bowtie index.'
    print '     -p              Cpus. Default:1.'
    print '     -n              The mismatch number for bowtie mapping. Default=2.'
    print '     -5              Trim bases number of 5 prime.'
    print '     -3              Trim bases number of 3 prime.'
    print '     --re            The resites for cut unmapped fastq files. It is the cut site from 5\'-> 3\'.'
    print '                     Default, DPNII:GATC.'
    print '     --genomic_dis   The length for remove genomic DNA.Default:500.'
    print '     --start_search_REsite       The starting step is search REsite.If set, you should set old_file_dir.'
    print '     --start_REsite_mapping      The starting step is search RE reads mapping.If set, you should set old_file_dir.'
    print '     --old_file_dir              If not start at the full length map step,set this to get former files.'
    print '     --keep_media_fq_files        If set, keep the median sam and fastq files. Default delete them.'

if __name__ == "__main__":
    R1_fq_files = ""
    R2_fq_files = ""
    output_dir = ""
    labels = []
    bowtie_index = ""
    cpu = 1
    mismatch = 2
    trim5 = 0
    trim3 = 0
    sites5 = ["GATC"]
    split_cutoff = 1000000
    genomic_dis = 500
    re = False
    nomap = False
    pets_cutoff = 4000000
    start_search_REsite = False
    start_REsite_mapping = False
    old_file_dir = False
    keep_media_fq_files = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h1:2:o:g:p:n:5:3:',
                                   ['label=','re=', 'genomic_dis=','split_cutoff=','start_search_REsite',
                                    'start_REsite_mapping','old_file_dir=','keep_media_fq_files'])
        for o, a in opts:
            if o == '-h':
                usage()
                sys.exit()
            if o == '-1':
                R1_fq_files = a.split(',')
            if o == '-2':
                R2_fq_files = a.split(',')
            if o == '-o':
                output_dir = a
            if o == '--label':
                labels = a.split(',')
            if o == '-g':
                bowtie_index = a
            if o == '-p':
                cpu = int(a)
            if o == '-n':
                mismatch = int(a)
            if o == '-5':
                trim5 = int(a)
            if o == '-3':
                trim3 = int(a)
            if o == '--re':
                sites5 = a.split(",")
            if o == '--split_cutoff':
                split_cutoff = int(a)
            if o == '--genomic_dis':
                genomic_dis = int(a)
            if o == '--old_file_dir':
                old_file_dir = a
            if o == '--start_search_REsite':
                start_search_REsite = True
            if o == '--start_REsite_mapping':
                start_REsite_mapping = True
            if o == '--keep_media_fq_files':
                keep_media_fq_files = True
    except getopt.GetoptError:
        print 'Error in getting parametres!'
        sys.exit()

    ### check parameter
    if R1_fq_files != "":
        if len(labels) != len(R1_fq_files):
            print "Fastq files number non-equal labels number!"
            sys.exit()

    command = 'hic_xuan.py ' + '-1 ' + ','.join(R1_fq_files) + ' -2 ' + ','.join(R2_fq_files) + ' -o ' + output_dir\
              + ' -g ' + bowtie_index + ' -p ' + str(cpu) + ' -n ' + str(mismatch) + ' -5 ' + str(trim5) +\
              ' -3 ' + str(trim3) + ' --re ' + ','.join(sites5) + ' --genome_dir ' + str(genomic_dis) +\
              ' --old_file_dir ' + str(old_file_dir) + ' --start_search_REsite ' + str(start_search_REsite) +\
              ' --start_REsite_mapping ' + str(start_REsite_mapping) + ' --keep_median_files ' + str(keep_media_fq_files)
    if start_search_REsite == False and start_REsite_mapping == False:
        start_full_mapping = True
        start_search_REsite = True
        start_REsite_mapping = True
        full_map_dir = output_dir + os.sep + 'map'
        RE_searched_fq_dir = output_dir + os.sep + 'map'
    elif start_search_REsite:
        start_full_mapping = False
        start_search_REsite = True
        start_REsite_mapping = True
        full_map_dir = old_file_dir + os.sep + 'map'
        RE_searched_fq_dir = output_dir + os.sep + 'map'
    elif start_REsite_mapping:
        start_full_mapping = False
        start_search_REsite = False
        start_REsite_mapping = True
        full_map_dir = old_file_dir + os.sep + 'map'
        RE_searched_fq_dir = old_file_dir + os.sep + 'map'

    ### make dirs
    RE_map_dir = output_dir + os.sep + 'map'
    pair_dir = output_dir + os.sep + 'pair'
    pair_sub_dirs = [pair_dir + os.sep + label for label in labels]
    pets_dir = output_dir + os.sep + 'pets'
    pets_dir = output_dir + os.sep + 'pets'
    pets_intra_dir = pets_dir + os.sep + 'intra'
    pets_inter_dir = pets_dir + os.sep + 'inter'
    if os.path.exists(output_dir):
        print "Output dir \"%s\" is exist!!! Please make sure the output dir is new maked!" % output_dir
        sys.exit()
    else:
        os.makedirs(output_dir)
    dirs = [full_map_dir, RE_searched_fq_dir, RE_map_dir,pair_dir] + pair_sub_dirs + [pets_intra_dir,pets_inter_dir]
    for one_dir in dirs:
        if not os.path.exists(one_dir):
            os.makedirs(one_dir)

    ### making file names
    hic_log_file = output_dir + os.sep + 'hic_log.txt'
    os.system('echo %s > %s ' % (command, hic_log_file))

    R1_fq_files_qname_changed = [full_map_dir + os.sep + label + '.R1.qname_changed.fq' for label in labels]
    R2_fq_files_qname_changed = [full_map_dir + os.sep + label + '.R2.qname_changed.fq' for label in labels]

    R1_full_sam_pathes = [full_map_dir + os.sep + label + '.R1.full.sam' for label in labels]
    R1_full_unmap_fq = [full_map_dir + os.sep + label + '.R1.full.un.fq' for label in labels]
    R1_full_bowtie_log = [full_map_dir + os.sep + 'bowtie.full.R1.' + label + '.txt' for label in labels]
    R2_full_sam_pathes = [full_map_dir + os.sep + label + '.R2.full.sam' for label in labels]
    R2_full_unmap_fq = [full_map_dir + os.sep + label + '.R2.full.un.fq' for label in labels]
    R2_full_bowtie_log = [full_map_dir + os.sep + 'bowtie.full.R2.' + label + '.txt' for label in labels]

    R1_RE_fq_pathes = [RE_searched_fq_dir + os.sep + label + '.R1.RE.fq' for label in labels]
    R2_RE_fq_pathes = [RE_searched_fq_dir + os.sep + label + '.R2.RE.fq' for label in labels]
    R1_RE_summary_log = [RE_searched_fq_dir + os.sep + 'REsite.R1.' + label + '.txt' for label in labels]
    R2_RE_summary_log = [RE_searched_fq_dir + os.sep + 'REsite.R2.' + label + '.txt' for label in labels]

    R1_RE_sam_pathes = [RE_map_dir + os.sep + label + '.R1.RE.sam' for label in labels]
    R1_RE_bowtie_log = [RE_map_dir + os.sep + 'bowtie.RE.R1.' + label + '.txt' for label in labels]
    R2_RE_sam_pathes = [RE_map_dir + os.sep + label + '.R2.RE.sam' for label in labels]
    R2_RE_bowtie_log = [RE_map_dir + os.sep + 'bowtie.RE.R2.' + label + '.txt' for label in labels]


    fq_pathes = R1_fq_files + R2_fq_files
    fq_files_qname_changed = R1_fq_files_qname_changed + R2_fq_files_qname_changed
    full_unmapped_fq_path = R1_full_unmap_fq + R2_full_unmap_fq
    full_sam_path = R1_full_sam_pathes + R2_full_sam_pathes
    RE_fq_pathes = R1_RE_fq_pathes + R2_RE_fq_pathes
    RE_sam_pathes = R1_RE_sam_pathes + R2_RE_sam_pathes
    RE_summary_log = R1_RE_summary_log + R2_RE_summary_log
    map_cpu = int(cpu/2)

    pair_intra_pets_pathes = [pair_dir + os.sep + label + '.intra.pets' for label in labels]
    pair_inter_pets_pathes = [pair_dir + os.sep + label + '.inter.pets' for label in labels]
    pair_log_pathes = [pair_dir + os.sep + label + '.paired.log' for label in labels]
    pair_log_summary = output_dir + os.sep + 'paired_summay.txt'

    ### full length mapping
    if start_full_mapping:
        ### change fq qnames
        time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Fastq file change qname!"
        print time
        os.system('echo %s >> %s ' % (time, hic_log_file))
        Parallel(n_jobs=len(fq_pathes))(delayed(qname)(fq_pathes[x], fq_files_qname_changed[x]) for x in range(len(fq_pathes)))

        ### full length mapping
        time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Full length mapping!"
        print time
        os.system('echo %s >> %s ' % (time, hic_log_file))
        for i in range(len(R1_fq_files)):
            fq_path = [R1_fq_files_qname_changed[i], R2_fq_files_qname_changed[i]]
            unmapped_fq = [R1_full_unmap_fq[i], R2_full_unmap_fq[i]]
            sam_path = [R1_full_sam_pathes[i], R2_full_sam_pathes[i]]
            bowtie_log = [R1_full_bowtie_log[i], R2_full_bowtie_log[i]]
            Parallel(n_jobs=len(fq_path))(delayed(full_mapping)(fq_path[x], bowtie_index, map_cpu, mismatch, trim5, trim3,
                                                     unmapped_fq[x], sam_path[x], bowtie_log[x]) for x in range(len(fq_path)))
        os.system('rm %s -rf' % (' '.join(fq_files_qname_changed)))
    if start_search_REsite:
    ### unmaped reads find REsites and cutted by REsite
        time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Searching REsite!"
        print time
        os.system('echo %s >> %s ' % (time, hic_log_file))
        Parallel(n_jobs=len(full_unmapped_fq_path))(delayed(REsite)(full_unmapped_fq_path[x], sites5, RE_fq_pathes[x], trim5,trim3, RE_summary_log[x])
                                                for x in range(len(full_unmapped_fq_path)))
        os.system('rm %s -rf' % (' '.join(full_unmapped_fq_path)))

    if start_REsite_mapping:
    ### reads has REsite mapping
        time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " RE cutted mapping!"
        print time
        os.system('echo %s >> %s ' % (time, hic_log_file))
        for i in range(len(R1_RE_fq_pathes)):
            RE_fq_path = [R1_RE_fq_pathes[i], R2_RE_fq_pathes[i]]
            sam_path = [R1_RE_sam_pathes[i], R2_RE_sam_pathes[i]]
            bowtie_log = [R1_RE_bowtie_log[i], R2_RE_bowtie_log[i]]
            Parallel(n_jobs=2)(delayed(RE_mapping)(RE_fq_path[x], bowtie_index, map_cpu, mismatch, sam_path[x], bowtie_log[x])
                               for x in range(2))
    os.system('rm %s -rf' % (' '.join(RE_fq_pathes)))

    ### divided reads into small data set files
    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Splitting mapped reads!"
    print time
    os.system('echo %s >> %s ' % (time, hic_log_file))
    Parallel(n_jobs=len(R1_full_sam_pathes))(delayed(sam_split)(R1_full_sam_pathes[i], "full", "R1", pair_sub_dirs[i], split_cutoff)
                                                 for i in range(len(R1_full_sam_pathes)))
    Parallel(n_jobs=len(R2_full_sam_pathes))(delayed(sam_split)(R2_full_sam_pathes[i], "full", "R2", pair_sub_dirs[i], split_cutoff)
                                                 for i in range(len(R2_full_sam_pathes)))
    Parallel(n_jobs=len(R1_RE_sam_pathes))(delayed(sam_split)(R1_RE_sam_pathes[i], "RE", "R1", pair_sub_dirs[i], split_cutoff)
                                                 for i in range(len(R1_RE_sam_pathes)))
    Parallel(n_jobs=len(R2_RE_sam_pathes))(delayed(sam_split)(R2_RE_sam_pathes[i], "RE", "R2", pair_sub_dirs[i], split_cutoff)
                                                 for i in range(len(R2_RE_sam_pathes)))

    ### pairing
    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Pairing!"
    print time
    os.system('echo %s >> %s ' % (time, hic_log_file))
    paired_type_all_log = {}
    for i in range(len(pair_sub_dirs)):
        label = labels[i]
        pair_sub_dir = pair_sub_dirs[i]
        pair_intra_pets_path = pair_intra_pets_pathes[i]
        pair_inter_pets_path = pair_inter_pets_pathes[i]
        pair_log_path =pair_log_pathes[i]
        hit_files = os.listdir(pair_sub_dir)
        hit_files = [pair_sub_dir + os.sep + hit_file for hit_file in hit_files if hit_file.startswith('hits.')]
        paired_type_each_file = Parallel(n_jobs=cpu)(delayed(pairing)(hit_file, genomic_dis, pair_sub_dir)
                             for hit_file in hit_files)
        ### sum each paired log in each sub dir
        paired_type_each_label = paired_type_each_file[0]
        for x in range(1,len(paired_type_each_file)):
            for key in paired_type_each_label.keys():
                paired_type_each_label[key] += paired_type_each_file[x][key]
        paired_type_all_log[label] = paired_type_each_label
        ### merge intra or inter files for each sample
        pair_files = os.listdir(pair_sub_dir)
        pair_intra_files = [pair_sub_dir + os.sep + file for file in pair_files if file.startswith("paired.") and file.endswith('.intra.pets')]
        pair_inter_files = [pair_sub_dir + os.sep + file for file in pair_files if file.startswith("paired.") and file.endswith('.inter.pets')]
        os.system('cat %s > %s' % (' '.join(pair_intra_files), pair_intra_pets_path))
        os.system('cat %s > %s' % (' '.join(pair_inter_files), pair_inter_pets_path))
        #os.system('rm %s' % pair_dir)
    ### sum each paired log for each sample
    paired_type_names = ["intra_pets", "inter_pets", "genomic",
                         "R1full_R2full_intra", "R1full_R2RE_intra", "R1RE_R2full_intra",
                          "R1RE_R2RE_intra", 'R1RE_R1RE_intra', 'R2RE_R2RE_intra',
                          "R1full_R2full_inter", "R1full_R2RE_inter", "R1RE_R2full_inter",
                          "R1RE_R2RE_inter", 'R1RE_R1RE_inter', 'R2RE_R2RE_inter',
                          "R1full_unpaired", "R1RE_unpaired", "R2full_unpaired", "R2RE_unpaired"] \
                        + ["reads_has_%d_pets" % x for x in range(8)]
    paired_type_all_sum = {}
    pair_log_summary = open(pair_log_summary, 'w')
    header = 'paired_type\ttotal\t' + '\t'.join(labels) + '\n'
    pair_log_summary.write(header)
    for key in paired_type_names:
        each_label_nums = []
        for label in labels:
            each_label_nums.append(paired_type_all_log[label][key])
        each_label_sum = sum(each_label_nums)
        paired_type_all_sum[key] = each_label_sum
        pair_log_summary.write(key + '\t' + str(each_label_sum) + '\t' + '\t'.join([str(x) for x in each_label_nums]) + '\n')
    pair_log_summary.close()

    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Split pets by chr!"
    print time
    os.system('echo %s >> %s ' % (time, hic_log_file))
    for pets_path in pair_intra_pets_pathes:
        split_by_chr(pets_path, pets_intra_dir)
    for pets_path in pair_inter_pets_pathes:
        split_by_chr(pets_path, pets_inter_dir)

    time = datetime.strftime(datetime.now(), '%Y-%m-%d %H:%M:%S') + " Remove duplicate!"
    print time
    os.system('echo %s >> %s ' % (time, hic_log_file))
    intra_pets_files = os.listdir(pets_intra_dir)
    intra_pets_files = [pets_intra_dir + os.sep + file for file in intra_pets_files if file.endswith('.pets.no_rmdup')]
    intra_keep_numbers = Parallel(n_jobs=cpu)(delayed(pets_remove_dumplicate)(pets_file) for pets_file in intra_pets_files)
    os.system('rm %s' % (' '.join(intra_pets_files)))
    inter_pets_files = os.listdir(pets_inter_dir)
    inter_pets_files = [pets_inter_dir + os.sep + file for file in inter_pets_files if file.endswith('.pets.no_rmdup')]
    inter_keep_numbers =Parallel(n_jobs=cpu)(delayed(pets_remove_dumplicate)(pets_file) for pets_file in inter_pets_files)
    os.system('rm %s' % (' '.join(inter_pets_files)))
    intra_keep_numbers = sum(intra_keep_numbers)
    inter_keep_numbers = sum(inter_keep_numbers)

    ### write the full mapping, RE_site for hic_log_file
    R1_raw_rds = os.popen('cat %s|grep processed|tr ":" "\t"|cut -f 2 |awk \'{a+=$1}END{print a}\''
                          % (' '.join(R1_full_bowtie_log))).readlines()[0].rstrip('\n')
    R2_raw_rds = os.popen('cat %s|grep processed|tr ":" "\t"|cut -f 2 |awk \'{a+=$1}END{print a}\''
                          % (' '.join(R2_full_bowtie_log))).readlines()[0].rstrip('\n')
    R1_full_mapped_rds = os.popen('cat %s|grep least|tr ":" "\t"|cut -f 2|tr " " "\t" |cut -f 2|awk \'{a+=$1}END{print a}\''
                          % (' '.join(R1_full_bowtie_log))).readlines()[0].rstrip('\n')
    R2_full_mapped_rds = os.popen('cat %s|grep least|tr ":" "\t"|cut -f 2|tr " " "\t" |cut -f 2|awk \'{a+=$1}END{print a}\''
                          % (' '.join(R2_full_bowtie_log))).readlines()[0].rstrip('\n')
    R1_raw_rds = int(R1_raw_rds)
    R2_raw_rds = int(R2_raw_rds)
    R1_full_mapped_rds = int(R1_full_mapped_rds)
    R2_full_mapped_rds = int(R2_full_mapped_rds)
    R1_full_mapped_ratio = float(R1_full_mapped_rds) * 100 / R1_raw_rds
    R2_full_mapped_ratio = float(R2_full_mapped_rds) * 100 / R2_raw_rds

    R1_REsite_summary = os.popen('cat %s|grep -v REsite|awk \'BEGIN{FS=OFS="\t"}{a[$1]+=$2;b[$1]+=$3;c[$1]+=$4}'
                                 'END{for (i in a) print i,a[i],b[i],c[i]}\''% (' '.join(R1_RE_summary_log))).readlines()
    R1_REsite_summary = [x.rstrip('\n').split('\t') for x in R1_REsite_summary]
    R1_REsite_summary = {x[0] : [int(x[1]), int(x[2]), int(x[3])] for x in R1_REsite_summary}
    R2_REsite_summary = os.popen('cat %s|grep -v REsite|awk \'BEGIN{FS=OFS="\t"}{a[$1]+=$2;b[$1]+=$3;c[$1]+=$4}'
                                 'END{for (i in a) print i,a[i],b[i],c[i]}\''% (' '.join(R2_RE_summary_log))).readlines()
    R2_REsite_summary = [x.rstrip('\n').split('\t') for x in R2_REsite_summary]
    R2_REsite_summary = {x[0] : [int(x[1]), int(x[2]), int(x[3])] for x in R2_REsite_summary}

    R1_unmaped_rds_has_RE = sum([x[0] for x in R1_REsite_summary.values()])
    R2_unmaped_rds_has_RE = sum([x[0] for x in R2_REsite_summary.values()])
    R1_RE_cutted_rds_keep = sum([x[1] for x in R1_REsite_summary.values()]) + sum([x[2] for x in R1_REsite_summary.values()])
    R2_RE_cutted_rds_keep = sum([x[1] for x in R2_REsite_summary.values()]) + sum([x[2] for x in R2_REsite_summary.values()])
    #R1_RE_mapped
    R1_RE_mapped_rds = os.popen('cat %s|grep least|tr ":" "\t"|cut -f 2|tr " " "\t" |cut -f 2|'
                                'awk \'{a+=$1}END{print a}\'' % ' '.join(R1_RE_bowtie_log)).readlines()[0].rstrip('\n')
    R2_RE_mapped_rds = os.popen('cat %s|grep least|tr ":" "\t"|cut -f 2|tr " " "\t" |cut -f 2|'
                                'awk \'{a+=$1}END{print a}\'' % ' '.join(R2_RE_bowtie_log)).readlines()[0].rstrip('\n')
    R1_RE_mapped_rds = int(R1_RE_mapped_rds)
    R2_RE_mapped_rds = int(R2_RE_mapped_rds)
    R1_rds_has_RE_ratio = float(R1_unmaped_rds_has_RE) * 100 / R1_raw_rds
    R2_rds_has_RE_ratio = float(R2_unmaped_rds_has_RE) * 100 / R2_raw_rds
    R1_RE_mapped_ratio = float(R1_RE_mapped_rds) * 100 / R1_RE_cutted_rds_keep
    R2_RE_mapped_ratio = float(R2_RE_mapped_rds) * 100 / R2_RE_cutted_rds_keep
    intra_rm_dup_ratio = float(intra_keep_numbers) * 100 / paired_type_all_sum['intra_pets']
    inter_rm_dup_ratio = float(inter_keep_numbers) * 100 / paired_type_all_sum['inter_pets']

    hic_log_file = open(hic_log_file,'a')
    hic_log_file.write('Mapping info:\n')
    hic_log_file.write('reads\traw\tfull_mapped\trds_has_RE\trds_RE_cutted\tRE_mapped\n')
    hic_log_file.write('R1\t%d\t%d(%.2f%%)\t%d(%.2f%%)\t%d\t%d(%.2f%%)\n'
                       %(R1_raw_rds,R1_full_mapped_rds,R1_full_mapped_ratio,R1_unmaped_rds_has_RE,R1_rds_has_RE_ratio,R1_RE_cutted_rds_keep,R1_RE_mapped_rds,R1_RE_mapped_ratio))
    hic_log_file.write('R2\t%d\t%d(%.2f%%)\t%d(%.2f%%)\t%d\t%d(%.2f%%)\n\n'
                       %(R2_raw_rds,R2_full_mapped_rds,R2_full_mapped_ratio,R2_unmaped_rds_has_RE,R2_rds_has_RE_ratio,R2_RE_cutted_rds_keep,R2_RE_mapped_rds,R2_RE_mapped_ratio))
    hic_log_file.write('REsite summary:\n')
    hic_log_file.write('REsite\tR1_rds\tR1_5_longer_20bp\tR1_3_longer_20bp\tR2_rds\tR2_5_longer_20bp\tR2_3_longer_20bp\n')
    for key in R1_REsite_summary:
        hic_log_file.write(key + '\t' + '\t'.join([str(x) for x in R1_REsite_summary[key]]) +
                           '\t' + '\t'.join([str(x) for x in R2_REsite_summary[key]]) + '\n')
    hic_log_file.write('\n')
    hic_log_file.write('Pair summary:\n')
    hic_log_file.write('genomic\tintra_pets\tinter_pets\tintra_rm_dup\tinter_rm_dup\n')
    hic_log_file.write('%d\t%d\t%d\t%d(%.2f%%)\t%d(%.2f%%)\n\n' % (paired_type_all_sum['genomic'], paired_type_all_sum['intra_pets'],
                    paired_type_all_sum['inter_pets'], intra_keep_numbers, intra_rm_dup_ratio, inter_keep_numbers,inter_rm_dup_ratio))
    #os.system('rm %s -rf' % (' '.join(pair_intra_pets_pathes + pair_inter_pets_pathes + pair_sub_dirs)))

    hic_log_file.write('or unit: million\n')
    hic_log_file.write('Mapping info:\n')
    hic_log_file.write('reads\traw/M\tfull_mapped/M\trds_has_RE/M\trds_RE_cutted/M\tRE_mapped/M\n')
    hic_log_file.write('R1\t%d\t%d(%.2f%%)\t%d(%.2f%%)\t%d\t%d(%.2f%%)\n'
                       % (int(round(float(R1_raw_rds)/1000000)), int(round(float(R1_full_mapped_rds)/1000000)), R1_full_mapped_ratio,
                        int(round(float(R1_unmaped_rds_has_RE)/1000000)), R1_rds_has_RE_ratio,
                        int(round(float(R1_RE_cutted_rds_keep)/1000000)), int(round(float(R1_RE_mapped_rds)/1000000)) ,R1_RE_mapped_ratio))
    hic_log_file.write('R2\t%d\t%d(%.2f%%)\t%d(%.2f%%)\t%d\t%d(%.2f%%)\n'
                       % (int(round(float(R2_raw_rds)/1000000)), int(round(float(R2_full_mapped_rds)/1000000)), R2_full_mapped_ratio,
                        int(round(float(R2_unmaped_rds_has_RE)/1000000)), R2_rds_has_RE_ratio,
                        int(round(float(R2_RE_cutted_rds_keep)/1000000)), int(round(float(R2_RE_mapped_rds)/1000000)) ,R2_RE_mapped_ratio))

    hic_log_file.write('REsite summary:\n')
    hic_log_file.write('REsite\tR1_rds\tR1_5_longer_20bp\tR1_3_longer_20bp\tR2_rds\tR2_5_longer_20bp\tR2_3_longer_20bp\n')
    for key in R1_REsite_summary:
        hic_log_file.write(key + '\t' + '\t'.join([str(int(round(float(x)/1000000))) for x in R1_REsite_summary[key]]) +
                           '\t' + '\t'.join([str(int(round(float(x)/1000000))) for x in R2_REsite_summary[key]]) + '\n')
    hic_log_file.write('\n')
    hic_log_file.write('Pair summary:\n')
    hic_log_file.write('genomic/M\tintra_pets/M\tinter_pets/M\tintra_rm_dup/M\tinter_rm_dup/M\n')
    hic_log_file.write('%d\t%d\t%d\t%d(%.2f%%)\t%d(%.2f%%)\n' % (int(round(paired_type_all_sum['genomic']/1000000)),
                        int(round(paired_type_all_sum['intra_pets']/1000000)),int(round(paired_type_all_sum['inter_pets']/1000000)),
                        int(round(intra_keep_numbers/1000000)),intra_rm_dup_ratio,int(round(inter_keep_numbers/1000000)),inter_rm_dup_ratio))
    hic_log_file.close()

    ### remove files
    os.system('rm %s -rf' % (pair_dir))
