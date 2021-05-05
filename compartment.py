#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy import linalg as LA
from joblib import Parallel, delayed
import os
import getopt
import sys
from datetime import datetime

def bins_num_and_header(chrsize_path, bin_size):
    headers = {}
    chr_list = []
    f = open(chrsize_path)
    for line in f:
        line = line.rstrip('\n').split('\t')
        chr = line[0]
        one_chr_size = int(line[1])
        bin_num = one_chr_size / bin_size
        bin_num += 1
        one_chr_header = [chr + ':1-' + str(bin_size - 1)]
        one_chr_header = one_chr_header + \
                         [chr + ':' + str(i * bin_size) + '-'
                          + str((i + 1) * bin_size - 1) for i in xrange(1, bin_num - 1)]
        one_chr_header = one_chr_header + \
                         [chr + ':' +
                          str((bin_num - 1) * bin_size) + '-' + str(one_chr_size)]
        headers[chr] = one_chr_header
        chr_list.append(chr)
    f.close()
    return headers, chr_list

def read_matrix(matrix_path, weather_header):
    matrix = []
    f = open(matrix_path)
    for line in f:
        line = line.rstrip('\n').split('\t')
        matrix.append(line)
    f.close()
    matrix = np.array(matrix)
    if weather_header:
        matrix = matrix[1:, 1:]
    matrix = matrix.astype(float)
    matrix = matrix.tolist()
    return matrix

def observed_expected_intra(matrix):
    bin_num = len(matrix[0])
    expect_nums = [float(0) for x in range(bin_num)]
    for i in range(bin_num):
        for j in range(i, bin_num):
            idx = j - i
            expect_nums[idx] += float(matrix[i][j])
    for i in range(bin_num):
        d = bin_num - i
        expect_num_d = expect_nums[i] / d
        expect_nums[i] = expect_num_d
    for i in range(bin_num):
        for j in range(bin_num):
            if i == j:
                matrix[i][j] = 0
            else:
                d = abs(j - i)
                expect_num_d = expect_nums[d]
                if expect_num_d == 0:
                    expect_num_d = 1
                matrix[i][j] = matrix[i][j] / expect_num_d
    return matrix

def pcc(matrix):
    matrix = np.array(matrix, dtype='float32')
    np.seterr(divide='ignore', invalid='ignore')
    pcc_matrix = np.corrcoef(matrix)
    pcc_matrix = np.nan_to_num(pcc_matrix)
    pcc_matrix = pcc_matrix.tolist()
    return pcc_matrix

def PCA(data):
    """
    returns: data transformed in 2 dims/columns + regenerated original data
    pass in: data as 2D NumPy array
    """
    m, n = data.shape
    # mean center the data
    data -= data.mean(axis=0)
    # calculate the covariance matrix
    R = np.cov(data, rowvar=False)
    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since R is symmetric,
    # the performance gain is substantial
    evals, evecs = LA.eigh(R)
    # sort eigenvalue in decreasing order
    idx = np.argsort(evals)[::-1]
    evecs = evecs[:, idx]
    # sort eigenvectors according to same index
    evals = evals[idx]
    # select the first n eigenvectors (n is desired dimension
    # of rescaled data array, or dims_rescaled_data)
    evecs_1 = evecs[:, :1]
    evecs_2 = evecs[:, 1:2]
    # carry out the transformation on the data using eigenvectors
    # and return the re-scaled data, eigenvalues, and eigenvectors
    pca_1_data = np.dot(evecs_1.T, data.T).T
    pca_2_data = np.dot(evecs_2.T, data.T).T
    return pca_1_data, pca_2_data

def pca_cal_for_matrix_remove_zero_row_cols(pcc_matrix):
    pcc_matrix = np.array(pcc_matrix, dtype=np.float)
    pc1_data = np.zeros(len(pcc_matrix[0]))
    pc2_data = np.zeros(len(pcc_matrix[0]))
    col_sum = np.sum(pcc_matrix, axis=0)
    idxes_zero = np.where(col_sum == 0)
    idxes_nozero = np.where(col_sum != 0)
    pcc_matrix = np.delete(pcc_matrix, idxes_zero, axis=0)
    pcc_matrix = np.delete(pcc_matrix, idxes_zero, axis=1)
    pca_pc1_pc2_data = PCA(pcc_matrix)
    pca_1_data = pca_pc1_pc2_data[0]
    pca_2_data = pca_pc1_pc2_data[1]
    pc1_data[idxes_nozero] = pca_1_data.T
    pc2_data[idxes_nozero] = pca_2_data.T
    pc1_data = pc1_data.tolist()
    pc2_data = pc2_data.tolist()
    return pc1_data, pc2_data

def wirte_matix(matrix, header_one_chr, save_path):
    f = open(save_path, 'w')
    f.write('region\t' + '\t'.join(header_one_chr) + '\n')
    for i in xrange(len(matrix)):
        f.write(header_one_chr[i] + '\t' + '\t'.join([str(x) for x in matrix[i]]) + '\n')
    f.close()

def wirte_pca(pca_data, header_one_chr, save_path):
    f = open(save_path, 'w')
    for i in xrange(len(pca_data)):
        f.write(header_one_chr[i].replace("-", "\t").replace(":", "\t") + '\t' + str(pca_data[i]) + '\n')
    f.close()

def compartment_main(matrix_path, weather_header, header_one_chr, save_dir):
    file_name = os.path.basename(matrix_path).replace(".matrix", "")
    chr = file_name.split('.')[0]

    matrix = read_matrix(matrix_path, weather_header)
    ob_ex_matrix = observed_expected_intra(matrix)
    pcc_matrix = pcc(ob_ex_matrix)
    pca_data = pca_cal_for_matrix_remove_zero_row_cols(pcc_matrix)
    pca_1_data = pca_data[0]
    pca_2_data = pca_data[1]

    ob_ex_path = save_dir + os.sep + "ob_ex" + os.sep + chr + '.mat'
    pcc_path = save_dir + os.sep + "pcc" + os.sep + chr + '.mat'
    pca_1_path = save_dir + os.sep + "compartment_pc1" + os.sep + chr + '.pca.bdg'
    pca_2_path = save_dir + os.sep + "compartment_pc2" + os.sep + chr + '.pca.bdg'

    wirte_matix(ob_ex_matrix, header_one_chr, ob_ex_path)
    wirte_matix(pcc_matrix, header_one_chr, pcc_path)
    wirte_pca(pca_1_data, header_one_chr, pca_1_path)
    wirte_pca(pca_2_data, header_one_chr, pca_2_path)

def usage():
    print ' Usage: Calculate compartment.Calculate observed frequency and expected frequency.'
    print '     -i      The input matrix dir.Matirx files must endswith .matrix or .mat.'
    print '     -o      The output dir.Has 2 subdir:ob_ex,pcc.'
    print '     -c      The chr_size_file to make header.'
    print '     -d      The input matrix weather has header.Default:True. If set: False.'
    print '     -p      The cpu number.'
    print '     -b      The active histone peaks bed file to define compartment. Default: No this step.'
    print '             If set, do this step!'
    print '     -s      The binsize. Default:40000.'

if __name__ == '__main__':
    start_time = datetime.now()
    weather_header = True
    cpu = 1
    bin_size = 40000
    bed_file = False
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:o:c:dp:b:s:')
        for o, a in opts:
            if o == '-h':
                usage()
                sys.exit()
            if o == '-i':
                input_dir = a
            if o == '-o':
                save_dir = a
            if o == '-c':
                chr_size_path = a
            if o == '-d':
                weather_header = False
            if o == '-p':
                cpu = int(a)
            if o == '-b':
                bed_file = a
            if o == '-s':
                bin_size = int(a)
    except getopt.GetoptError:
        print 'Error in getting parametres!'
        sys.exit()

    ob_ex_dir = save_dir + os.sep + "ob_ex"
    pcc_dir = save_dir + os.sep + "pcc"
    pca_1_dir = save_dir + os.sep + "compartment_pc1"
    pca_2_dir = save_dir + os.sep + "compartment_pc2"

    if len(save_dir) == 0:
        save_dir = '.'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    if not os.path.exists(ob_ex_dir):
        os.makedirs(ob_ex_dir)
    if not os.path.exists(pcc_dir):
        os.makedirs(pcc_dir)
    if not os.path.exists(pca_1_dir):
        os.makedirs(pca_1_dir)
    if not os.path.exists(pca_2_dir):
        os.makedirs(pca_2_dir)

    headers, chr_list = bins_num_and_header(chr_size_path, bin_size)

    all_files = os.listdir(input_dir)
    matrix_files_name = [file for file in all_files if file.endswith('.matrix') or file.endswith(".mat")]
    matrix_chrs = [f.split('.')[0] for f in matrix_files_name]
    chrs = [chr for chr in matrix_chrs if chr in chr_list]
    chrs = sorted(chrs)
    input_matrix_pathes = []
    pca_1_pathes = []
    pca_2_pathes = []
    for chr in chrs:
        file_name = [f for f in matrix_files_name if chr + '.' in f][0]
        file = input_dir + os.sep + file_name
        input_matrix_pathes.append(file)

        pca_1_path = save_dir + os.sep + "compartment_pc1" + os.sep + chr + '.pca.bdg'
        pca_1_pathes.append(pca_1_path)
        pca_2_path = save_dir + os.sep + "compartment_pc2" + os.sep + chr + '.pca.bdg'
        pca_2_pathes.append(pca_2_path)
    pca_pathes = [pca_1_pathes, pca_2_pathes]
    Parallel(n_jobs=cpu)(delayed(compartment_main)(input_matrix_pathes[x], weather_header, headers[chrs[x]], save_dir)
                         for x in range(len(input_matrix_pathes)))

    compartment_pc1_pc2_no_define_AB_pathes = [save_dir + os.sep + 'compartment_pc' + str(x) + '_no_define_AB.bdg'
                                         for x in [1, 2]]
    for f in compartment_pc1_pc2_no_define_AB_pathes:
        if os.path.isfile(f):
            os.remove(f)

    ### Merge pca files
    commands = ['cat %s > %s' % (' '.join(pca_pathes[x]), compartment_pc1_pc2_no_define_AB_pathes[x])
                for x in range(2)]
    Parallel(n_jobs=len(commands))(delayed(os.system)(command) for command in commands)

    os.system("rm %s -r" % (pca_1_dir))
    os.system("rm %s -r" % (pca_2_dir))

    if bed_file:
        compartment_pathes = [save_dir + os.sep + 'compartment_pc' + str(x) + '.bdg' for x in [1, 2]]
        for f in compartment_pathes:
            if os.path.isfile(f):
                os.remove(f)
        define_commands = []
        for i in range(2):
            com_in = compartment_pc1_pc2_no_define_AB_pathes[i]
            com_out = compartment_pathes[i]
            command = "chrs=( %s ); for chr in ${chrs[@]};do abstract=`grep -w $chr %s|" \
                  "bedtools intersect -a - -b %s -wao|" \
                  "perl -alne 'print \"$F[0]\t$F[1]\t$F[2]\t$F[3]\" if $F[5] ne -1'|" \
                  "awk 'BEGIN{a=0;b=0}{if ($4 > 0){a+=1}}{if ($4 < 0){b+=1}}END{print a,b}'`;" \
                  "a=`echo $abstract|tr ' ' '\t'|cut -f 1`;b=`echo $abstract|tr ' ' '\t'|cut -f 2`;" \
                  "if [ $a -gt $b ]; then grep -w $chr %s >> %s;" \
                  "else grep -w $chr %s |perl -alne 'print \"$F[0]\t$F[1]\t$F[2]\t\".($F[3]*(-1))' >> %s;fi;done" \
                  % (' '.join(chrs), com_in, bed_file, com_in, com_out, com_in, com_out)
            define_commands.append(command)
        Parallel(n_jobs=len(commands))(delayed(os.system)(command) for command in define_commands)
        os.system("rm %s" % ' '.join(compartment_pc1_pc2_no_define_AB_pathes))

        compartment_bw = [save_dir + os.sep + 'compartment_pc' + str(x)+'.bw' for x in [1,2]]
        commands = ["bedGraphToBigWig %s %s %s" % (compartment_pathes[x], chr_size_path, compartment_bw[x]) for x in range(2)]
        Parallel(n_jobs=len(commands))(delayed(os.system)(command) for command in commands)
    else:
        compartment_bw = [x.replace(".bdg", ".bw") for x in compartment_pc1_pc2_no_define_AB_pathes]
        commands = ["bedGraphToBigWig %s %s %s"
                    % (compartment_pc1_pc2_no_define_AB_pathes[x], chr_size_path, compartment_bw[x])
                    for x in range(2)]
        Parallel(n_jobs=len(commands))(delayed(os.system)(command) for command in commands)

    time_caused = datetime.now() - start_time
    print "The time caused is : ", time_caused, '\n'
