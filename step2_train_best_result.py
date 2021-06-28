import glob
import multiprocessing
import os
from itertools import combinations
import pandas as pd
import numpy as np
import re

from ml_auc import start
from parameter import KMER_K, HASH_SIZE, RANDOM, METHOD, methods, CORE
from prepare_ncbi_homolog_neg_pair import get_protein_list, neg_homologous


def extract_pos_neg(protein_a,protein_b, train_dir,pos_path,rd =123):

    raw_data = glob.glob(pos_path + "*.raw")
    all_name = []
    all_pos = []
    all_a_list = get_protein_list(protein_a)
    all_b_list = get_protein_list(protein_b)
    all_a = set(all_a_list)

    all_name.extend(all_a_list)
    all_name.extend(all_b_list)
    all_name = set(all_name)

    for each in raw_data:
        database_name = os.path.basename(each).split("_")[0]
        name_a = os.path.basename(protein_a).split("_")[1].split(".")[0]
        name_b = os.path.basename(protein_b).split("_")[1].split(".")[0]
        w_name = os.path.join(train_dir,database_name + "-" + name_a + "_" + name_b +'.txt')
        w_path = open(w_name, "w")
        with open(each) as f:
            for line in f:
                tem = line.strip().split("\t")
                if (tem[0] in all_name and tem[1] in all_name):
                    # order must be: protein_a ---------> protein_b
                    if (tem[0] in all_a):
                        w_path.write(str(1) + "\t" + tem[0] + '\t' + tem[1] + '\n')
                    else:
                        w_path.write(str(1) + "\t" + tem[1] + '\t' + tem[0] + '\n')
        w_path.close()
        all_pos.append(w_name)
    intersection_res = intersection_raw(all_pos)
    union_res = union_raw(all_pos)
    all_pos.append(intersection_res)
    all_pos.append(union_res)

    all_combine = []
    for line in all_pos:
        # print(line)
        pos_filter = line
        # must protein_a -------> protein_b
        neg_file_name = neg_homologous(protein_a, protein_b, pos_filter,rd)
        combine_name = neg_file_name.split(".")[0] + "_pos.tsv"
        cmd = "cat " + pos_filter + " " + neg_file_name + " > " + combine_name
        os.system(cmd)
        all_combine.append(combine_name)

    return all_combine



def combine(feature_name,label_name):
    feature_db = {}
    # combine result


    #feature file
    for line in feature_name:
        with open(line) as f:
            for line in f:
                tem = line.strip().split("\t")
                feature_db[tem[0]] = tem[1]

    # pos and neg train data
    all_train_file =[]
    for each in label_name:
        tem = each.split("-")
        train_file_name = tem[0] + "-" + tem[1] + ".train"
        w_ = open(train_file_name,'w')
        all_train_file.append(train_file_name)
        with open(each) as f:
            for line in f:
                tem = line.strip().split('\t')
                w_.write(str(tem[0]) +',' + feature_db[tem[1]] +','+ feature_db[tem[2]] + '\n')
        w_.close()
    return all_train_file

def intersection_raw(file):

    base_dir = os.path.dirname(file[0])
    base_name = os.path.basename(file[0]).split("-")[1]
    w_name = os.path.join(base_dir,"intersection-" + base_name)
    w_path = open(w_name, 'w')

    first = []
    with open(file[0]) as f:
        for line in f:
            first.append(line)
    file.remove(file[0])

    for each in file:
        each_pairs = set()
        with open(each) as f:
            for line in f:
                each_pairs.add(line)
        res = [x for x in first if x in each_pairs]
        first = res

    for line2 in first:
        w_path.write(line2)
    w_path.close()
    return w_name

def union_raw(file):

    base_dir = os.path.dirname(file[0])
    base_name = os.path.basename(file[0]).split("-")[1]

    w_name = os.path.join(base_dir,"union-" + base_name)
    w_path = open(w_name, 'w')

    union = set()

    for each in file:
        with open(each) as f:
            for line in f:
                if line not in union:
                    w_path.write(line)
                    union.add(line)
    w_path.close()
    return w_name

def run(pairs):
    # get pos and neg train
    all_combine = extract_pos_neg(pairs[0],pairs[1],train_dir,pos_path,RANDOM)
    feature_name = [var.split(".")[0] + "_k" + str(KMER_K) + "_sim" + str(HASH_SIZE) + ".feature" for var in pairs]
    # combine feature and label
    train_file = combine(feature_name,all_combine)

    # train auc
    all_auc = []
    method = METHOD
    method_name = "_".join(method.split(" "))
    w_names = re.split("-|\.",train_file[0])[1]
    res = open(res_dir+w_names+"_ROC.auc","w")
    for each_train in train_file:
        for each_method in methods[-2:]:
            train_data = each_train
            data = pd.read_csv(train_data,sep=',', header=None)
            dataCol = data.columns.values.tolist()[1:]
            data_x = np.array(data[dataCol])
            data_y = np.array(data[0])
            w_path = res_dir
            auc_plot_name = os.path.basename(each_train).split(".")[0]
            mean_acu = start(each_method,w_path,data_x,data_y,auc_plot_name)
            print(each_train + "\t"+ str(mean_acu)[:5] + "\t" + each_method+ "\n")
            res.write(each_train + "\t"+ str(mean_acu)[:5] + "\t" + each_method + "\n")
    res.close()

if __name__ == '__main__':
    pos_path = "/home/yangfang/Sash/QFO/test_ML/data_sort/"
    train_dir = "/home/yangfang/Sash/QFO/test_ML/test_train/"
    res_dir = "/home/yangfang/Sash/QFO/test_ML/test_res/"
    files = glob.glob("/home/yangfang/Sash/QFO/QfO_release_2018_04/protein_all_trans_test/*.fasta")
    os.makedirs(train_dir,exist_ok=True)
    os.makedirs(res_dir,exist_ok=True)
    combination_file = list(combinations(files,2))
    #filter already run
    a_run = set()
    all_value_file = glob.glob(res_dir+"*.auc")
    print(len(all_value_file))
    for line in all_value_file:
        tem_name = os.path.basename(line)
        a_run.add(tem_name.split("_ROC")[0])
    all_pairs = []
    for line2 in combination_file:
        name_aa = os.path.basename(line2[0]).split("_")[1].split(".")[0]
        name_bb = os.path.basename(line2[1]).split("_")[1].split(".")[0]
        tem = name_aa + "_" + name_bb
        if tem not in a_run:
            all_pairs.append(line2)
    print("total: " + str(len(combination_file)))
    print("already run: " + str(len(a_run)))
    print("will run: " + str(len(all_pairs)))
    # protein_a = "/home/yangfang/Sash/QFO/QfO_release_2018_04/protein_all_trans_test/UP000000425_122586.fasta"
    # protein_b = "/home/yangfang/Sash/QFO/QfO_release_2018_04/protein_all_trans_test/UP000000429_85962.fasta"
    # Do with multiprocess

    # all_pairs = list(combinations(files,2))
    pool = multiprocessing.Pool(processes=CORE)
    pool.map(run,all_pairs)