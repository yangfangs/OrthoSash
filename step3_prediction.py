import glob
import multiprocessing
import os

import numpy as np
import pandas as pd

from sklearn import svm, neighbors, tree
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier

from sklearn.naive_bayes import BernoulliNB
from parameter import RANDOM, CORE


def start(method, save_path, X, y, pairs):
    all_feature_A = []
    all_feature_B = []
    P_A = pairs.split("_")[0]
    P_B = pairs.split("_")[1]
    all_feature_file = glob.glob(feature_path)
    P_A_feaure = [var for var in all_feature_file if os.path.basename(var).split("_")[1] == P_A][0]
    P_B_feaure = [var for var in all_feature_file if os.path.basename(var).split("_")[1] == P_B][0]
    with open(P_A_feaure) as f:
        for line in f:
            all_feature_A.append(line.strip())

    with open(P_B_feaure) as f2:
        for line2 in f2:
            all_feature_B.append(line2.strip())

    if method == "SVC":
        clf = svm.SVC(probability=True, gamma='auto', random_state=RANDOM, max_iter=10000)
    elif method == "Random_Forest_Classifier":
        clf = RandomForestClassifier(n_estimators=1000, max_features="sqrt", random_state=RANDOM, n_jobs=1)
    elif method == "Gradient_Boosting_Classifier":
        clf = GradientBoostingClassifier(n_estimators=1000, learning_rate=1.0, max_depth=1, random_state=RANDOM)
    elif method == "K_Neighbors_Classifier":
        clf = neighbors.KNeighborsClassifier(n_neighbors=50, n_jobs=1)
    elif method == "Bernoulli_NB":
        clf = BernoulliNB()
    elif method == "AdaBoost_Classifier":
        clf = AdaBoostClassifier(n_estimators=1000, random_state=RANDOM)
    elif method == "Decision_Tree_Classifier":
        clf = tree.DecisionTreeClassifier(random_state=RANDOM)

    pred = clf.fit(X, y)
    w_path = open(save_path, "w")
    for f_a in all_feature_A:
        for f_b in all_feature_B:
            a = f_a.split("\t")
            b = f_b.split("\t")
            pair_name = a[0] + "\t" + b[0]
            feature_ab = a[1] + "," + b[1]
            feature_ab = np.array([list(map(int, feature_ab.split(",")))])
            pre_label = pred.predict(feature_ab)[[0]][0]
            if pre_label == 1:
                probas_ = pred.predict_proba(feature_ab)[[0]][0][1]
                print(pair_name + "\t" + str(probas_) + "\n")
                w_path.write(pair_name + "\t" + str(probas_) + "\n")
    w_path.close()


def choose_best_auc(auc_file):
    max = 0
    res = []
    with open(auc_file) as f:
        for line in f:
            tem = line.strip().split("\t")
            score = float(tem[1])
            if score > max:
                max = score
                res = tem
    return res


def run(best_train):
    # best_train = all_best_train[0]
    os.makedirs(res_dir, exist_ok=True)
    train_file_path = best_train[0]
    best_auc = float(best_train[1])
    method = best_train[2]

    # write best choose
    pairs = os.path.basename(train_file_path).split("-")[1].split(".")[0]
    best_choose_path = os.path.join(res_dir, pairs + ".best")
    best_choose_w = open(best_choose_path, "w")
    best_choose_w.write("\t".join(best_train))
    best_choose_w.close()

    best_predict_path = os.path.join(res_dir, pairs + ".result")

    data = pd.read_csv(train_data, sep=',', header=None)
    dataCol = data.columns.values.tolist()[1:]
    data_x = np.array(data[dataCol])
    data_y = np.array(data[0])
    if best_auc > 0.6:
        start(method=method,
              save_path=best_predict_path,
              X=data_x,
              y=data_y,
              pairs=pairs)


if __name__ == '__main__':
    all_trains = "/home/yangfang/Sash/QFO/test_ML/test_res/*.auc"
    trains_res = glob.glob(all_trains)
    train_data = "/home/yangfang/Sash/QFO/test_ML/9606_10090/9606_10090_test_k5.txt"
    feature_path = "/home/yangfang/Sash/QFO/QfO_release_2018_04/protein_all_trans_test/*.feature"
    res_dir = "/home/yangfang/Sash/QFO/test_ML/predict_result/"

    all_best_train = [choose_best_auc(var) for var in trains_res]
    pool = multiprocessing.Pool(processes=CORE)
    pool.map(run,all_best_train)