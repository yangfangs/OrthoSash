import os

import numpy as np
import pandas as pd
from scipy import interp
import matplotlib.pyplot as plt

from sklearn import svm, datasets, neighbors, tree
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold

# #############################################################################
from sklearn.naive_bayes import BernoulliNB

# #############################################################################
# Classification and ROC analysis
from PDS.QFO.test_all_server.parameter import RANDOM


def start(method, save_path, X, y, num):
    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    cv = StratifiedKFold(n_splits=10, random_state=RANDOM,shuffle=True)

    plt.figure(figsize=(20, 10))

    if method == "SVC":
        clf = svm.SVC(probability=True, gamma='auto', random_state=RANDOM,max_iter=10000)
    elif method == "Random_Forest_Classifier":
        clf = RandomForestClassifier(n_estimators=1000, max_features="sqrt", random_state=RANDOM, n_jobs=1)
    elif method == "Gradient_Boosting_Classifier":
        clf = GradientBoostingClassifier(n_estimators=1000, learning_rate=1.0, max_depth=1, random_state=RANDOM)
    elif method == "K_Neighbors_Classifier":
        clf = neighbors.KNeighborsClassifier(n_neighbors= 50, n_jobs=1)
    elif method == "Bernoulli_NB":
        clf = BernoulliNB()
    elif method == "AdaBoost_Classifier":
        clf = AdaBoostClassifier(n_estimators=1000, random_state=RANDOM)
    elif method == "Decision_Tree_Classifier":
        clf = tree.DecisionTreeClassifier(random_state=RANDOM)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    i = 0
    for train, test in cv.split(X, y):
        probas_ = clf.fit(X[train], y[train]).predict_proba(X[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,
                 label='ROC fold %d (AUC = %0.3f)' % (i, roc_auc))

        i += 1
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='k',
             label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.3f $\pm$ %0.3f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    # plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
    #                  label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('10-fold cross-validation ({0})'.format(method))
    plt.legend(loc="lower right")
    # plt.show()
    w_value = save_path + num + "-" + method + ".value"
    df = pd.DataFrame({"tpr": mean_tpr, "fpr": mean_fpr})
    df.to_csv(w_value, header=True, index=None)
    w_path = save_path + num + "-" + method + ".pdf"
    plt.savefig(w_path)
    plt.close()
    return mean_auc


if __name__ == '__main__':
    train_data = "/home/yangfang/Sash/QFO/test_ML/9606_10090/9606_10090_test_k5.txt"
    res_dir = "/home/yangfang/Sash/QFO/test_ML/9606_10090/result/"
    os.makedirs(res_dir, exist_ok=True)
    # method =["SVC","Random Forest Classifier",
    #          "Gradient Boosting Classifier",
    #          "K Neighbors Classifier",
    #          "Bernoulli NB",
    #          "AdaBoost Classifier",
    #          "Decision Tree Classifier"]
    # method =["Random Forest Classifier"]
    method = [
        "Gradient Boosting Classifier",
        "K Neighbors Classifier",
        "Bernoulli NB",
        "AdaBoost Classifier",
        "Decision Tree Classifier"]
    for line in method:
        print(line)
        res = open(res_dir + line + "_ROC.txt", "w")
        all_auc = []

        data = pd.read_csv(train_data, sep=',', header=None)
        dataCol = data.columns.values.tolist()[1:]
        data_x = np.array(data[dataCol])
        data_y = np.array(data[0])
        w_path = res_dir
        mean_acu = start(line, w_path, data_x, data_y, line)
        all_auc.append(str(mean_acu)[:5])
        res.write(line + "," + ",".join(all_auc) + "\n")
        res.close()
