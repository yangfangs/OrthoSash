#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~get neg pairs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import random

def get_protein_list(path):
    res = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                tem = line.strip()[1:]
                res.append(tem)
    return res

def neg_homologous(portein_A: str, protein_B: str, pos_filter: str,rd = 123) -> str:
    hsa_gene = get_protein_list(portein_A)
    mmu_gene = get_protein_list(protein_B)

    random.seed(rd)

    pos_pair = set()

    with open(pos_filter) as f:
        for line in f:
            tem = line.strip().split("\t")
            pos_pair.add(tem[1]+"\t"+tem[2])

    # with open("/home/yangfang/Sash/QFO/QfO_release_2018_04/othFinder_9606_272561.txt") as f:
    #     for line in f:
    #         tem = line.strip().split("\t")
    #         pos_pair.add(tem[1] + "\t" + tem[2])



    def get_random():
        a = random.choice(hsa_gene)
        m = random.choice(mmu_gene)
        pair1 = a + "\t" + m
        pair2 = m + "\t" + a
        return pair1,pair2
    w_name = pos_filter.split(".")[0] + "-neg" + ".txt"
    w_path = open(w_name,"w")
    set_res = set()
    while(len(set_res) != len(pos_pair)*10):
        # print(len(set_res))
        pair1, pair2 = get_random()
        # if pair1 not in pos_pair and pair2 not in pos_pair:
        #     if pair1 not in set_res and pair2 not in set_res:
        #         set_res.add(pair1)
        if pair1 not in pos_pair:
            if pair1 not in set_res:
                set_res.add(pair1)
    for line in set_res:
        w_path.write("0" + "\t" + line +"\n")
    w_path.close()
    return w_name

if __name__ == '__main__':
    protein_a =""
    protein_b =""
    pos_filter =""
    rd = 123
    neg_homologous(portein_A= protein_a,
                   protein_B= protein_b,
                   pos_filter= pos_filter,
                   rd=rd)