import collections
import mmh3
from Bio import SeqIO
import glob

from PDS.QFO.test_all_server.parameter import KMER_K, HASH_SIZE


def kmers(seq, k):
    for i in range(len(seq) - k + 1):
        yield seq[i:i + k]


def simhash(kmers,d,sim_size = 64):
    total = sum(d.values())
    v = [0] * sim_size
    if sim_size == 64:
        for t in kmers:
            for i in range(sim_size):
                bitmask = 1 << i
                if mmh3.hash64(t, signed=False)[0] & bitmask:
                    v[i] += 1
                    # v[i] += d[t]/total
                else:
                    # v[i] = d[t] / total
                    v[i] -=1
    elif sim_size == 128:
        for t in kmers:
            for i in range(sim_size):
                bitmask = 1 << i
                if mmh3.hash128(t, signed=False) & bitmask:
                    v[i] += 1
                    # v[i] += d[t]/total
                else:
                    # v[i] = d[t] / total
                    v[i] -=1
    fingerprint = 0
    each_v = [0] *sim_size
    for i in range(sim_size):
        if v[i] > 0:
            each_v[i] = 1
            fingerprint += 1 << i
    return fingerprint, each_v

def get_kmers(seq, k=11):
    kmers_list = []
    d = collections.defaultdict(int)
    for km in kmers(seq, k):
        d[km] += 1
        kmers_list.append(km)
    return d, kmers_list


def get_feature(genome,write_path,kmer_K=5,min_size=128):
    write_file = open(write_path, "w")
    seq_db = {}
    reads = SeqIO.parse(genome, 'fasta')
    for read in reads:
        seq = str(read.seq)
        seq_name = str(read.name)
        seq_db[seq_name] = seq

    for k,v in seq_db.items():
        d,km = get_kmers(v,kmer_K)

        f, bi = simhash(km,d,min_size)
        write_file.write(k + "\t" + ",".join(map(str,bi)) + "\n")
    write_file.close()

if __name__ == '__main__':
    kmer_k = KMER_K
    min_size = HASH_SIZE
    genome = "/home/yangfang/Sash/QFO/QfO_release_2018_04/protein_all_trans_test2/*.fasta"
    files = glob.glob(genome)
    for line in files:
        print(line)
        name = line.split(".")
        write_path = name[0] + "_k" + str(kmer_k) + "_sim" + str(min_size) + ".feature"
        get_feature(line,write_path,kmer_k,min_size)