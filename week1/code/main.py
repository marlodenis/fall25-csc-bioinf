from dbg import DBG
from utils import read_data

import sys
from python import os


def get_n50(arr) -> int:
    mid_pt = sum(arr) / 2
    loc = 0
    for ctg in arr:
        curr = ctg
        loc += curr
        if loc >= mid_pt:
            return curr


def main():
    argv = sys.argv
    short1, short2, long1 = read_data(os.path.join('./', argv[1]))

    k = 25
    dbg = DBG(k=k, data_list=short1.append(short2).append(long1))
    ctg_info = []
    with open(os.path.join('./', argv[1], 'contig.fasta'), 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write('>contig_%d\n' % i)
            f.write(c + '\n')
            ctg_info.append(len(c))
    print(ctg_info)
    print("n50:", get_n50(ctg_info))

if __name__ == "__main__":
    main()