from dbg import DBG
from utils import read_data

import sys


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
    data_list = read_data(argv[1])

    k = 25
    dbg = DBG(k=k, data_list=data_list)
    ctg_info = []
    out_path = argv[1] + '/' + 'contig.fasta'
    with open(out_path, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write('>contig_'+ str(i) +'\n')
            f.write(c + '\n')
            ctg_info.append(len(c))
    print(ctg_info)
    print("n50:", get_n50(ctg_info))

if __name__ == "__main__":
    main()