from dbg import DBG
from utils import read_data

import sys
import time


def get_n50(arr) -> int:
    mid_pt = sum(arr) / 2
    loc = 0
    for ctg in arr:
        curr = ctg
        loc += curr
        if loc >= mid_pt:
            return curr


def main():
    start_time = time.time()

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
            f.write('>contig_'+ str(i) +'\n')
            f.write(c + '\n')
            ctg_info.append(len(c))
    total_time = time.time() - start_time
    print(f"{total_time:.2f}    {get_n50(ctg_info)}")

if __name__ == "__main__":
    main()