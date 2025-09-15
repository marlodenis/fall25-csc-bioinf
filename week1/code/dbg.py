from typing import Optional
# from python import copy

def reverse_complement(s: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    try:
        return ''.join(complement[c] for c in reversed(s))
    except KeyError as e:
        print(f"Invalid character in sequence: {e}")
        return ""

class Node:
    kmer: str
    children: list[str]
    count: int
    visited: bool
    depth: int
    max_depth_child: Optional[str]

    def __init__(self, kmer: str):
        self.kmer = kmer
        self.children = []
        self.count = 0 ### Was 0, changed to 1
        self.visited = False
        self.depth = 0 ### TODO: investigate
        self.max_depth_child = None

    def add_child(self, kmer: str) -> None:
        self.children.append(kmer)

    def increase(self) -> None:
        self.count += 1

    def reset(self) -> None:
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self.count

    def get_children(self) -> list:
        return self.children

    def remove_children(self, target: list) -> None:
        self.children = [itm for itm in self.children if itm not in target]


class DBG:
    k: int
    nodes: dict[int, Node]
    kmer2idx: dict[str, int]
    kmer_count: int
    def __init__(self, k, data_list):
        self.k = k
        self.nodes = {}
        # private
        self.kmer2idx = {}
        self.kmer_count = 0
        # build
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        # check data list
        try:
            assert len(data_list) > 0
            assert self.k <= len(data_list[0])
        except Exception as e:
            print(f"Error in data_list or k: {e}")
            raise e

    def _build(self, data_list):
        for original in data_list:
            rc = reverse_complement(original)
            for i in range(len(original) - self.k): ###
                self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def show_count_distribution(self):
        count = [0] * 30
        for idx in self.nodes:
            count[self.nodes[idx].get_count()] += 1
        print(count[0:10])
        # plt.plot(count)
        # plt.show()

    def _add_node(self, kmer):
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1, kmer2):
        idx1 = self._add_node(kmer1)
        self._add_node(kmer2)
        self.nodes[idx1].add_child(kmer2)
        # idx2 = self._add_node(kmer2)
        # self.nodes[idx1].add_child(idx2)

    def _get_count(self, child):
        return self.nodes[self.kmer2idx[child]].get_count()

    def _get_sorted_children(self, idx):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx):
        if not self.nodes[idx].visited:
            self.nodes[idx].visited = True
            children = self._get_sorted_children(idx)
            max_depth, max_child = 0, None
            for child in children:
                depth = self._get_depth(self.kmer2idx[child])
                if depth > max_depth:
                    max_depth, max_child = depth, child
            self.nodes[idx].depth, self.nodes[idx].max_depth_child = max_depth + 1, max_child
        return self.nodes[idx].depth

    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()

    def _get_longest_path(self):
        max_depth, max_idx = 0, None
        for idx in self.nodes.keys():
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth, max_idx = depth, idx

        path = []
        while max_idx is not None:
            path.append(max_idx)
            if self.nodes[max_idx].max_depth_child is None:
                max_idx = None
            else:
                 max_idx = self.kmer2idx[self.nodes[max_idx].max_depth_child]
        return path

    def _delete_path(self, path):
        kmers = [self.nodes[idx].kmer for idx in path]
        for idx in path:
            del self.nodes[idx]
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(kmers)

    def _concat_path(self, path):
        if len(path) < 1:
            return None
        # concat = copy.copy(self.nodes[path[0]].kmer)
        concat = self.nodes[path[0]].kmer
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self):
        # reset params in nodes for getting longest path
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
