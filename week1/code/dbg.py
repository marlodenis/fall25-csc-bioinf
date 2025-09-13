import numpy as np


def reverse_complement(s: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[c] for c in reversed(s))

class Node:
    def __init__(self, kmer: str) -> None:
        self.kmer: str = kmer
        self.children: set[str] = set()
        self.count: int = 0
        self.visited: bool = False
        self.depth: int = 0
        self.max_depth_child: Optional[str] = None

    def add_child(self, kmer: str) -> None:
        self.children.add(kmer)

    def increase(self) -> None:
        self.count += 1

    def reset(self) -> None:
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self.count

    def get_children(self) -> list:
        return list(self.children)

    def remove_children(self, target: set) -> None:
        self.children -= target

class DBG:
    def __init__(self, k: int, data_list: list[str]) -> None:
        self.k: int = k
        self.nodes: dict[str, Node] = {}
        self.kmer2idx: dict[str, int] = {}
        self.kmer_count: int = 0
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list: list[str]) -> None:
        assert len(data_list) > 0
        assert self.k <= len(data_list[0])

    def _build(self, data_list: list[str]) -> None:
        for original in data_list:
            rc = reverse_complement(original)
            for i in range(len(original) - self.k - 1):
                self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def _add_node(self, kmer):
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[kmer] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[kmer].increase()
        return idx #### return kmer instead of idx???

    def _add_arc(self, kmer1, kmer2):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child):
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, idx: int) -> int:
        stack = [(idx, False)]
        depth_map = {}

        while stack:
            node_idx, processed = stack.pop()
            if processed:
                # Process children and calculate depth
                children = self._get_sorted_children(node_idx)
                max_depth, max_child = 0, None
                for child in children:
                    depth = depth_map.get(child, 0)
                    if depth > max_depth:
                        max_depth, max_child = depth, child
                self.nodes[node_idx].depth = max_depth + 1
                self.nodes[node_idx].max_depth_child = max_child
                depth_map[node_idx] = max_depth + 1
            else:
                self.nodes[node_idx].visited = True
                stack.append((node_idx, True))
                # Add children to stack
                for child in self._get_sorted_children(node_idx):
                    if not self.nodes[child].visited:
                        stack.append((child, False))

        return depth_map.get(idx, 0)

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
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path):
        for idx in path:
            del self.nodes[idx]
        path_set = set(path)
        for idx in self.nodes.keys():
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path):
        if len(path) < 1:
            return None
        concat = copy.copy(self.nodes[path[0]].kmer)
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
