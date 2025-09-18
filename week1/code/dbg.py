from typing import Optional
import copy

def reverse_complement(key: str):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}

    key = list(key[::-1])
    for i in range(len(key)):
        key[i] = complement[key[i]]
    return ''.join(key)

class Node:
    kmer: str
    children: set[int]
    count: int
    visited: bool
    depth: int
    max_depth_child: Optional[int]
    children_sorted : bool

    def __init__(self, kmer: str):
        self.kmer = kmer
        self.children = set()
        self.count = 0
        self.visited = False
        self.depth = 0
        self.max_depth_child = None
        self.children_sorted = False

    def add_child(self, kmer_idx: int) -> None:
        self.children.add(kmer_idx)

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
        self.children = set([itm for itm in self.children if itm not in target])
        # self.children = self.children - target


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

    def _check(self, data_list: list[str]):
        # check data list
        try:
            assert len(data_list) > 0
            assert self.k <= len(data_list[0])
        except Exception as e:
            print(f"Error in data_list or k: {e}")
            raise e

    def _build(self, data_list: list[str]):
        for original in data_list:
            rc = reverse_complement(original)
            for i in range(len(original) - self.k):
                self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def _add_node(self, kmer: str) -> int:
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1: str, kmer2: str):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child: int):
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)
        return children

    def _get_depth(self, start_idx: int) -> int:
        """
        Ultra-optimized iterative depth calculation for very large, deep graphs.
        Minimizes all possible overhead for maximum performance.
        """
        # Early return if already visited
        start_node = self.nodes[start_idx]
        if start_node.visited:
            return start_node.depth

        # Pre-allocate stack with estimated size to avoid reallocations
        stack: list[tuple[int, int]] = []
        stack.append((start_idx, 0))

        # Single allocation for children lists - reuse for better memory performance
        temp_children: list[int] = []

        while stack:
            current_idx, phase = stack.pop()
            current_node = self.nodes[current_idx]

            if phase == 0:
                # First visit
                if current_node.visited:
                    continue

                current_node.visited = True

                # Clear and reuse temp list instead of creating new ones
                temp_children.clear()
                temp_children.extend(current_node.children)

                if not temp_children:
                    # Leaf node - immediate computation
                    current_node.depth = 1
                    current_node.max_depth_child = None
                    continue

                # Sort in-place for memory efficiency
                temp_children.sort(key=lambda child_idx: self.nodes[child_idx].count, reverse=True)

                # Add back for post-processing
                stack.append((current_idx, 1))

                # Add children in reverse order, only unvisited ones
                for i in range(len(temp_children) - 1, -1, -1):
                    child_idx = temp_children[i]
                    if not self.nodes[child_idx].visited:
                        stack.append((child_idx, 0))

            else:  # phase == 1
                # Post-process: compute depth from children
                # Reuse temp_children list
                temp_children.clear()
                temp_children.extend(current_node.children)
                temp_children.sort(key=lambda child_idx: self.nodes[child_idx].count, reverse=True)

                max_depth = 0
                max_child: Optional[int] = None

                # Optimized loop with direct access
                for child_idx in temp_children:
                    child_depth = self.nodes[child_idx].depth
                    if child_depth > max_depth:
                        max_depth = child_depth
                        max_child = child_idx

                current_node.depth = max_depth + 1
                current_node.max_depth_child = max_child

        return start_node.depth

    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()

    def _get_longest_path(self):
        """
        Optimized longest path that avoids redundant work by skipping visited nodes.
        """
        max_depth = 0
        max_idx: Optional[int] = None

        for idx in self.nodes.keys():
            # Skip if already visited
            if self.nodes[idx].visited:
                # Still need to check if this is maximum depth
                if self.nodes[idx].depth > max_depth:
                    max_depth = self.nodes[idx].depth
                    max_idx = idx
                continue

            # Calculate depth for unvisited node
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth = depth
                max_idx = idx

        # Build path
        path: list[int] = []
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
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
