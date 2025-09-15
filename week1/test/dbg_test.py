from week1.code.dbg import reverse_complement
from week1.code.dbg import Node, DBG
import time

def test_performance():
    """Test assembly performance"""
    data = ["ATGC" * 100]  # 400bp sequence
    start = time.time()
    dbg = DBG(k=10, data_list=data)
    c = dbg.get_longest_contig()
    assert c is not None
    assert time.time() - start < 1.0  # Should complete in <1s

# Method tests
def test_complement():
    # arrange
    cases = ["", "A", "T", "C", "G", "ATGC", "AAAACCCGGT"]
    responses = ["", "T", "A", "G", "C", "GCAT", "ACCGGGTTTT"]
    # act
    results = [reverse_complement(c) for c in cases]
    # assert
    assert results == responses

# Node tests

def test_node_init():
    # arrange
    kmer = "ATGC"
    # act
    node = Node(kmer)
    # assert
    assert node.kmer == kmer
    assert node.get_count() == 0
    assert node.get_children() == []
    assert node.max_depth_child is None

def test_node_children():
    # arrange
    kmer = "AAAA"
    child1 = "AAAT"
    child2 = "AATA"

    # act
    parent = Node(kmer)

    # assert
    assert parent.kmer == kmer
    assert parent.get_children() == []

    # act
    parent.add_child(child1)

    # assert
    assert parent.kmer == kmer
    assert parent.get_children() == [child1]

    # act
    parent.add_child(child2)

    # assert
    assert parent.kmer == kmer
    assert parent.get_children() == [child1, child2]

    # act
    parent.remove_children({child1, child2})

    # assert
    assert parent.kmer == kmer
    assert parent.get_children() == []

# chatbot tests

from week1.code.dbg import reverse_complement, Node, DBG
import pytest


# Existing tests remain here...

# New DBG class tests
def test_dbg_init_valid():
    """Test DBG initialization with valid data"""
    data = ["ATGC", "GCAT"]
    dbg = DBG(k=2, data_list=data)
    assert dbg.k == 2
    assert len(dbg.nodes) > 0
    assert len(dbg.kmer2idx) > 0


def test_dbg_init_invalid_data():
    """Test DBG initialization with empty data"""
    with pytest.raises(AssertionError):
        DBG(k=2, data_list=[])


def test_dbg_init_invalid_k():
    """Test DBG initialization with k larger than read length"""
    with pytest.raises(AssertionError):
        DBG(k=10, data_list=["ATGC"])


def test_add_node():
    """Test node addition to DBG"""
    dbg = DBG(k=2, data_list=["ATGC"])

    # Add existing node
    idx1 = dbg._add_node("AT")
    assert dbg.nodes[idx1].count == 3  # Already added during init

    # Add new node
    idx2 = dbg._add_node("GC")
    assert dbg.kmer_count == 4
    assert dbg.nodes[idx2].kmer == "GC"


def test_add_arc():
    """Test arc addition to DBG"""
    dbg = DBG(k=2, data_list=["ATG"])
    dbg._reset()
    dbg._add_arc("AT", "TG")

    # Check nodes exist
    assert "AT" in dbg.kmer2idx
    assert "TG" in dbg.kmer2idx

    # Check child relationship
    at_idx = dbg.kmer2idx["AT"]
    assert "TG" in dbg.nodes[at_idx].children


def test_get_depth():
    """Test depth calculation in DBG"""
    dbg = DBG(k=2, data_list=["ATGC"])
    dbg._reset()

    # Linear path: AT -> TG -> GC
    at_idx = dbg.kmer2idx["AT"]
    tg_idx = dbg.kmer2idx["TG"]
    gc_idx = dbg.kmer2idx["GC"]

    # Manually set children for controlled test
    dbg.nodes[at_idx].children = ["TG"]
    dbg.nodes[tg_idx].children = ["GC"]

    # Calculate depths
    depth_gc = dbg._get_depth(gc_idx)  # Should be 1 (leaf)
    depth_tg = dbg._get_depth(tg_idx)  # Should be 2
    depth_at = dbg._get_depth(at_idx)  # Should be 3

    assert depth_gc == 1
    assert depth_tg == 2
    assert depth_at == 3


def test_get_longest_path():
    """Test longest path detection"""
    dbg = DBG(k=2, data_list=["ATGC"])
    dbg._reset()

    # Create a linear path: AT -> TG -> GC
    at_idx = dbg.kmer2idx["AT"]
    tg_idx = dbg.kmer2idx["TG"]
    gc_idx = dbg.kmer2idx["GC"]

    dbg.nodes[at_idx].children = ["TG"]
    dbg.nodes[tg_idx].children = ["GC"]

    path = dbg._get_longest_path()
    assert path == [at_idx, tg_idx, gc_idx]


def test_delete_path():
    """Test path deletion from DBG"""
    dbg = DBG(k=2, data_list=["ATGC"])
    dbg._reset()

    # Create a path: AT -> TG -> GC
    at_idx = dbg.kmer2idx["AT"]
    tg_idx = dbg.kmer2idx["TG"]
    gc_idx = dbg.kmer2idx["GC"]

    dbg.nodes[at_idx].children = ["TG"]
    dbg.nodes[tg_idx].children = ["GC"]

    # Delete path
    path = [at_idx, tg_idx, gc_idx]
    dbg._delete_path(path)

    # Verify nodes are deleted
    assert at_idx not in dbg.nodes
    assert tg_idx not in dbg.nodes
    assert gc_idx not in dbg.nodes


def test_concat_path():
    """Test path concatenation"""
    dbg = DBG(k=2, data_list=["ATGC"])
    dbg._reset()

    # Create a path: AT -> TG -> GC
    at_idx = dbg.kmer2idx["AT"]
    tg_idx = dbg.kmer2idx["TG"]
    gc_idx = dbg.kmer2idx["GC"]

    path = [at_idx, tg_idx, gc_idx]
    contig = dbg._concat_path(path)
    assert contig == "ATGC"


def test_get_longest_contig():
    """Test full contig extraction process"""
    dbg = DBG(k=2, data_list=["ATGC"])
    contig = dbg.get_longest_contig()

    # Should return the full sequence
    assert contig == "ATGC"

    # Verify nodes are deleted
    assert len(dbg.nodes) == 0


def test_reverse_complement_handling():
    """Test reverse complement handling in DBG"""
    data = ["ATGC"]
    dbg = DBG(k=2, data_list=data)

    # Should contain both original and reverse complement kmers
    assert "AT" in dbg.kmer2idx
    assert "TG" in dbg.kmer2idx  # Original
    assert "GC" in dbg.kmer2idx  # Reverse complement of AT
    assert "CA" in dbg.kmer2idx  # Reverse complement of GC


# def test_multiple_contigs():
#     """Test extraction of multiple contigs"""
#     data = ["ATGC", "GCTA"]  # Two separate sequences
#     dbg = DBG(k=2, data_list=data)
#
#     # Extract first contig
#     contig1 = dbg.get_longest_contig()
#     assert len(contig1) == 7
#
#     # Extract second contig
#     contig2 = dbg.get_longest_contig()
#     assert len(contig2) == 4
#
#     # No more contigs should be available
#     assert dbg.get_longest_contig() is None


def test_circular_sequence():
    """Test handling of circular sequences"""
    data = ["ATGCA"]  # Circular: AT->TG->GC->CA->AT
    dbg = DBG(k=2, data_list=data)

    # Should handle circularity without infinite loops
    contig = dbg.get_longest_contig()
    assert len(contig) >= 2  # At least one k-mer pair


def test_kmer_count_accuracy():
    """Test k-mer counting accuracy"""
    data = ["ATATAT"]  # Contains AT multiple times
    dbg = DBG(k=2, data_list=data)

    at_idx = dbg.kmer2idx["AT"]
    count = dbg.nodes[at_idx].get_count()

    # AT added 4 times in original and 4 times in reverse complement
    assert count == 8
