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
    assert parent.get_children().sort() == [child1, child2].sort()

    # act
    parent.remove_children({child1, child2})

    # assert
    assert parent.kmer == kmer
    assert parent.get_children() == []


def test_reverse_complement():
    # Test single nucleotides
    assert reverse_complement("A") == "T"
    assert reverse_complement("T") == "A"
    assert reverse_complement("G") == "C"
    assert reverse_complement("C") == "G"

    # Test dinucleotides
    assert reverse_complement("AT") == "AT"
    assert reverse_complement("CG") == "CG"
    assert reverse_complement("AG") == "CT"
    assert reverse_complement("TC") == "GA"

    # Test longer sequences
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("GGATCC") == "GGATCC"
    assert reverse_complement("ACGTACGT") == "ACGTACGT"
    print("✓ test_reverse_complement passed")


def test_node():
    # Test initialization
    node = Node("ATCG")
    assert node.kmer == "ATCG"
    assert node.children == set()
    assert node.count == 0
    assert not node.visited
    assert node.depth == 0
    assert node.max_depth_child is None

    # Test add_child
    node.add_child(1)
    assert node.children == {1}
    node.add_child(2)
    assert node.children == {1, 2}
    node.add_child(1)  # Adding duplicate
    assert node.children == {1, 2}

    # Test increase
    assert node.count == 0
    node.increase()
    assert node.count == 1
    node.increase()
    assert node.count == 2

    # Test reset
    node.visited = True
    node.depth = 5
    node.max_depth_child = 3
    node.reset()
    assert not node.visited
    assert node.depth == 0
    assert node.max_depth_child is None
    # Ensure children and count are not reset
    assert node.children == {1, 2}
    assert node.count == 2

    # Test get_count
    assert node.get_count() == 2

    # Test get_children
    children = node.get_children()
    assert 1 in children
    assert 2 in children
    assert len(children) == 2

    # Test remove_children
    node.remove_children({2, 4})
    assert node.children == {1}
    print("✓ test_node passed")


def test_dbg_init():
    # Test valid data
    data = ["ATCG", "GCTA"]
    dbg = DBG(2, data)
    assert dbg.k == 2
    assert len(dbg.nodes) == 8  # 4 kmers per sequence
    assert dbg.kmer_count == 8

    # Test invalid data
    try:
        DBG(3, [])  # Empty data list
        assert False, "Expected exception for empty data list"
    except Exception:
        pass

    try:
        DBG(10, ["ATCG"])  # k > sequence length
        assert False, "Expected exception for k > sequence length"
    except Exception:
        pass
    print("✓ test_dbg_init passed")


def test_dbg_add_node():
    dbg = DBG(2, ["ATCG"])
    # First addition

    idx1 = dbg._add_node("AT")
    assert idx1 == 0
    assert dbg.nodes[0].kmer == "AT"
    assert dbg.nodes[0].count == 2
    # Adding existing node
    assert dbg.nodes[1].count == 1

    idx2 = dbg._add_node("AT")
    assert idx2 == 0
    assert dbg.nodes[0].count == 3
    # Adding new node
    idx3 = dbg._add_node("TC")
    assert idx3 == 1
    assert dbg.nodes[1].kmer == "TC"
    assert dbg.nodes[1].count == 2
    print("✓ test_dbg_add_node passed")


def test_dbg_add_arc():
    dbg = DBG(2, ["ATCG"])
    # Add arc AT->TC
    dbg._add_arc("AT", "TC")
    # Check nodes were created
    assert 0 in dbg.nodes
    assert 1 in dbg.nodes
    assert dbg.nodes[0].kmer == "AT"
    assert dbg.nodes[1].kmer == "TC"
    # Check child relationship
    assert 1 in dbg.nodes[0].children
    print("✓ test_dbg_add_arc passed")


def test_dbg_get_sorted_children():
    # Create DBG with minimal data and reset it
    dbg = DBG(2, ["AT"])
    dbg.nodes = {}
    dbg.kmer2idx = {}
    dbg.kmer_count = 0

    # Now add nodes manually
    node1 = dbg.nodes[dbg._add_node("AA")]
    node2 = dbg.nodes[dbg._add_node("TT")]
    node3 = dbg.nodes[dbg._add_node("CC")]

    node1.increase()  # count=1
    node2.increase()  # count=1
    node2.increase()  # count=2
    node3.increase()  # count=1
    node3.increase()  # count=2
    node3.increase()  # count=3

    # Create a parent node and add children
    parent_idx = dbg._add_node("AT")
    parent = dbg.nodes[parent_idx]
    parent.add_child(dbg.kmer2idx["AA"])
    parent.add_child(dbg.kmer2idx["TT"])
    parent.add_child(dbg.kmer2idx["CC"])

    # Get sorted children (should be sorted by count descending)
    sorted_children = dbg._get_sorted_children(parent_idx)
    expected_order = [
        dbg.kmer2idx["CC"],  # count=3
        dbg.kmer2idx["TT"],  # count=2
        dbg.kmer2idx["AA"]  # count=1
    ]
    assert sorted_children == expected_order
    print("✓ test_dbg_get_sorted_children passed")

def test_dbg_get_depth():
    dbg = DBG(2, ["ATCG"])
    # Create a linear chain: A -> B -> C
    a_idx = dbg._add_node("AA")
    b_idx = dbg._add_node("BB")
    c_idx = dbg._add_node("CC")

    dbg.nodes[a_idx].add_child(b_idx)
    dbg.nodes[b_idx].add_child(c_idx)

    # Calculate depths
    depth_a = dbg._get_depth(a_idx)
    depth_b = dbg._get_depth(b_idx)
    depth_c = dbg._get_depth(c_idx)

    # C has no children -> depth=1
    assert depth_c == 1
    # B has C as child -> depth=1 + depth(C)=2
    assert depth_b == 2
    # A has B as child -> depth=1 + depth(B)=3
    assert depth_a == 3

    # Check max_depth_child pointers
    assert dbg.nodes[a_idx].max_depth_child == b_idx
    assert dbg.nodes[b_idx].max_depth_child == c_idx
    assert dbg.nodes[c_idx].max_depth_child is None
    print("✓ test_dbg_get_depth passed")


def test_dbg_reset():
    dbg = DBG(2, ["ATCG"])
    # Modify some nodes
    for node in dbg.nodes.values():
        node.visited = True
        node.depth = 5
        node.max_depth_child = 1

    dbg._reset()

    # Check all nodes are reset
    for node in dbg.nodes.values():
        assert not node.visited
        assert node.depth == 0
        assert node.max_depth_child is None
    print("✓ test_dbg_reset passed")


def test_dbg_get_longest_path():
    dbg = DBG(2, ["ATCG"])
    # Create two paths: A->B->C (length 3) and D->E (length 2)
    a_idx = dbg._add_node("AA")
    b_idx = dbg._add_node("BB")
    c_idx = dbg._add_node("CC")
    d_idx = dbg._add_node("DD")
    e_idx = dbg._add_node("EE")

    dbg.nodes[a_idx].add_child(b_idx)
    dbg.nodes[b_idx].add_child(c_idx)
    dbg.nodes[d_idx].add_child(e_idx)

    # Get longest path
    path = dbg._get_longest_path()
    expected_path = [a_idx, b_idx, c_idx]
    assert path == expected_path
    print("✓ test_dbg_get_longest_path passed")


def test_dbg_delete_path():
    dbg = DBG(2, ["ATCG"])
    # Create nodes and relationships
    a_idx = dbg._add_node("AA")
    b_idx = dbg._add_node("BB")
    c_idx = dbg._add_node("CC")
    d_idx = dbg._add_node("DD")

    dbg.nodes[a_idx].add_child(b_idx)
    dbg.nodes[b_idx].add_child(c_idx)
    dbg.nodes[d_idx].add_child(b_idx)  # D also points to B

    # Delete path A->B->C
    path = [a_idx, b_idx, c_idx]
    dbg._delete_path(path)

    # Check nodes are deleted
    assert a_idx not in dbg.nodes
    assert b_idx not in dbg.nodes
    assert c_idx not in dbg.nodes
    assert d_idx in dbg.nodes

    # Check D's children are updated
    assert dbg.nodes[d_idx].children == set()
    print("✓ test_dbg_delete_path passed")


def test_dbg_concat_path():
    dbg = DBG(2, ["ATCG"])
    # Create nodes with kmers
    a_idx = dbg._add_node("AT")
    b_idx = dbg._add_node("TC")
    c_idx = dbg._add_node("CG")

    path = [a_idx, b_idx, c_idx]
    contig = dbg._concat_path(path)
    assert contig == "ATCG"

    # Test empty path
    assert dbg._concat_path([]) is None
    print("✓ test_dbg_concat_path passed")


def test_dbg_get_longest_contig(): ###
    # Test with simple linear sequence
    dbg = DBG(3, ["ATCG"])
    contig = dbg.get_longest_contig()
    assert contig in ["ATCG", "CGAT"]  # Could be original or reverse complement

    # Test with branching
    dbg = DBG(2, ["ATCG", "ATGG"])  # Creates branches at AT
    contig1 = dbg.get_longest_contig()
    # Should return either ATCG or ATGG (both length 4)
    assert contig1 in ["ATCG", "ATGG"]

    # After first contig is removed, next call should return the other
    contig2 = dbg.get_longest_contig()
    assert contig2 in ["ATCG", "ATGG"]
    assert contig1 != contig2

    # Test with reverse complement
    dbg = DBG(3, ["ATCG"])  # Reverse complement is CGAT
    contig = dbg.get_longest_contig()
    assert contig in ["ATCG", "CGAT"]
    print("✓ test_dbg_get_longest_contig passed")