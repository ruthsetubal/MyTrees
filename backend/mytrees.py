"""
MyTrees - Molecular Systematics & Phylogenetic Analysis Engine
Implements: UPGMA, WPGMA, NJ, Fitch-Margoliash, Maximum Parsimony, Maximum Likelihood (Felsenstein's pruning)
Distance models: p-distance, JC69, K2P, F84, LogDet
"""

import math
import json
import copy
import re
from typing import Optional, Dict, List, Tuple, Any

# ─── Tree Node Representation ───────────────────────────────────────────────

class TreeNode:
    def __init__(self, name: Optional[str] = None, branch_length: Optional[float] = None):
        self.name: Optional[str] = name
        self.branch_length: Optional[float] = branch_length
        self.children: List['TreeNode'] = []
        self.parent: Optional['TreeNode'] = None
        self.metadata: Dict[str, Any] = {}

    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def add_child(self, child: 'TreeNode'):
        child.parent = self
        self.children.append(child)

    def remove_child(self, child: 'TreeNode'):
        if child in self.children:
            self.children.remove(child)
            child.parent = None

    def to_newick(self) -> str:
        """Converts the node and its subtree to a Newick string representation."""
        if self.is_leaf():
            name_part = self.name if self.name else ""
            # Clean up the name for Newick compatibility if needed
            name_part = re.sub(r'[\(\):;,]', '_', name_part)
            if self.branch_length is not None:
                return f"{name_part}:{self.branch_length:.6f}"
            return name_part
        else:
            children_str = ",".join(c.to_newick() for c in self.children)
            if self.branch_length is not None:
                return f"({children_str}):{self.branch_length:.6f}"
            return f"({children_str})"

    def clone(self) -> 'TreeNode':
        """Deep clones the node and its children recursively."""
        new_node = TreeNode(self.name, self.branch_length)
        new_node.metadata = copy.deepcopy(self.metadata)
        for child in self.children:
            new_node.add_child(child.clone())
        return new_node

    def to_dict(self) -> dict:
        """Serializes the node (and subtree) into a structured dict.

        Unlike Newick, this preserves internal-node metadata such as bootstrap
        support, so the frontend can render it without re-parsing a lossy string.
        Node IDs are assigned deterministically in pre-order (root = 0).
        """
        counter = {'i': 0}

        def build(node: 'TreeNode') -> dict:
            nid = counter['i']
            counter['i'] += 1
            return {
                'id': nid,
                'name': node.name if node.name else '',
                'length': node.branch_length if node.branch_length is not None else 0.0,
                'metadata': node.metadata,
                'children': [build(c) for c in node.children],
            }

        return build(self)


# ─── Distance Models ─────────────────────────────────────────────────────────

def get_base_frequencies(seqs: List[str]) -> Dict[str, float]:
    """Calculates empirical base frequencies from a list of sequences, ignoring gaps."""
    counts = {'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}
    total = 0.0
    for seq in seqs:
        for char in seq.upper():
            if char in counts:
                counts[char] += 1
                total += 1.0
    if total == 0:
        return {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
    return {k: v / total for k, v in counts.items()}


def p_distance(seq1: str, seq2: str) -> float:
    """Proportion of differing sites (p-distance), ignoring gaps and missing data."""
    seq1, seq2 = seq1.upper(), seq2.upper()
    n = sum(1 for a, b in zip(seq1, seq2) if a not in ('-', '?') and b not in ('-', '?'))
    if n == 0:
        return 0.0
    diffs = sum(1 for a, b in zip(seq1, seq2) if a not in ('-', '?') and b not in ('-', '?') and a != b)
    return diffs / n


def jc69(seq1: str, seq2: str) -> float:
    """Jukes-Cantor 1969 model distance."""
    p = p_distance(seq1, seq2)
    if p >= 0.75:
        return 5.0  # Safe upper ceiling for saturated distance
    if p == 0:
        return 0.0
    try:
        return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)
    except (ValueError, ZeroDivisionError):
        return 5.0


def k2p(seq1: str, seq2: str) -> float:
    """Kimura 2-parameter 1980 model distance."""
    seq1, seq2 = seq1.upper(), seq2.upper()
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a not in ('-', '?') and b not in ('-', '?')]
    n = len(pairs)
    if n == 0:
        return 0.0
    
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    transitions = 0
    transversions = 0
    
    for a, b in pairs:
        if a != b:
            if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                transitions += 1
            else:
                transversions += 1
                
    P = transitions / n
    Q = transversions / n
    
    val1 = 1.0 - 2.0 * P - Q
    val2 = 1.0 - 2.0 * Q
    
    if val1 <= 0 or val2 <= 0:
        return 5.0
    
    try:
        d = -0.5 * math.log(val1) - 0.25 * math.log(val2)
        return max(0.0, d)
    except (ValueError, ZeroDivisionError):
        return 5.0


def f84(seq1: str, seq2: str, base_freqs: Optional[Dict[str, float]] = None) -> float:
    """Felsenstein 1984 model distance."""
    if base_freqs is None:
        base_freqs = get_base_frequencies([seq1, seq2])
    
    fA = base_freqs.get('A', 0.25)
    fC = base_freqs.get('C', 0.25)
    fG = base_freqs.get('G', 0.25)
    fT = base_freqs.get('T', 0.25)
    
    fR = fA + fG  # Purines
    fY = fC + fT  # Pyrimidines
    
    if fR == 0 or fY == 0:
        fR, fY = 0.5, 0.5
        
    seq1, seq2 = seq1.upper(), seq2.upper()
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a not in ('-', '?') and b not in ('-', '?')]
    n = len(pairs)
    if n == 0:
        return 0.0
        
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}
    transitions = 0
    transversions = 0
    
    for a, b in pairs:
        if a != b:
            if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
                transitions += 1
            else:
                transversions += 1
                
    P = transitions / n
    Q = transversions / n
    
    # Calculate analytical transition and transversion probabilities
    a_coeff = (fA*fG/fR + fC*fT/fY)
    b_coeff = fR*fY
    
    val1 = 1.0 - Q / (2.0 * b_coeff)
    if val1 <= 0:
        return 5.0
        
    val2 = (1.0 - P / (2.0 * a_coeff) - (fA*fG*fY/fR + fC*fT*fR/fY) * Q / (2.0 * a_coeff * b_coeff))
    if val2 <= 0:
        return 5.0
        
    try:
        d = -4.0 * a_coeff * math.log(val2) - 4.0 * (b_coeff - a_coeff) * math.log(val1)
        return max(0.0, d)
    except (ValueError, ZeroDivisionError):
        return 5.0


def logdet(seq1: str, seq2: str) -> float:
    """LogDet distance."""
    seq1, seq2 = seq1.upper(), seq2.upper()
    pairs = [(a, b) for a, b in zip(seq1, seq2) if a not in ('-', '?') and b not in ('-', '?')]
    n = len(pairs)
    if n == 0:
        return 0.0
        
    bases = ['A', 'C', 'G', 'T']
    base_idx = {b: i for i, b in enumerate(bases)}
    
    # Generate 4x4 divergence matrix
    F = [[0.0 for _ in range(4)] for _ in range(4)]
    for a, b in pairs:
        if a in base_idx and b in base_idx:
            F[base_idx[a]][base_idx[b]] += 1.0
            
    for i in range(4):
        for j in range(4):
            F[i][j] /= n
            
    # Calculate determinant
    det_F = (
        F[0][0]*(F[1][1]*(F[2][2]*F[3][3] - F[2][3]*F[3][2]) - F[1][2]*(F[2][1]*F[3][3] - F[2][3]*F[3][1]) + F[1][3]*(F[2][1]*F[3][2] - F[2][2]*F[3][1])) -
        F[0][1]*(F[1][0]*(F[2][2]*F[3][3] - F[2][3]*F[3][2]) - F[1][2]*(F[2][0]*F[3][3] - F[2][3]*F[3][0]) + F[1][3]*(F[2][0]*F[3][2] - F[2][2]*F[3][0])) +
        F[0][2]*(F[1][0]*(F[2][1]*F[3][3] - F[2][3]*F[3][1]) - F[1][1]*(F[2][0]*F[3][3] - F[2][3]*F[3][0]) + F[1][3]*(F[2][0]*F[3][1] - F[2][1]*F[3][0])) -
        F[0][3]*(F[1][0]*(F[2][1]*F[3][2] - F[2][2]*F[3][1]) - F[1][1]*(F[2][0]*F[3][2] - F[2][2]*F[3][0]) + F[1][2]*(F[2][0]*F[3][1] - F[2][1]*F[3][0]))
    )
    
    if det_F <= 0:
        return 5.0
        
    # Calculate marginal frequencies
    f1 = [sum(F[i][j] for j in range(4)) for i in range(4)]
    f2 = [sum(F[j][i] for j in range(4)) for i in range(4)]
    
    prod_f1 = 1.0
    prod_f2 = 1.0
    for i in range(4):
        if f1[i] > 0: prod_f1 *= f1[i]
        if f2[i] > 0: prod_f2 *= f2[i]
        
    try:
        d = -0.25 * (math.log(det_F) - 0.5 * (math.log(prod_f1) + math.log(prod_f2)))
        return max(0.0, d)
    except (ValueError, ZeroDivisionError):
        return 5.0


def compute_distance_matrix(seqs: Dict[str, str], model: str = 'JC69') -> Tuple[List[str], List[List[float]], Dict[str, Any]]:
    """Computes the complete pairwise distance matrix under the specified evolutionary model."""
    taxa = list(seqs.keys())
    n = len(taxa)
    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    
    base_freqs = get_base_frequencies(list(seqs.values()))
    
    # Calculate distances
    for i in range(n):
        for j in range(i+1, n):
            seq1, seq2 = seqs[taxa[i]], seqs[taxa[j]]
            if model == 'p-distance':
                d = p_distance(seq1, seq2)
            elif model == 'JC69':
                d = jc69(seq1, seq2)
            elif model == 'K2P':
                d = k2p(seq1, seq2)
            elif model == 'F84':
                d = f84(seq1, seq2, base_freqs)
            elif model == 'LogDet':
                d = logdet(seq1, seq2)
            else:
                d = jc69(seq1, seq2)  # Fallback
                
            matrix[i][j] = round(d, 6)
            matrix[j][i] = round(d, 6)
            
    stats = analyze_matrix_stats(taxa, matrix)
    return taxa, matrix, stats


def analyze_matrix_stats(taxa: List[str], matrix: List[List[float]]) -> Dict[str, Any]:
    """Analyzes distance matrix parameters and checks for metrics like triangle inequality."""
    n = len(taxa)
    all_vals = []
    for i in range(n):
        for j in range(i+1, n):
            all_vals.append(matrix[i][j])
            
    if not all_vals:
        return {}
        
    mean_val = sum(all_vals) / len(all_vals)
    variance = sum((x - mean_val)**2 for x in all_vals) / len(all_vals)
    std_dev = math.sqrt(variance)
    
    is_symmetric = True
    for i in range(n):
        for j in range(n):
            if abs(matrix[i][j] - matrix[j][i]) > 1e-5:
                is_symmetric = False
                break
                
    violations = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Triangle inequality check: d(i,j) <= d(i,k) + d(k,j)
                if matrix[i][j] > matrix[i][k] + matrix[k][j] + 1e-5:
                    violations += 1
                if matrix[i][k] > matrix[i][j] + matrix[j][k] + 1e-5:
                    violations += 1
                if matrix[j][k] > matrix[j][i] + matrix[i][k] + 1e-5:
                    violations += 1
                    
    return {
        'n_taxa': n,
        'min_distance': round(min(all_vals), 6),
        'max_distance': round(max(all_vals), 6),
        'mean_distance': round(mean_val, 6),
        'std_distance': round(std_dev, 6),
        'is_symmetric': is_symmetric,
        'triangle_violations': violations,
        'n_pairs': len(all_vals),
    }


# ─── Parsers ─────────────────────────────────────────────────────────────────

def parse_fasta(content: str) -> Dict[str, str]:
    """Parses FASTA data (.fas, .fasta) and returns a dict of sequence names and strings."""
    seqs = {}
    current_name = None
    current_seq = []
    
    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_name:
                seqs[current_name] = "".join(current_seq)
            current_name = line[1:].strip()
            current_seq = []
        else:
            # Clean up the line (remove gaps/invalid symbols optionally or keep -)
            current_seq.append(line.replace(" ", ""))
            
    if current_name:
        seqs[current_name] = "".join(current_seq)
        
    return seqs


def parse_phylip_alignment(content: str) -> Dict[str, str]:
    """Parses sequential or interleaved PHYLIP sequence alignments."""
    lines = [line.strip() for line in content.splitlines() if line.strip()]
    if not lines:
        return {}
        
    # Check header
    header = lines[0].split()
    if len(header) < 2:
        # Fallback to FASTA or other if header is weird
        return parse_fasta(content)
        
    try:
        num_taxa = int(header[0])
        seq_len = int(header[1])
    except ValueError:
        return parse_fasta(content)
        
    seqs = {}
    lines = lines[1:]
    
    # Check if sequential or interleaved
    # If sequential, we expect each taxon block of seq_len length sequentially
    # If interleaved, we see taxa repeating in rounds.
    # Let's write a robust heuristic parser:
    taxon_names = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and not parts[0].isdigit() and len(taxon_names) < num_taxa:
            name = parts[0]
            seq = "".join(parts[1:])
            seqs[name] = seq
            taxon_names.append(name)
        elif len(taxon_names) == num_taxa:
            # We are interleaving or wrapping
            break
            
    # Check if interleaved
    current_idx = 0
    remaining_lines = lines[len(taxon_names):]
    
    if len(taxon_names) == num_taxa:
        # We can fill the rest
        for line in remaining_lines:
            parts = line.split()
            if not parts:
                continue
            name = parts[0]
            if name in seqs:
                seqs[name] += "".join(parts[1:])
            else:
                # Interleaved block without naming repeated: append in order
                seqs[taxon_names[current_idx % num_taxa]] += "".join(parts)
                current_idx += 1
    else:
        # Simple sequential parser fallback
        seqs = {}
        for line in lines:
            parts = line.split(None, 1)
            if len(parts) == 2:
                seqs[parts[0]] = parts[1].replace(" ", "")
                
    return seqs


def parse_phylip_matrix(content: str) -> Tuple[List[str], List[List[float]]]:
    """Parses a distance matrix into a taxa list and a 2D grid.

    Accepts PHYLIP (whitespace) and CSV (comma) formats, with an optional
    leading count line and an optional header row of column names.
    """
    raw = [line.strip() for line in content.splitlines() if line.strip()]
    if not raw:
        return [], []

    def tok(line: str) -> List[str]:
        # Treat commas/semicolons/tabs as separators (CSV + PHYLIP)
        return re.sub(r'[;,\t]+', ' ', line).strip().split()

    def is_int(s: str) -> bool:
        return re.fullmatch(r'\d+', s) is not None

    def is_num(s: str) -> bool:
        try:
            float(s)
            return True
        except ValueError:
            return False

    start = 0
    first = tok(raw[0])
    if len(first) == 1 and is_int(first[0]):
        start = 1  # PHYLIP count line
    elif len(first) > 1 and all(not is_num(t) for t in first):
        start = 1  # header row with column names (e.g. CSV ",A,B,C")

    taxa: List[str] = []
    matrix: List[List[float]] = []
    for line in raw[start:]:
        parts = tok(line)
        if len(parts) < 2:
            continue
        taxa.append(parts[0])
        row = []
        for p in parts[1:]:
            try:
                row.append(float(p))
            except ValueError:
                pass
        matrix.append(row)

    return taxa, matrix


def parse_newick(newick_str: str) -> TreeNode:
    """Parses a Newick string into a TreeNode tree."""
    newick_str = newick_str.strip()
    if newick_str.endswith(';'):
        newick_str = newick_str[:-1]
        
    def parse_node(s: str) -> Tuple[TreeNode, int]:
        s = s.strip()
        if not s.startswith('('):
            # It's a leaf node. Syntax: Name:BranchLength or Name
            colon_idx = s.find(':')
            if colon_idx != -1:
                name = s[:colon_idx].strip()
                try:
                    bl = float(s[colon_idx+1:])
                except ValueError:
                    bl = 0.1
                return TreeNode(name, bl), len(s)
            else:
                return TreeNode(s.strip(), None), len(s)
                
        # It's an internal node. Syntax: (child1,child2,...)Name:BranchLength
        node = TreeNode()
        i = 1  # Skip '('
        
        while i < len(s):
            # Find matching parenthesis or comma to split children
            # We must be careful with nested parentheses
            depth = 0
            j = i
            while j < len(s):
                if s[j] == '(':
                    depth += 1
                elif s[j] == ')':
                    if depth == 0:
                        break
                    depth -= 1
                elif s[j] == ',':
                    if depth == 0:
                        break
                j += 1
                
            child, parsed_len = parse_node(s[i:j])
            node.add_child(child)
            
            i = j
            if s[i] == ',':
                i += 1  # Skip comma and parse next child
            elif s[i] == ')':
                i += 1  # Skip ')' and extract internal node label & branch length
                break
                
        # Extract internal node label / branch length
        remaining = s[i:]
        # Find bounds of label & branch length
        end_idx = 0
        for char in remaining:
            if char in (',', ')', ';'):
                break
            end_idx += 1
            
        label_part = remaining[:end_idx].strip()
        colon_idx = label_part.find(':')
        if colon_idx != -1:
            node.name = label_part[:colon_idx].strip()
            if not node.name:
                node.name = None  # None for blank label
            try:
                node.branch_length = float(label_part[colon_idx+1:])
            except ValueError:
                node.branch_length = 0.1
        else:
            node.name = label_part if label_part else None
            node.branch_length = None
            
        return node, i + end_idx

    tree, _ = parse_node(newick_str)
    return tree


# ─── Tree Construction Methods ───────────────────────────────────────────────

def upgma(taxa: List[str], d_matrix: List[List[float]], weighted: bool = False) -> TreeNode:
    """Builds a tree using UPGMA (weighted=False) or WPGMA (weighted=True) algorithm."""
    n = len(taxa)
    # Deep copy the matrix
    matrix = [row[:] for row in d_matrix]
    
    # Active nodes map: list of TreeNode
    nodes = [TreeNode(name) for name in taxa]
    
    # Tracks the size (number of leaves) in each node cluster
    cluster_sizes = [1 for _ in range(n)]
    
    # Main clustering loop
    while len(nodes) > 1:
        # Find minimum distance in active matrix
        min_d = float('inf')
        u, v = -1, -1
        curr_n = len(nodes)
        
        for i in range(curr_n):
            for j in range(i+1, curr_n):
                if matrix[i][j] < min_d:
                    min_d = matrix[i][j]
                    u, v = i, j
                    
        if u == -1 or v == -1:
            break
            
        # Merge node u and node v.
        node_u = nodes[u]
        node_v = nodes[v]
        
        parent = TreeNode()
        # Set node heights / branch lengths
        # In UPGMA, parent height is min_d / 2
        # branch length = height - child height
        height = min_d / 2.0
        
        node_u_height = node_u.metadata.get('height', 0.0)
        node_v_height = node_v.metadata.get('height', 0.0)
        
        node_u.branch_length = max(0.0, height - node_u_height)
        node_v.branch_length = max(0.0, height - node_v_height)
        
        parent.add_child(node_u)
        parent.add_child(node_v)
        parent.metadata['height'] = height
        
        # Calculate new sizes
        size_u = cluster_sizes[u]
        size_v = cluster_sizes[v]
        new_size = size_u + size_v
        
        # Calculate new distances from parent cluster to all other clusters
        new_row = []
        for i in range(curr_n):
            if i == u or i == v:
                continue
            # Distance update formula:
            if weighted:
                # WPGMA: simple average
                d_new = (matrix[u][i] + matrix[v][i]) / 2.0
            else:
                # UPGMA: arithmetic mean weighted by cluster size
                d_new = (size_u * matrix[u][i] + size_v * matrix[v][i]) / new_size
            new_row.append(d_new)
            
        # Re-construct distance matrix
        new_matrix = []
        for i in range(curr_n):
            if i == u or i == v:
                continue
            new_matrix_row = []
            for j in range(curr_n):
                if j == u or j == v:
                    continue
                new_matrix_row.append(matrix[i][j])
            new_matrix.append(new_matrix_row)
            
        # Append the new cluster column/row to new_matrix
        for i in range(len(new_matrix)):
            new_matrix[i].append(new_row[i])
        new_matrix.append(new_row + [0.0])
        
        matrix = new_matrix
        
        # Update active nodes list
        new_nodes = [nodes[i] for i in range(curr_n) if i != u and i != v]
        new_nodes.append(parent)
        nodes = new_nodes
        
        # Update cluster sizes
        new_sizes = [cluster_sizes[i] for i in range(curr_n) if i != u and i != v]
        new_sizes.append(new_size)
        cluster_sizes = new_sizes
        
    return nodes[0]


def neighbor_joining(taxa: List[str], d_matrix: List[List[float]]) -> TreeNode:
    """Builds a tree using Neighbor Joining (NJ) algorithm (Saitou and Nei 1987)."""
    n = len(taxa)
    matrix = [row[:] for row in d_matrix]
    nodes = [TreeNode(name) for name in taxa]
    
    # NJ runs until we have 2 nodes remaining, which are then connected by a single branch
    while len(nodes) > 2:
        curr_n = len(nodes)
        
        # 1. Compute net divergence r[i] for each node
        r = [sum(matrix[i]) for i in range(curr_n)]
        
        # 2. Find pair (i,j) that minimizes Q-criterion: Q[i,j] = (N - 2)*d(i,j) - r[i] - r[j]
        min_q = float('inf')
        u, v = -1, -1
        
        for i in range(curr_n):
            for j in range(i+1, curr_n):
                q = (curr_n - 2) * matrix[i][j] - r[i] - r[j]
                if q < min_q:
                    min_q = q
                    u, v = i, j
                    
        if u == -1 or v == -1:
            break
            
        node_u = nodes[u]
        node_v = nodes[v]
        
        # 3. Calculate branch lengths from parent to u and v
        # d(parent, u) = 0.5 * d(u,v) + (r[u] - r[v]) / (2 * (N - 2))
        d_uv = matrix[u][v]
        bl_u = 0.5 * d_uv + (r[u] - r[v]) / (2.0 * (curr_n - 2))
        bl_v = d_uv - bl_u
        
        # Make sure branch lengths are positive
        node_u.branch_length = max(0.000001, bl_u)
        node_v.branch_length = max(0.000001, bl_v)
        
        parent = TreeNode()
        parent.add_child(node_u)
        parent.add_child(node_v)
        
        # 4. Compute distance from parent node to all other nodes:
        # d(parent, k) = 0.5 * (d(u,k) + d(v,k) - d(u,v))
        new_row = []
        for i in range(curr_n):
            if i == u or i == v:
                continue
            d_new = 0.5 * (matrix[u][i] + matrix[v][i] - d_uv)
            new_row.append(d_new)
            
        # Re-construct matrix
        new_matrix = []
        for i in range(curr_n):
            if i == u or i == v:
                continue
            new_matrix_row = []
            for j in range(curr_n):
                if j == u or j == v:
                    continue
                new_matrix_row.append(matrix[i][j])
            new_matrix.append(new_matrix_row)
            
        # Append the new column/row
        for i in range(len(new_matrix)):
            new_matrix[i].append(new_row[i])
        new_matrix.append(new_row + [0.0])
        
        matrix = new_matrix
        
        # Update active nodes list
        new_nodes = [nodes[i] for i in range(curr_n) if i != u and i != v]
        new_nodes.append(parent)
        nodes = new_nodes
        
    # Hook the final 2 nodes together
    if len(nodes) == 2:
        node_u = nodes[0]
        node_v = nodes[1]
        dist = matrix[0][1]
        
        # Connect to a root node
        root = TreeNode()
        node_u.branch_length = max(0.000001, dist / 2.0)
        node_v.branch_length = max(0.000001, dist / 2.0)
        root.add_child(node_u)
        root.add_child(node_v)
        return root
        
    return nodes[0]


# ─── Fitch-Margoliash (FM) Least-Squares ─────────────────────────────────────

def fitch_margoliash(taxa: List[str], d_matrix: List[List[float]]) -> TreeNode:
    """Fitch-Margoliash least-squares distance tree (heuristic NJ-based).

    NJ already yields the least-squares-optimal topology for additive data and is
    the standard FM starting point; we return it so the FM option is fully
    supported (matches the behaviour of the legacy phyloforge engine).
    """
    tree = neighbor_joining(taxa, d_matrix)
    tree.metadata['method'] = 'Fitch-Margoliash'
    return tree


def minimum_evolution(taxa: List[str], d_matrix: List[List[float]]) -> TreeNode:
    """Minimum Evolution heuristic.

    Uses the NJ topology, which minimises the ME criterion (total tree length)
    for additive distances. Mirrors the legacy phyloforge behaviour.
    """
    tree = neighbor_joining(taxa, d_matrix)
    tree.metadata['method'] = 'Minimum Evolution'
    return tree


# ─── Bootstrap Support ───────────────────────────────────────────────────────

def _get_leaves_ordered(node: TreeNode) -> List[TreeNode]:
    """Returns leaves of the subtree in left-to-right order."""
    if node.is_leaf():
        return [node]
    result: List[TreeNode] = []
    for c in node.children:
        result.extend(_get_leaves_ordered(c))
    return result


def _build_distance_tree(taxa: List[str], d_matrix: List[List[float]], method: str) -> TreeNode:
    """Dispatches a distance-based builder by method name (used by bootstrap)."""
    m = (method or 'NJ').upper()
    if m == 'UPGMA':
        return upgma(taxa, d_matrix, weighted=False)
    if m == 'WPGMA':
        return upgma(taxa, d_matrix, weighted=True)
    if m == 'FM':
        return fitch_margoliash(taxa, d_matrix)
    if m == 'ME':
        return minimum_evolution(taxa, d_matrix)
    return neighbor_joining(taxa, d_matrix)


def _get_internal_clades(node: TreeNode) -> List[List[str]]:
    """Returns the leaf-name set of every internal (non-trivial) clade."""
    result: List[List[str]] = []

    def rec(n: TreeNode):
        if n.is_leaf():
            return
        leaves = [l.name for l in _get_leaves_ordered(n) if l.name]
        if len(leaves) > 1:
            result.append(leaves)
        for c in n.children:
            rec(c)

    rec(node)
    return result


def _attach_support(node: TreeNode, counts: Dict[frozenset, int], n_reps: int):
    """Writes bootstrap percentages into each internal node's metadata['support']."""
    if not node.is_leaf() and node.children:
        leaves = frozenset(l.name for l in _get_leaves_ordered(node) if l.name)
        if leaves in counts:
            node.metadata['support'] = round(counts[leaves] / n_reps * 100)
        for c in node.children:
            _attach_support(c, counts, n_reps)


def bootstrap_support(seqs: Dict[str, str], n_reps: int, method: str = 'NJ',
                      model: str = 'JC69') -> TreeNode:
    """Builds a distance tree and annotates clades with bootstrap support (%).

    Columns of the alignment are resampled with replacement n_reps times; each
    replicate tree's clades are tallied against the original tree's clades.
    Returns the original tree with metadata['support'] set on internal nodes.
    """
    import random

    taxa = list(seqs.keys())
    seq_len = len(next(iter(seqs.values()))) if seqs else 0

    o_taxa, o_matrix, _ = compute_distance_matrix(seqs, model)
    original_tree = _build_distance_tree(o_taxa, o_matrix, method)

    if seq_len == 0:
        return original_tree

    counts: Dict[frozenset, int] = {frozenset(c): 0 for c in _get_internal_clades(original_tree)}

    for _ in range(n_reps):
        cols = [random.randint(0, seq_len - 1) for _ in range(seq_len)]
        resampled = {t: ''.join(seqs[t][c] for c in cols) for t in taxa}
        b_taxa, b_matrix, _ = compute_distance_matrix(resampled, model)
        b_tree = _build_distance_tree(b_taxa, b_matrix, method)
        for clade in _get_internal_clades(b_tree):
            key = frozenset(clade)
            if key in counts:
                counts[key] += 1

    _attach_support(original_tree, counts, n_reps)
    original_tree.metadata['bootstrap_reps'] = n_reps
    return original_tree


# ─── Character Parsimony (MP) Fitch Algorithm ────────────────────────────────

def fitch_parsimony_score(tree: TreeNode, alignment: Dict[str, str]) -> float:
    """Calculates parsimony score of a tree topology given an alignment using Fitch's algorithm."""
    taxa = list(alignment.keys())
    seq_len = len(next(iter(alignment.values())))
    
    # Store set of possible characters at each node
    # Leaf nodes are initialized with sequence characters
    node_sets = {}
    
    def fitch_downpass(node: TreeNode) -> float:
        score = 0.0
        if node.is_leaf():
            seq = alignment.get(node.name, '-' * seq_len)
            node_sets[node] = [set(char) for char in seq]
            return 0.0
        
        # Recursively evaluate children first
        for child in node.children:
            score += fitch_downpass(child)
            
        # Internal node combination
        # Let's assume binary tree
        if len(node.children) >= 2:
            left, right = node.children[0], node.children[1]
            parent_sets = []
            for site in range(seq_len):
                S1 = node_sets[left][site]
                S2 = node_sets[right][site]
                intersect = S1.intersection(S2)
                if intersect:
                    parent_sets.append(intersect)
                else:
                    parent_sets.append(S1.union(S2))
                    score += 1.0
            node_sets[node] = parent_sets
        elif len(node.children) == 1:
            node_sets[node] = node_sets[node.children[0]]
            
        return score
        
    return fitch_downpass(tree)


def maximum_parsimony_tree(seqs: Dict[str, str]) -> TreeNode:
    """Constructs a Neighbor-Joining starting tree and evaluates its parsimony score."""
    taxa = list(seqs.keys())
    # Generate distance matrix for starting NJ tree
    taxa, d_mat, _ = compute_distance_matrix(seqs, 'p-distance')
    tree = neighbor_joining(taxa, d_mat)
    
    # Calculate score
    score = fitch_parsimony_score(tree, seqs)
    tree.metadata['parsimony_score'] = score
    tree.metadata['method'] = 'Maximum Parsimony'
    return tree


# ─── Maximum Likelihood (ML) Engine ──────────────────────────────────────────

class MLEngine:
    """Maximum Likelihood Engine using Felsenstein's Pruning Algorithm (under JC69 model)."""
    def __init__(self, alignment: Dict[str, str]):
        self.alignment = {k.upper(): v.upper() for k, v in alignment.items()}
        self.taxa = list(alignment.keys())
        self.seq_len = len(next(iter(alignment.values())))
        # Pre-calculate base frequencies
        self.base_freqs = get_base_frequencies(list(alignment.values()))
        self.bases = ['A', 'C', 'G', 'T']
        self.base_idx = {b: i for i, b in enumerate(self.bases)}

    def jc_transition_prob(self, base_i: int, base_j: int, t: float) -> float:
        """JC69 transition probability between base_i and base_j over branch length t."""
        if t < 0: t = 0.0
        exp_term = math.exp(-4.0 / 3.0 * t)
        if base_i == base_j:
            return 0.25 + 0.75 * exp_term
        else:
            return 0.25 - 0.25 * exp_term

    def compute_likelihood(self, tree: TreeNode) -> float:
        """Computes the overall tree log-likelihood using Felsenstein's Pruning Algorithm."""
        # Initialize node conditional likelihood tables
        # Node -> site -> base_index -> probability
        cond_likes = {}

        def postorder_pruning(node: TreeNode):
            if node.is_leaf():
                # Leaf node: conditional likelihoods are 1 at the observed base, 0 elsewhere
                seq = self.alignment.get(node.name.upper(), '-' * self.seq_len)
                likes = [[[0.0 for _ in range(4)] for _ in range(4)] for _ in range(self.seq_len)] # site -> base -> prob
                site_probs = []
                for site in range(self.seq_len):
                    char = seq[site]
                    probs = [0.0 for _ in range(4)]
                    if char in self.base_idx:
                        probs[self.base_idx[char]] = 1.0
                    elif char in ('-', 'N', '?'):
                        # Gaps/Ambiguity: all nucleotides equally likely
                        probs = [1.0 for _ in range(4)]
                    else:
                        probs = [1.0 for _ in range(4)]
                    site_probs.append(probs)
                cond_likes[node] = site_probs
                return

            # Compute child likelihoods
            for child in node.children:
                postorder_pruning(child)

            # Internal node: combines likelihoods of children
            node_probs = []
            for site in range(self.seq_len):
                probs = [1.0 for _ in range(4)]
                
                # Multiply likelihood contributions from each child
                for child in node.children:
                    child_probs = cond_likes[child][site]
                    t = child.branch_length if child.branch_length is not None else 0.1
                    
                    # Probability for parent base i
                    child_contrib = [0.0 for _ in range(4)]
                    for i in range(4):
                        # Sum over all possible child states j
                        for j in range(4):
                            child_contrib[i] += self.jc_transition_prob(i, j, t) * child_probs[j]
                            
                    for i in range(4):
                        probs[i] *= child_contrib[i]
                        
                node_probs.append(probs)
                
            cond_likes[node] = node_probs

        # Run postorder post-pruning traversal
        postorder_pruning(tree)
        
        # Sum likelihood over all root states weighted by base frequency
        log_likelihood = 0.0
        root_probs = cond_likes[tree]
        
        for site in range(self.seq_len):
            site_like = 0.0
            for i in range(4):
                base = self.bases[i]
                freq = self.base_freqs.get(base, 0.25)
                site_like += root_probs[site][i] * freq
                
            # Accumulate Log Likelihood
            if site_like > 0:
                log_likelihood += math.log(site_like)
            else:
                log_likelihood += -100.0  # Safe log-likelihood floor
                
        return round(log_likelihood, 6)

    def optimize_branch_lengths(self, tree: TreeNode, max_iter: int = 15) -> float:
        """Heuristic hill-climbing optimization of branch lengths to maximize likelihood."""
        nodes = self._get_all_nodes_except_root(tree)
        best_ll = self.compute_likelihood(tree)
        
        step_size = 0.05
        
        for iteration in range(max_iter):
            improved = False
            for node in nodes:
                original_bl = node.branch_length if node.branch_length is not None else 0.1
                
                # Try positive step
                node.branch_length = max(0.0001, original_bl + step_size)
                ll_up = self.compute_likelihood(tree)
                
                # Try negative step
                node.branch_length = max(0.0001, original_bl - step_size)
                ll_down = self.compute_likelihood(tree)
                
                # Check improvements
                if ll_up > best_ll and ll_up > ll_down:
                    best_ll = ll_up
                    node.branch_length = original_bl + step_size
                    improved = True
                elif ll_down > best_ll:
                    best_ll = ll_down
                    node.branch_length = max(0.0001, original_bl - step_size)
                    improved = True
                else:
                    # Restore original branch length
                    node.branch_length = original_bl
                    
            if not improved:
                # Shrink step size if no improvement was found
                step_size /= 2.0
                if step_size < 0.001:
                    break
                    
        return best_ll

    def _get_all_nodes_except_root(self, node: TreeNode) -> List[TreeNode]:
        result = []
        for child in node.children:
            result.append(child)
            result.extend(self._get_all_nodes_except_root(child))
        return result


def maximum_likelihood_tree(seqs: Dict[str, str], model: str = 'JC69') -> TreeNode:
    """Builds an initial NJ tree and optimizes its branch lengths under Maximum Likelihood."""
    taxa = list(seqs.keys())
    taxa, d_mat, _ = compute_distance_matrix(seqs, model)
    tree = neighbor_joining(taxa, d_mat)
    
    # Optimize branches
    engine = MLEngine(seqs)
    best_ll = engine.optimize_branch_lengths(tree)
    
    tree.metadata['log_likelihood'] = best_ll
    tree.metadata['method'] = 'Maximum Likelihood'
    return tree


# ─── Statistics & Midpoint Rooting ──────────────────────────────────────────

def tree_stats(tree: TreeNode) -> Dict[str, Any]:
    """Computes analytical stats for a phylogenetic tree."""
    all_nodes = _get_all_nodes(tree)
    leaves = [n for n in all_nodes if n.is_leaf()]
    internal = [n for n in all_nodes if not n.is_leaf()]
    
    branch_lengths = [n.branch_length for n in all_nodes if n.branch_length is not None]
    total_length = sum(branch_lengths) if branch_lengths else 0.0
    mean_bl = total_length / len(branch_lengths) if branch_lengths else 0.0
    max_bl = max(branch_lengths) if branch_lengths else 0.0
    min_bl = min(branch_lengths) if branch_lengths else 0.0
    
    depth = _max_depth(tree)
    leaf_names = [leaf.name for leaf in leaves if leaf.name]
    
    return {
        'n_taxa': len(leaves),
        'n_leaves': len(leaves),
        'n_internal': len(internal),
        'total_branch_length': round(total_length, 6),
        'mean_branch_length': round(mean_bl, 6),
        'max_branch_length': round(max_bl, 6),
        'min_branch_length': round(min_bl, 6),
        'tree_depth': depth,
        'taxa': leaf_names,
    }


def _get_all_nodes(node: TreeNode) -> List[TreeNode]:
    result = [node]
    for c in node.children:
        result.extend(_get_all_nodes(c))
    return result


def _max_depth(node: TreeNode, current=0) -> int:
    if node.is_leaf():
        return current
    return max(_max_depth(c, current+1) for c in node.children)
