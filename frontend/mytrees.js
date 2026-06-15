/**
 * MyTrees JS - Client-Side Phylogenetic Analysis Engine
 * Implements: UPGMA, WPGMA, NJ, Fitch Parsimony, Maximum Likelihood (Felsenstein's Pruning)
 * Distance Models: p-distance, JC69, K2P, F84, LogDet
 */

(function(window) {
    'use strict';

    const MyTreesEngine = {};

    // ─── TreeNode Data Structure ───────────────────────────────────────────────

    class TreeNode {
        constructor(name = null, branchLength = null) {
            this.name = name;
            this.branchLength = branchLength;
            this.children = [];
            this.parent = null;
            this.metadata = {};
        }

        isLeaf() {
            return this.children.length === 0;
        }

        addChild(child) {
            child.parent = this;
            this.children.push(child);
        }

        removeChild(child) {
            const index = this.children.indexOf(child);
            if (index !== -1) {
                this.children.splice(index, 1);
                child.parent = null;
            }
        }

        toNewick() {
            if (this.isLeaf()) {
                const namePart = this.name ? this.name.replace(/[\(\):;,]/g, '_') : "";
                if (this.branchLength !== null) {
                    return `${namePart}:${this.branchLength.toFixed(6)}`;
                }
                return namePart;
            } else {
                const childrenStr = this.children.map(c => c.toNewick()).join(",");
                if (this.branchLength !== null) {
                    return `(${childrenStr}):${this.branchLength.toFixed(6)}`;
                }
                return `(${childrenStr})`;
            }
        }

        clone() {
            const newNode = new TreeNode(this.name, this.branchLength);
            newNode.metadata = JSON.parse(JSON.stringify(this.metadata));
            for (const child of this.children) {
                newNode.addChild(child.clone());
            }
            return newNode;
        }
    }

    MyTreesEngine.TreeNode = TreeNode;

    /**
     * Rebuilds a TreeNode from the structured dict produced by the Python
     * backend's TreeNode.to_dict(). Unlike re-parsing Newick, this preserves
     * internal-node metadata such as bootstrap support. The backend field
     * `length` maps to this engine's `branchLength`.
     */
    MyTreesEngine.dictToTreeNode = function dictToTreeNode(d) {
        if (!d) return null;
        const len = (d.length !== undefined && d.length !== null) ? d.length : null;
        const node = new TreeNode(d.name ? d.name : null, len);
        node.metadata = d.metadata ? JSON.parse(JSON.stringify(d.metadata)) : {};
        if (Array.isArray(d.children)) {
            for (const childDict of d.children) {
                const child = dictToTreeNode(childDict);
                if (child) node.addChild(child);
            }
        }
        return node;
    };

    // ─── Distance Models ─────────────────────────────────────────────────────────

    function getBaseFrequencies(seqs) {
        const counts = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 };
        let total = 0;
        for (const seq of seqs) {
            for (let i = 0; i < seq.length; i++) {
                const char = seq[i].toUpperCase();
                if (char in counts) {
                    counts[char]++;
                    total++;
                }
            }
        }
        if (total === 0) {
            return { 'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25 };
        }
        return {
            'A': counts['A'] / total,
            'C': counts['C'] / total,
            'G': counts['G'] / total,
            'T': counts['T'] / total
        };
    }

    function pDistance(seq1, seq2) {
        let n = 0;
        let diffs = 0;
        const len = Math.min(seq1.length, seq2.length);
        for (let i = 0; i < len; i++) {
            const a = seq1[i].toUpperCase();
            const b = seq2[i].toUpperCase();
            if (a !== '-' && a !== '?' && b !== '-' && b !== '?') {
                n++;
                if (a !== b) {
                    diffs++;
                }
            }
        }
        return n === 0 ? 0.0 : diffs / n;
    }

    function jc69(seq1, seq2) {
        const p = pDistance(seq1, seq2);
        if (p >= 0.75) return 5.0;
        if (p === 0) return 0.0;
        return -0.75 * Math.log(1.0 - (4.0 / 3.0) * p);
    }

    function k2p(seq1, seq2) {
        let n = 0;
        let transitions = 0;
        let transversions = 0;
        const purines = new Set(['A', 'G']);
        const pyrimidines = new Set(['C', 'T']);
        const len = Math.min(seq1.length, seq2.length);

        for (let i = 0; i < len; i++) {
            const a = seq1[i].toUpperCase();
            const b = seq2[i].toUpperCase();
            if (a !== '-' && a !== '?' && b !== '-' && b !== '?') {
                n++;
                if (a !== b) {
                    const aIsPurine = purines.has(a);
                    const bIsPurine = purines.has(b);
                    const aIsPyrimidine = pyrimidines.has(a);
                    const bIsPyrimidine = pyrimidines.has(b);
                    
                    if ((aIsPurine && bIsPurine) || (aIsPyrimidine && bIsPyrimidine)) {
                        transitions++;
                    } else {
                        transversions++;
                    }
                }
            }
        }

        if (n === 0) return 0.0;
        const P = transitions / n;
        const Q = transversions / n;

        const val1 = 1.0 - 2.0 * P - Q;
        const val2 = 1.0 - 2.0 * Q;

        if (val1 <= 0 || val2 <= 0) return 5.0;
        return -0.5 * Math.log(val1) - 0.25 * Math.log(val2);
    }

    function f84(seq1, seq2, baseFreqs = null) {
        if (!baseFreqs) {
            baseFreqs = getBaseFrequencies([seq1, seq2]);
        }
        const fA = baseFreqs['A'] || 0.25;
        const fC = baseFreqs['C'] || 0.25;
        const fG = baseFreqs['G'] || 0.25;
        const fT = baseFreqs['T'] || 0.25;

        const fR = fA + fG; // Purines
        const fY = fC + fT; // Pyrimidines

        let n = 0;
        let transitions = 0;
        let transversions = 0;
        const purines = new Set(['A', 'G']);
        const pyrimidines = new Set(['C', 'T']);
        const len = Math.min(seq1.length, seq2.length);

        for (let i = 0; i < len; i++) {
            const a = seq1[i].toUpperCase();
            const b = seq2[i].toUpperCase();
            if (a !== '-' && a !== '?' && b !== '-' && b !== '?') {
                n++;
                if (a !== b) {
                    const aIsPurine = purines.has(a);
                    const bIsPurine = purines.has(b);
                    const aIsPyrimidine = pyrimidines.has(a);
                    const bIsPyrimidine = pyrimidines.has(b);

                    if ((aIsPurine && bIsPurine) || (aIsPyrimidine && bIsPyrimidine)) {
                        transitions++;
                    } else {
                        transversions++;
                    }
                }
            }
        }

        if (n === 0) return 0.0;
        const P = transitions / n;
        const Q = transversions / n;

        const a_coeff = (fA * fG / fR) + (fC * fT / fY);
        const b_coeff = fR * fY;

        const val1 = 1.0 - Q / (2.0 * b_coeff);
        if (val1 <= 0) return 5.0;

        const val2 = 1.0 - P / (2.0 * a_coeff) - ((fA * fG * fY / fR) + (fC * fT * fR / fY)) * Q / (2.0 * a_coeff * b_coeff);
        if (val2 <= 0) return 5.0;

        return -4.0 * a_coeff * Math.log(val2) - 4.0 * (b_coeff - a_coeff) * Math.log(val1);
    }

    function logdet(seq1, seq2) {
        let n = 0;
        const bases = ['A', 'C', 'G', 'T'];
        const baseIdx = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
        const F = Array(4).fill(0).map(() => Array(4).fill(0));
        const len = Math.min(seq1.length, seq2.length);

        for (let i = 0; i < len; i++) {
            const a = seq1[i].toUpperCase();
            const b = seq2[i].toUpperCase();
            if (a in baseIdx && b in baseIdx) {
                F[baseIdx[a]][baseIdx[b]]++;
                n++;
            }
        }

        if (n === 0) return 0.0;
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                F[i][j] /= n;
            }
        }

        // Calculate 4x4 determinant
        const det_F = (
            F[0][0]*(F[1][1]*(F[2][2]*F[3][3] - F[2][3]*F[3][2]) - F[1][2]*(F[2][1]*F[3][3] - F[2][3]*F[3][1]) + F[1][3]*(F[2][1]*F[3][2] - F[2][2]*F[3][1])) -
            F[0][1]*(F[1][0]*(F[2][2]*F[3][3] - F[2][3]*F[3][2]) - F[1][2]*(F[2][0]*F[3][3] - F[2][3]*F[3][0]) + F[1][3]*(F[2][0]*F[3][2] - F[2][2]*F[3][0])) +
            F[0][2]*(F[1][0]*(F[2][1]*F[3][3] - F[2][3]*F[3][1]) - F[1][1]*(F[2][0]*F[3][3] - F[2][3]*F[3][0]) + F[1][3]*(F[2][0]*F[3][1] - F[2][1]*F[3][0])) -
            F[0][3]*(F[1][0]*(F[2][1]*F[3][2] - F[2][2]*F[3][1]) - F[1][1]*(F[2][0]*F[3][2] - F[2][2]*F[3][0]) + F[1][2]*(F[2][0]*F[3][1] - F[2][1]*F[3][0]))
        );

        if (det_F <= 0) return 5.0;

        const f1 = Array(4).fill(0);
        const f2 = Array(4).fill(0);
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                f1[i] += F[i][j];
                f2[i] += F[j][i];
            }
        }

        let prod_f1 = 1.0;
        let prod_f2 = 1.0;
        for (let i = 0; i < 4; i++) {
            if (f1[i] > 0) prod_f1 *= f1[i];
            if (f2[i] > 0) prod_f2 *= f2[i];
        }

        return -0.25 * (Math.log(det_F) - 0.5 * (Math.log(prod_f1) + Math.log(prod_f2)));
    }

    MyTreesEngine.computeDistanceMatrix = function(seqs, model = 'JC69') {
        const taxa = Object.keys(seqs);
        const n = taxa.length;
        const matrix = Array(n).fill(0).map(() => Array(n).fill(0.0));
        const baseFreqs = getBaseFrequencies(Object.values(seqs));

        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                const s1 = seqs[taxa[i]];
                const s2 = seqs[taxa[j]];
                let d = 0.0;
                if (model === 'p-distance') {
                    d = pDistance(s1, s2);
                } else if (model === 'JC69') {
                    d = jc69(s1, s2);
                } else if (model === 'K2P') {
                    d = k2p(s1, s2);
                } else if (model === 'F84') {
                    d = f84(s1, s2, baseFreqs);
                } else if (model === 'LogDet') {
                    d = logdet(s1, s2);
                } else {
                    d = jc69(s1, s2);
                }
                
                matrix[i][j] = parseFloat(d.toFixed(6));
                matrix[j][i] = parseFloat(d.toFixed(6));
            }
        }

        const stats = analyzeMatrixStats(taxa, matrix);
        return { taxa, matrix, stats };
    };

    function analyzeMatrixStats(taxa, matrix) {
        const n = taxa.length;
        const vals = [];
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                vals.push(matrix[i][j]);
            }
        }
        if (vals.length === 0) return {};
        const sum = vals.reduce((a, b) => a + b, 0);
        const mean = sum / vals.length;
        const variance = vals.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / vals.length;
        const stdDev = Math.sqrt(variance);

        let violations = 0;
        for (let i = 0; i < n; i++) {
            for (let j = i + 1; j < n; j++) {
                for (let k = j + 1; k < n; k++) {
                    if (matrix[i][j] > matrix[i][k] + matrix[k][j] + 1e-5) violations++;
                    if (matrix[i][k] > matrix[i][j] + matrix[j][k] + 1e-5) violations++;
                    if (matrix[j][k] > matrix[j][i] + matrix[i][k] + 1e-5) violations++;
                }
            }
        }

        return {
            n_taxa: n,
            min_distance: parseFloat(Math.min(...vals).toFixed(6)),
            max_distance: parseFloat(Math.max(...vals).toFixed(6)),
            mean_distance: parseFloat(mean.toFixed(6)),
            std_distance: parseFloat(stdDev.toFixed(6)),
            is_symmetric: true,
            triangle_violations: violations,
            n_pairs: vals.length
        };
    }

    MyTreesEngine.analyzeMatrixStats = analyzeMatrixStats;

    // ─── Parsers ─────────────────────────────────────────────────────────────────

    MyTreesEngine.parseFASTA = function(text) {
        const seqs = {};
        let currentTaxon = null;
        const lines = text.split('\n');
        
        for (let line of lines) {
            line = line.trim();
            if (!line) continue;
            if (line.startsWith('>')) {
                currentTaxon = line.substring(1).trim();
                seqs[currentTaxon] = "";
            } else if (currentTaxon) {
                seqs[currentTaxon] += line.replace(/\s/g, "");
            }
        }
        return seqs;
    };

    MyTreesEngine.parsePHYLIPAlignment = function(text) {
        const lines = text.split('\n').map(l => l.trim()).filter(l => l.length > 0);
        if (lines.length === 0) return {};
        const header = lines[0].split(/\s+/);
        if (header.length < 2) return MyTreesEngine.parseFASTA(text);

        const numTaxa = parseInt(header[0]);
        if (isNaN(numTaxa)) return MyTreesEngine.parseFASTA(text);

        const seqs = {};
        const parsedTaxa = [];
        let i = 1;
        
        // Read initial round of sequences
        while (i < lines.length && parsedTaxa.length < numTaxa) {
            const parts = lines[i].split(/\s+/);
            if (parts.length >= 2) {
                const name = parts[0];
                const seq = parts.slice(1).join("");
                seqs[name] = seq;
                parsedTaxa.push(name);
            }
            i++;
        }

        // Interleaved contents continuation
        let idx = 0;
        while (i < lines.length) {
            const line = lines[i];
            const parts = line.split(/\s+/);
            if (parts.length > 0) {
                const targetTaxon = parsedTaxa[idx % numTaxa];
                // Check if interleaving rows carry names
                if (parts[0] === targetTaxon && parts.length >= 2) {
                    seqs[targetTaxon] += parts.slice(1).join("");
                } else {
                    seqs[targetTaxon] += parts.join("");
                }
                idx++;
            }
            i++;
        }

        return seqs;
    };

    MyTreesEngine.parsePHYLIPMatrix = function(text) {
        // Aceita PHYLIP (espaços) e CSV (vírgulas), com linha de contagem e/ou
        // cabeçalho de colunas opcionais.
        const lines = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length > 0);
        if (lines.length === 0) return { taxa: [], matrix: [] };

        const tok = l => l.replace(/[;,\t]+/g, ' ').trim().split(/\s+/);
        const isInt = s => /^\d+$/.test(s);
        const isNum = s => s !== '' && !isNaN(Number(s));

        let start = 0;
        const first = tok(lines[0]);
        if (first.length === 1 && isInt(first[0])) {
            start = 1; // linha de contagem PHYLIP
        } else if (first.length > 1 && first.every(s => !isNum(s))) {
            start = 1; // cabeçalho com nomes de colunas (ex.: CSV ",A,B,C")
        }

        const taxa = [];
        const matrix = [];
        for (let i = start; i < lines.length; i++) {
            const parts = tok(lines[i]);
            if (parts.length < 2) continue;
            taxa.push(parts[0]);
            matrix.push(parts.slice(1).map(Number).filter(x => !isNaN(x)));
        }

        return { taxa, matrix };
    };

    MyTreesEngine.parseNewick = function(str) {
        str = str.trim();
        if (str.endsWith(';')) str = str.slice(0, -1);

        function parseNode(s) {
            s = s.trim();
            if (!s.startsWith('(')) {
                const colonIdx = s.indexOf(':');
                if (colonIdx !== -1) {
                    const name = s.substring(0, colonIdx).trim();
                    const bl = parseFloat(s.substring(colonIdx + 1));
                    return [new TreeNode(name, isNaN(bl) ? 0.1 : bl), s.length];
                }
                return [new TreeNode(s, null), s.length];
            }

            const node = new TreeNode();
            let i = 1;
            
            while (i < s.length) {
                let depth = 0;
                let j = i;
                while (j < s.length) {
                    if (s[j] === '(') depth++;
                    else if (s[j] === ')') {
                        if (depth === 0) break;
                        depth--;
                    } else if (s[j] === ',') {
                        if (depth === 0) break;
                    }
                    j++;
                }

                const [child, parsedLen] = parseNode(s.substring(i, j));
                node.addChild(child);

                i = j;
                if (s[i] === ',') i++;
                else if (s[i] === ')') {
                    i++;
                    break;
                }
            }

            // Extract labels/branch length for internal node
            const remaining = s.substring(i);
            let endIdx = 0;
            for (let k = 0; k < remaining.length; k++) {
                if (remaining[k] === ',' || remaining[k] === ')' || remaining[k] === ';') break;
                endIdx++;
            }

            const labelPart = remaining.substring(0, endIdx).trim();
            const colonIdx = labelPart.indexOf(':');
            if (colonIdx !== -1) {
                const name = labelPart.substring(0, colonIdx).trim();
                node.name = name === "" ? null : name;
                const bl = parseFloat(labelPart.substring(colonIdx + 1));
                node.branchLength = isNaN(bl) ? 0.1 : bl;
            } else {
                node.name = labelPart === "" ? null : labelPart;
                node.branchLength = null;
            }

            return [node, i + endIdx];
        }

        const [tree] = parseNode(str);
        return tree;
    };

    // ─── Algorithms ─────────────────────────────────────────────────────────────

    MyTreesEngine.upgma = function(taxa, dMatrix, weighted = false) {
        const n = taxa.length;
        let matrix = dMatrix.map(row => row.slice());
        let nodes = taxa.map(name => new TreeNode(name));
        let clusterSizes = Array(n).fill(1);

        while (nodes.length > 1) {
            let minD = Infinity;
            let u = -1, v = -1;
            const currN = nodes.length;

            for (let i = 0; i < currN; i++) {
                for (let j = i + 1; j < currN; j++) {
                    if (matrix[i][j] < minD) {
                        minD = matrix[i][j];
                        u = i;
                        v = j;
                    }
                }
            }

            if (u === -1 || v === -1) break;

            const nodeU = nodes[u];
            const nodeV = nodes[v];

            const parent = new TreeNode();
            const height = minD / 2.0;

            const uHeight = nodeU.metadata.height || 0.0;
            const vHeight = nodeV.metadata.height || 0.0;

            nodeU.branchLength = Math.max(0.0, height - uHeight);
            nodeV.branchLength = Math.max(0.0, height - vHeight);

            parent.addChild(nodeU);
            parent.addChild(nodeV);
            parent.metadata.height = height;

            const sizeU = clusterSizes[u];
            const sizeV = clusterSizes[v];
            const newSize = sizeU + sizeV;

            // Recalculate cluster distances
            const newRow = [];
            for (let i = 0; i < currN; i++) {
                if (i === u || i === v) continue;
                let dNew = 0.0;
                if (weighted) {
                    dNew = (matrix[u][i] + matrix[v][i]) / 2.0;
                } else {
                    dNew = (sizeU * matrix[u][i] + sizeV * matrix[v][i]) / newSize;
                }
                newRow.push(dNew);
            }

            // Slice out columns & rows
            const newMatrix = [];
            for (let i = 0; i < currN; i++) {
                if (i === u || i === v) continue;
                const newMatrixRow = [];
                for (let j = 0; j < currN; j++) {
                    if (j === u || j === v) continue;
                    newMatrixRow.push(matrix[i][j]);
                }
                newMatrix.push(newMatrixRow);
            }

            for (let i = 0; i < newMatrix.length; i++) {
                newMatrix[i].push(newRow[i]);
            }
            newMatrix.push([...newRow, 0.0]);

            matrix = newMatrix;

            nodes = nodes.filter((_, idx) => idx !== u && idx !== v);
            nodes.push(parent);

            clusterSizes = clusterSizes.filter((_, idx) => idx !== u && idx !== v);
            clusterSizes.push(newSize);
        }

        return nodes[0];
    };

    MyTreesEngine.neighborJoining = function(taxa, dMatrix) {
        const n = taxa.length;
        let matrix = dMatrix.map(row => row.slice());
        let nodes = taxa.map(name => new TreeNode(name));

        while (nodes.length > 2) {
            const currN = nodes.length;
            const r = matrix.map(row => row.reduce((a, b) => a + b, 0.0));

            let minQ = Infinity;
            let u = -1, v = -1;

            for (let i = 0; i < currN; i++) {
                for (let j = i + 1; j < currN; j++) {
                    const q = (currN - 2) * matrix[i][j] - r[i] - r[j];
                    if (q < minQ) {
                        minQ = q;
                        u = i;
                        v = j;
                    }
                }
            }

            if (u === -1 || v === -1) break;

            const nodeU = nodes[u];
            const nodeV = nodes[v];
            const d_uv = matrix[u][v];

            const blU = 0.5 * d_uv + (r[u] - r[v]) / (2.0 * (currN - 2));
            const blV = d_uv - blU;

            nodeU.branchLength = Math.max(0.000001, blU);
            nodeV.branchLength = Math.max(0.000001, blV);

            const parent = new TreeNode();
            parent.addChild(nodeU);
            parent.addChild(nodeV);

            // Compute distance from parent to others
            const newRow = [];
            for (let i = 0; i < currN; i++) {
                if (i === u || i === v) continue;
                const dNew = 0.5 * (matrix[u][i] + matrix[v][i] - d_uv);
                newRow.push(dNew);
            }

            // Build new matrix
            const newMatrix = [];
            for (let i = 0; i < currN; i++) {
                if (i === u || i === v) continue;
                const newMatrixRow = [];
                for (let j = 0; j < currN; j++) {
                    if (j === u || j === v) continue;
                    newMatrixRow.push(matrix[i][j]);
                }
                newMatrix.push(newMatrixRow);
            }

            for (let i = 0; i < newMatrix.length; i++) {
                newMatrix[i].push(newRow[i]);
            }
            newMatrix.push([...newRow, 0.0]);

            matrix = newMatrix;

            nodes = nodes.filter((_, idx) => idx !== u && idx !== v);
            nodes.push(parent);
        }

        if (nodes.length === 2) {
            const nodeU = nodes[0];
            const nodeV = nodes[1];
            const dist = matrix[0][1];

            const root = new TreeNode();
            nodeU.branchLength = Math.max(0.000001, dist / 2.0);
            nodeV.branchLength = Math.max(0.000001, dist / 2.0);
            root.addChild(nodeU);
            root.addChild(nodeV);
            return root;
        }

        return nodes[0];
    };

    // ─── Fitch Parsimony Score ──────────────────────────────────────────────────

    MyTreesEngine.fitchParsimony = function(tree, alignment) {
        const taxa = Object.keys(alignment);
        const seqLen = alignment[taxa[0]].length;
        const nodeSets = new Map();

        function downpass(node) {
            let score = 0;
            if (node.isLeaf()) {
                const seq = alignment[node.name] || "-".repeat(seqLen);
                const leafSets = [];
                for (let i = 0; i < seqLen; i++) {
                    leafSets.push(new Set([seq[i].toUpperCase()]));
                }
                nodeSets.set(node, leafSets);
                return 0;
            }

            for (const child of node.children) {
                score += downpass(child);
            }

            if (node.children.length >= 2) {
                const left = node.children[0];
                const right = node.children[1];
                const parentSets = [];

                for (let site = 0; site < seqLen; site++) {
                    const S1 = nodeSets.get(left)[site];
                    const S2 = nodeSets.get(right)[site];
                    
                    const intersect = new Set([...S1].filter(x => S2.has(x)));
                    if (intersect.size > 0) {
                        parentSets.push(intersect);
                    } else {
                        const union = new Set([...S1, ...S2]);
                        parentSets.push(union);
                        score++;
                    }
                }
                nodeSets.set(node, parentSets);
            } else if (node.children.length === 1) {
                nodeSets.set(node, nodeSets.get(node.children[0]));
            }

            return score;
        }

        return downpass(tree);
    };

    // ─── Maximum Likelihood (Pruning JC69) ────────────────────────────────────────

    class MLEngine {
        constructor(alignment) {
            this.alignment = {};
            for (const k of Object.keys(alignment)) {
                this.alignment[k.toUpperCase()] = alignment[k].toUpperCase();
            }
            this.taxa = Object.keys(alignment);
            this.seqLen = alignment[this.taxa[0]].length;
            this.baseFreqs = getBaseFrequencies(Object.values(alignment));
            this.bases = ['A', 'C', 'G', 'T'];
            this.baseIdx = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 };
        }

        jcTransitionProb(base_i, base_j, t) {
            if (t < 0) t = 0.0;
            const expTerm = Math.exp(-4.0 / 3.0 * t);
            if (base_i === base_j) {
                return 0.25 + 0.75 * expTerm;
            } else {
                return 0.25 - 0.25 * expTerm;
            }
        }

        computeLikelihood(tree) {
            const condLikes = new Map();
            const self = this;

            function postorder(node) {
                if (node.isLeaf()) {
                    const seq = self.alignment[node.name.toUpperCase()] || "-".repeat(self.seqLen);
                    const siteProbs = [];
                    for (let site = 0; site < self.seqLen; site++) {
                        const char = seq[site];
                        const probs = [0.0, 0.0, 0.0, 0.0];
                        if (char in self.baseIdx) {
                            probs[self.baseIdx[char]] = 1.0;
                        } else {
                            // Gaps/ambiguity
                            for (let idx = 0; idx < 4; idx++) probs[idx] = 1.0;
                        }
                        siteProbs.push(probs);
                    }
                    condLikes.set(node, siteProbs);
                    return;
                }

                for (const child of node.children) {
                    postorder(child);
                }

                const nodeProbs = [];
                for (let site = 0; site < self.seqLen; site++) {
                    const probs = [1.0, 1.0, 1.0, 1.0];
                    for (const child of node.children) {
                        const childProbs = condLikes.get(child)[site];
                        const t = child.branchLength !== null ? child.branchLength : 0.1;
                        const childContrib = [0.0, 0.0, 0.0, 0.0];

                        for (let i = 0; i < 4; i++) {
                            for (let j = 0; j < 4; j++) {
                                childContrib[i] += self.jcTransitionProb(i, j, t) * childProbs[j];
                            }
                        }

                        for (let i = 0; i < 4; i++) {
                            probs[i] *= childContrib[i];
                        }
                    }
                    nodeProbs.push(probs);
                }
                condLikes.set(node, nodeProbs);
            }

            postorder(tree);

            let logLikelihood = 0.0;
            const rootProbs = condLikes.get(tree);

            for (let site = 0; site < self.seqLen; site++) {
                let siteLike = 0.0;
                for (let i = 0; i < 4; i++) {
                    const base = self.bases[i];
                    const freq = self.baseFreqs[base] || 0.25;
                    siteLike += rootProbs[site][i] * freq;
                }
                if (siteLike > 0) {
                    logLikelihood += Math.log(siteLike);
                } else {
                    logLikelihood += -100.0;
                }
            }

            return parseFloat(logLikelihood.toFixed(6));
        }

        optimizeBranchLengths(tree, maxIter = 15) {
            const nodes = [];
            function gather(node) {
                for (const c of node.children) {
                    nodes.push(c);
                    gather(c);
                }
            }
            gather(tree);

            let bestLL = this.computeLikelihood(tree);
            let stepSize = 0.05;

            for (let iteration = 0; iteration < maxIter; iteration++) {
                let improved = false;
                for (const node of nodes) {
                    const originalBl = node.branchLength !== null ? node.branchLength : 0.1;

                    // Try positive step
                    node.branchLength = Math.max(0.0001, originalBl + stepSize);
                    const llUp = this.computeLikelihood(tree);

                    // Try negative step
                    node.branchLength = Math.max(0.0001, originalBl - stepSize);
                    const llDown = this.computeLikelihood(tree);

                    if (llUp > bestLL && llUp > llDown) {
                        bestLL = llUp;
                        node.branchLength = originalBl + stepSize;
                        improved = true;
                    } else if (llDown > bestLL) {
                        bestLL = llDown;
                        node.branchLength = Math.max(0.0001, originalBl - stepSize);
                        improved = true;
                    } else {
                        node.branchLength = originalBl;
                    }
                }

                if (!improved) {
                    stepSize /= 2.0;
                    if (stepSize < 0.001) break;
                }
            }

            return bestLL;
        }
    }

    MyTreesEngine.maximumLikelihood = function(seqs, model = 'JC69') {
        const taxa = Object.keys(seqs);
        const { matrix } = MyTreesEngine.computeDistanceMatrix(seqs, model);
        const tree = MyTreesEngine.neighborJoining(taxa, matrix);
        
        const engine = new MLEngine(seqs);
        const bestLL = engine.optimizeBranchLengths(tree);
        
        tree.metadata.log_likelihood = bestLL;
        tree.metadata.method = 'Maximum Likelihood';
        return tree;
    };

    MyTreesEngine.maximumParsimony = function(seqs) {
        const taxa = Object.keys(seqs);
        const { matrix } = MyTreesEngine.computeDistanceMatrix(seqs, 'p-distance');
        const tree = MyTreesEngine.neighborJoining(taxa, matrix);
        const score = MyTreesEngine.fitchParsimony(tree, seqs);
        tree.metadata.parsimony_score = score;
        tree.metadata.method = 'Maximum Parsimony';
        return tree;
    };

    // ─── Statistics helper ──────────────────────────────────────────────────────

    MyTreesEngine.treeStats = function(tree) {
        const nodes = [];
        function gather(node) {
            nodes.push(node);
            for (const child of node.children) {
                gather(child);
            }
        }
        gather(tree);

        const leaves = nodes.filter(n => n.isLeaf());
        const internal = nodes.filter(n => !n.isLeaf());
        
        const branchLengths = nodes.map(n => n.branchLength).filter(x => x !== null);
        const totalLength = branchLengths.reduce((a, b) => a + b, 0.0);
        const meanBl = branchLengths.length > 0 ? totalLength / branchLengths.length : 0.0;
        
        function getMaxDepth(node, current = 0) {
            if (node.isLeaf()) return current;
            return Math.max(...node.children.map(c => getMaxDepth(c, current + 1)));
        }

        const depth = getMaxDepth(tree);
        const leafNames = leaves.map(l => l.name).filter(x => x);

        return {
            n_taxa: leaves.length,
            n_leaves: leaves.length,
            n_internal: internal.length,
            total_branch_length: parseFloat(totalLength.toFixed(6)),
            mean_branch_length: parseFloat(meanBl.toFixed(6)),
            tree_depth: depth,
            taxa: leafNames
        };
    };

    // Register globally
    window.MyTreesEngine = MyTreesEngine;

})(window);
