/**
 * MyTrees Visualization Engine
 * Renders phylogenetic trees in SVG (Phylogram, Cladogram, Circular Radial)
 * Supports interactive node styling, swapping children, collapsing clades, and rerooting.
 */

(function(window) {
    'use strict';

    const MyTreesViz = {
        tree: null,
        layout: 'rectangular-phylogram', // rectangular-phylogram, rectangular-cladogram, circular-phylogram, circular-cladogram
        widthScale: 1.0,
        heightScale: 1.0,
        labelSizeScale: 1.0,
        selectedNode: null,
        nodeStyles: new Map(), // node -> { stroke, strokeWidth, strokeDasharray }
        collapsedNodes: new Set(), // Set of TreeNode objects that are collapsed
        labelsConfig: {
            showTaxa: true,
            taxaColor: '#e2f1f8',
            taxaSize: 12,
            taxaFont: 'Inter',
            taxaWeight: '400',
            
            showBootstrap: false,
            bootstrapColor: '#4eb1b1',
            bootstrapSize: 10,
            bootstrapFont: 'IBM Plex Mono',
            bootstrapWeight: '500',
            
            showBranchLength: false,
            branchLengthColor: '#5b7083',
            branchLengthSize: 9,
            branchLengthFont: 'IBM Plex Mono',
            branchLengthWeight: '400'
        },
        onNodeSelectedCallback: null
    };

    // Helper: Assign unique IDs to all nodes in the tree
    let nodeIndex = 0;
    function assignNodeIds(node) {
        if (!node.id) {
            node.id = 'node_' + (nodeIndex++);
        }
        for (const child of node.children) {
            assignNodeIds(child);
        }
    }

    MyTreesViz.initTree = function(tree) {
        nodeIndex = 0;
        assignNodeIds(tree);
        this.tree = tree;
        this.selectedNode = null;
        this.collapsedNodes.clear();
        this.nodeStyles.clear();
    };

    // Helper: Calculate tree depths and heights for layouts
    function computeDepthsAndHeights(node, depth = 0, dist = 0.0) {
        node.depth = depth;
        node.distanceFromRoot = dist;
        
        let maxDepth = depth;
        let maxDist = dist;
        let leavesCount = node.children.length === 0 ? 1 : 0;
        
        for (const child of node.children) {
            const bl = child.branchLength !== null ? child.branchLength : 0.1;
            const res = computeDepthsAndHeights(child, depth + 1, dist + bl);
            maxDepth = Math.max(maxDepth, res.maxDepth);
            maxDist = Math.max(maxDist, res.maxDist);
            leavesCount += res.leavesCount;
        }
        
        node.maxDepthBelow = maxDepth;
        node.maxDistBelow = maxDist;
        node.leavesCount = leavesCount;
        
        return { maxDepth, maxDist, leavesCount };
    }

    // Helper: Find Lowest Common Ancestor (LCA) of two nodes
    function findLCA(nodeA, nodeB) {
        const pathA = [];
        let curr = nodeA;
        while (curr) {
            pathA.push(curr);
            curr = curr.parent;
        }
        
        curr = nodeB;
        while (curr) {
            if (pathA.includes(curr)) {
                return curr;
            }
            curr = curr.parent;
        }
        return null;
    }

    // Helper: Get all leaves of a subtree
    function getLeaves(node, list = []) {
        if (node.children.length === 0 || MyTreesViz.collapsedNodes.has(node)) {
            list.push(node);
        } else {
            for (const child of node.children) {
                getLeaves(child, list);
            }
        }
        return list;
    }

    // Helper: Calculate absolute distance from a node to one of its leaves
    function getPathDistance(start, end) {
        const lca = findLCA(start, end);
        if (!lca) return Infinity;
        
        let dist = 0;
        let curr = start;
        while (curr !== lca) {
            dist += curr.branchLength || 0.0;
            curr = curr.parent;
        }
        curr = end;
        while (curr !== lca) {
            dist += curr.branchLength || 0.0;
            curr = curr.parent;
        }
        return dist;
    }

    // ─── Swap Children (Rotate Node) ───────────────────────────────────────────
    MyTreesViz.rotateNode = function(node) {
        if (node && node.children.length > 1) {
            node.children.reverse();
            return true;
        }
        return false;
    };

    // ─── Collapse/Expand Clade ─────────────────────────────────────────────────
    MyTreesViz.toggleCollapseClade = function(node) {
        if (!node || node.children.length === 0) return false;
        if (this.collapsedNodes.has(node)) {
            this.collapsedNodes.delete(node);
        } else {
            this.collapsedNodes.add(node);
            if (this.selectedNode === node || isDescendantOf(this.selectedNode, node)) {
                this.selectedNode = node;
            }
        }
        return true;
    };

    function isDescendantOf(childNode, parentNode) {
        if (!childNode) return false;
        let curr = childNode.parent;
        while (curr) {
            if (curr === parentNode) return true;
            curr = curr.parent;
        }
        return false;
    }

    // ─── Recursive Path Reversal for Rerooting ──────────────────────────────────
    function reversePath(node, newParent, newBranchLength) {
        const parent = node.parent;
        
        // Remove connection to old parent
        if (parent) {
            parent.removeChild(node);
        }
        
        // If there's an old parent, we reverse the link recursively
        if (parent) {
            reversePath(parent, node, node.branchLength);
        }
        
        // Add new parent as child
        if (newParent) {
            node.addChild(newParent);
            newParent.branchLength = newBranchLength;
        }
        
        // If we are at the old root and it now has only 1 child, we collapse it
        if (!parent && node.children.length === 1) {
            const onlyChild = node.children[0];
            node.removeChild(onlyChild);
            if (node.parent) {
                const nodeParent = node.parent;
                nodeParent.removeChild(node);
                nodeParent.addChild(onlyChild);
                onlyChild.branchLength = (onlyChild.branchLength || 0.0) + (node.branchLength || 0.0);
            } else {
                // It was the root, so the onlyChild is now the root
                onlyChild.parent = null;
                onlyChild.branchLength = null;
                MyTreesViz.tree = onlyChild;
            }
        }
    }

    // ─── Reroot Tree at Selected Node/Branch ────────────────────────────────────
    MyTreesViz.rerootAtNode = function(targetNode) {
        if (!targetNode || !targetNode.parent) return false; // Already root or invalid
        
        const parent = targetNode.parent;
        const totalBranchLength = targetNode.branchLength !== null ? targetNode.branchLength : 0.1;
        
        // We will insert a new root node directly on the branch between targetNode and parent
        const newRoot = new window.MyTreesEngine.TreeNode(null, null);
        newRoot.id = 'node_new_root_' + Math.floor(Math.random() * 10000);
        
        // Remove targetNode from parent
        parent.removeChild(targetNode);
        
        // Set new branch lengths (split original in half)
        targetNode.branchLength = totalBranchLength / 2.0;
        newRoot.addChild(targetNode);
        
        // Reverse parent path starting from parent
        reversePath(parent, newRoot, totalBranchLength / 2.0);
        
        // Update tree pointer
        this.tree = newRoot;
        assignNodeIds(this.tree);
        this.selectedNode = null;
        return true;
    };

    // ─── Midpoint Rooting ──────────────────────────────────────────────────────
    MyTreesViz.midpointRoot = function() {
        if (!this.tree) return false;
        
        // 1. Gather all leaves
        const leaves = [];
        function gatherLeaves(node) {
            if (node.children.length === 0) {
                leaves.push(node);
            } else {
                for (const c of node.children) gatherLeaves(c);
            }
        }
        gatherLeaves(this.tree);
        
        if (leaves.length < 2) return false;
        
        // 2. Find the pair of leaves with the maximum path distance (tree diameter)
        let maxDist = -1;
        let leafA = null;
        let leafB = null;
        
        for (let i = 0; i < leaves.length; i++) {
            for (let j = i + 1; j < leaves.length; j++) {
                const dist = getPathDistance(leaves[i], leaves[j]);
                if (dist > maxDist) {
                    maxDist = dist;
                    leafA = leaves[i];
                    leafB = leaves[j];
                }
            }
        }
        
        if (!leafA || !leafB || maxDist <= 0) return false;
        
        // 3. Find the exact midpoint of this path
        const targetMidpoint = maxDist / 2.0;
        
        // Trace path from leafA to leafB
        const lca = findLCA(leafA, leafB);
        const pathA = []; // from leafA up to LCA
        let curr = leafA;
        while (curr !== lca) {
            pathA.push(curr);
            curr = curr.parent;
        }
        
        const pathB = []; // from leafB up to LCA
        curr = leafB;
        while (curr !== lca) {
            pathB.push(curr);
            curr = curr.parent;
        }
        
        // Combine path: leafA -> ... -> LCA -> ... -> leafB
        const path = [...pathA, lca, ...pathB.reverse()];
        
        // Traverse path to find which branch contains the midpoint
        let accumulatedDist = 0.0;
        let midpointNode = null;
        let splitDist = 0.0;
        
        // The path array lists nodes. Path distance from leafA:
        // path[0] = leafA, path[1] = leafA's parent, etc.
        // The distance contribution is path[i].branchLength as we move from leafA towards leafB.
        // Wait, for nodes moving down from LCA to leafB, their branch lengths represent the branch leading to them.
        
        // Let's build a precise list of segments and their lengths
        const segments = []; // elements: { child, parent, length }
        
        // Upward path from leafA to LCA
        for (let i = 0; i < pathA.length; i++) {
            const child = pathA[i];
            const p = child.parent;
            segments.push({ child, parent: p, length: child.branchLength || 0.0 });
        }
        // Downward path from LCA to leafB
        const downPath = [lca, ...pathB];
        for (let i = 0; i < downPath.length - 1; i++) {
            const p = downPath[i];
            const child = downPath[i+1];
            segments.push({ child, parent: p, length: child.branchLength || 0.0 });
        }
        
        // Trace along segments from leafA (start of segments) to leafB
        let dist = 0.0;
        let chosenSeg = null;
        for (const seg of segments) {
            if (dist + seg.length >= targetMidpoint) {
                chosenSeg = seg;
                splitDist = targetMidpoint - dist; // distance from seg.child to midpoint
                break;
            }
            dist += seg.length;
        }
        
        if (!chosenSeg) {
            // Fallback: use the branch leading to LCA
            if (pathA.length > 0) {
                chosenSeg = { child: pathA[pathA.length-1], parent: lca, length: pathA[pathA.length-1].branchLength || 0.1 };
                splitDist = chosenSeg.length / 2;
            } else {
                return false;
            }
        }
        
        // 4. Reroot at the midpoint on the chosen branch!
        // We will insert a new root node on the branch between chosenSeg.child and chosenSeg.parent
        const childNode = chosenSeg.child;
        const parentNode = chosenSeg.parent;
        
        const newRoot = new window.MyTreesEngine.TreeNode(null, null);
        newRoot.id = 'node_midpoint_root_' + Math.floor(Math.random() * 10000);
        
        // Remove childNode from parentNode
        parentNode.removeChild(childNode);
        
        // Split distances
        childNode.branchLength = splitDist;
        newRoot.addChild(childNode);
        
        // Reverse parent path starting from parentNode
        const remainingLength = chosenSeg.length - splitDist;
        reversePath(parentNode, newRoot, remainingLength);
        
        this.tree = newRoot;
        assignNodeIds(this.tree);
        this.selectedNode = null;
        return true;
    };


    // ─── Tree Layout Coordinates Calculation ────────────────────────────────────
    MyTreesViz.calculateCoordinates = function(width, height) {
        if (!this.tree) return;

        // Reset details
        const info = computeDepthsAndHeights(this.tree);
        const maxDepth = info.maxDepth;
        const maxDist = info.maxDist || 0.1;
        const totalLeaves = info.leavesCount;

        const isCircular = this.layout.startsWith('circular');
        const isCladogram = this.layout.endsWith('cladogram');

        if (!isCircular) {
            // ─── Rectangular Layouts ───
            const paddingLeft = 30;
            const paddingRight = 150 * this.labelSizeScale; // Spacing for leaf labels
            const paddingTop = 40;
            const paddingBottom = 40;

            const drawWidth = (width - paddingLeft - paddingRight) * this.widthScale;
            const drawHeight = (height - paddingTop - paddingBottom) * this.heightScale;

            let leafIndex = 0;

            function layoutPass(node, currX) {
                // X position: depth-based (Cladogram) or branch-length-based (Phylogram)
                if (isCladogram) {
                    node.x = paddingLeft + (node.depth / maxDepth) * drawWidth;
                } else {
                    node.x = paddingLeft + (node.distanceFromRoot / maxDist) * drawWidth;
                }

                if (node.children.length === 0 || MyTreesViz.collapsedNodes.has(node)) {
                    // Leaf node (or collapsed internal node treated as leaf)
                    node.y = paddingTop + (leafIndex / Math.max(1, totalLeaves - 1)) * drawHeight;
                    leafIndex++;
                } else {
                    // Internal node
                    let sumY = 0;
                    for (const child of node.children) {
                        layoutPass(child, node.x);
                        sumY += child.y;
                    }
                    node.y = sumY / node.children.length;
                }
            }

            layoutPass(this.tree, paddingLeft);

        } else {
            // ─── Circular/Radial Layouts ───
            const centerX = width / 2;
            const centerY = height / 2;
            const maxRadius = (Math.min(width, height) / 2 - 100 * this.labelSizeScale) * Math.min(this.widthScale, this.heightScale);

            let leafIndex = 0;

            function circularLayoutPass(node) {
                // Radius (R)
                if (isCladogram) {
                    node.radius = (node.depth / maxDepth) * maxRadius;
                } else {
                    node.radius = (node.distanceFromRoot / maxDist) * maxRadius;
                }

                if (node.children.length === 0 || MyTreesViz.collapsedNodes.has(node)) {
                    // Leaf angle (distributed around the circle)
                    node.angle = (leafIndex / totalLeaves) * 2 * Math.PI - Math.PI / 2; // offset so first node is at top
                    leafIndex++;
                } else {
                    // Internal node angle
                    for (const child of node.children) {
                        circularLayoutPass(child);
                    }
                    // Compute average angle of children (handling angle wrapping)
                    let sinSum = 0;
                    let cosSum = 0;
                    for (const child of node.children) {
                        sinSum += Math.sin(child.angle);
                        cosSum += Math.cos(child.angle);
                    }
                    node.angle = Math.atan2(sinSum, cosSum);
                }

                // Convert Polar to Cartesian
                node.x = centerX + node.radius * Math.cos(node.angle);
                node.y = centerY + node.radius * Math.sin(node.angle);
            }

            circularLayoutPass(this.tree);
        }
    };

    // ─── Render SVG ────────────────────────────────────────────────────────────
    MyTreesViz.render = function(svgElement) {
        if (!this.tree) return;

        const parent = svgElement.parentElement;
        const baseWidth = parent ? parent.clientWidth : 800;
        const baseHeight = parent ? parent.clientHeight : 550;

        const scaledWidth = Math.max(baseWidth * this.widthScale, baseWidth);
        const scaledHeight = Math.max(baseHeight * this.heightScale, baseHeight);
        
        svgElement.style.width = `${scaledWidth}px`;
        svgElement.style.height = `${scaledHeight}px`;
        
        svgElement.setAttribute('viewBox', `0 0 ${scaledWidth} ${scaledHeight}`);
        svgElement.innerHTML = ''; // Clear SVG

        // Recompute positions. We pass the BASE width and height to calculateCoordinates 
        // because it multiplies by widthScale/heightScale internally for the tree paths.
        this.calculateCoordinates(baseWidth, baseHeight);

        // SVG Groups for logical order (background elements first, active paths, overlays, labels last)
        const gBranches = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const gCollapsedTriangles = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const gOverlays = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const gLabels = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const gNodes = document.createElementNS('http://www.w3.org/2000/svg', 'g');

        svgElement.appendChild(gBranches);
        svgElement.appendChild(gCollapsedTriangles);
        svgElement.appendChild(gOverlays);
        svgElement.appendChild(gLabels);
        svgElement.appendChild(gNodes);

        const self = this;

        // Recursive SVG elements rendering
        function draw(node) {
            const style = self.nodeStyles.get(node.id) || {
                stroke: '#e2f1f8',
                strokeWidth: 2,
                strokeDasharray: 'none'
            };

            // Selected node highlight coloring
            const isSelected = (self.selectedNode === node);
            const activeColor = isSelected ? '#f16e53' : style.stroke;
            const activeWidth = isSelected ? style.strokeWidth + 2.5 : style.strokeWidth;
            const activeDash = style.strokeDasharray;

            // Draw branch to children
            if (node.children.length > 0 && !self.collapsedNodes.has(node)) {
                for (const child of node.children) {
                    const childStyle = self.nodeStyles.get(child.id) || {
                        stroke: '#e2f1f8',
                        strokeWidth: 2,
                        strokeDasharray: 'none'
                    };
                    const isChildSelected = (self.selectedNode === child);
                    const branchColor = isChildSelected ? '#f16e53' : childStyle.stroke;
                    const branchWidth = isChildSelected ? childStyle.strokeWidth + 2.5 : childStyle.strokeWidth;
                    const branchDash = childStyle.strokeDasharray;

                    if (!self.layout.startsWith('circular')) {
                        // ─── Rectangular Branch Paths ───
                        // 1. Vertical split connector at parent X
                        const vLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        vLine.setAttribute('x1', node.x);
                        vLine.setAttribute('y1', node.y);
                        vLine.setAttribute('x2', node.x);
                        vLine.setAttribute('y2', child.y);
                        vLine.setAttribute('stroke', activeColor);
                        vLine.setAttribute('stroke-width', activeWidth);
                        if (activeDash !== 'none') vLine.setAttribute('stroke-dasharray', activeDash);
                        gBranches.appendChild(vLine);

                        // 2. Horizontal branch to child
                        const hLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        hLine.setAttribute('x1', node.x);
                        hLine.setAttribute('y1', child.y);
                        hLine.setAttribute('x2', child.x);
                        hLine.setAttribute('y2', child.y);
                        hLine.setAttribute('stroke', branchColor);
                        hLine.setAttribute('stroke-width', branchWidth);
                        if (branchDash !== 'none') hLine.setAttribute('stroke-dasharray', branchDash);
                        gBranches.appendChild(hLine);

                        // 3. Invisible clickable overlay spanning parent -> child (horizontal branch)
                        const hOverlay = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        hOverlay.setAttribute('x1', node.x);
                        hOverlay.setAttribute('y1', child.y);
                        hOverlay.setAttribute('x2', child.x);
                        hOverlay.setAttribute('y2', child.y);
                        hOverlay.setAttribute('stroke', 'transparent');
                        hOverlay.setAttribute('stroke-width', '10');
                        hOverlay.setAttribute('cursor', 'pointer');
                        hOverlay.addEventListener('click', (e) => {
                            e.stopPropagation();
                            self.selectNode(child);
                        });
                        gOverlays.appendChild(hOverlay);

                        // Vertical connector overlay (maps clicks to the parent split)
                        const vOverlay = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        vOverlay.setAttribute('x1', node.x);
                        vOverlay.setAttribute('y1', node.y);
                        vOverlay.setAttribute('x2', node.x);
                        vOverlay.setAttribute('y2', child.y);
                        vOverlay.setAttribute('stroke', 'transparent');
                        vOverlay.setAttribute('stroke-width', '10');
                        vOverlay.setAttribute('cursor', 'pointer');
                        vOverlay.addEventListener('click', (e) => {
                            e.stopPropagation();
                            self.selectNode(node);
                        });
                        gOverlays.appendChild(vOverlay);

                    } else {
                        // ─── Circular/Radial Branch Paths ───
                        const centerX = width / 2;
                        const centerY = height / 2;

                        // 1. Radial line from parent radius to child radius at child's angle
                        const rLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        const startX = centerX + node.radius * Math.cos(child.angle);
                        const startY = centerY + node.radius * Math.sin(child.angle);
                        rLine.setAttribute('x1', startX);
                        rLine.setAttribute('y1', startY);
                        rLine.setAttribute('x2', child.x);
                        rLine.setAttribute('y2', child.y);
                        rLine.setAttribute('stroke', branchColor);
                        rLine.setAttribute('stroke-width', branchWidth);
                        if (branchDash !== 'none') rLine.setAttribute('stroke-dasharray', branchDash);
                        gBranches.appendChild(rLine);

                        // Invisible clickable overlay for the radial branch
                        const rOverlay = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                        rOverlay.setAttribute('x1', startX);
                        rOverlay.setAttribute('y1', startY);
                        rOverlay.setAttribute('x2', child.x);
                        rOverlay.setAttribute('y2', child.y);
                        rOverlay.setAttribute('stroke', 'transparent');
                        rOverlay.setAttribute('stroke-width', '12');
                        rOverlay.setAttribute('cursor', 'pointer');
                        rOverlay.addEventListener('click', (e) => {
                            e.stopPropagation();
                            self.selectNode(child);
                        });
                        gOverlays.appendChild(rOverlay);
                    }

                    draw(child);
                }

                // If circular layout, we also draw the circular arc uniting all child angles at the parent's radius
                if (self.layout.startsWith('circular')) {
                    const centerX = width / 2;
                    const centerY = height / 2;
                    const angles = node.children.map(c => c.angle);
                    const minAngle = Math.min(...angles);
                    const maxAngle = Math.max(...angles);

                    if (maxAngle - minAngle > 0.001) {
                        const pathData = describeArc(centerX, centerY, node.radius, minAngle, maxAngle);
                        const arcPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                        arcPath.setAttribute('d', pathData);
                        arcPath.setAttribute('fill', 'none');
                        arcPath.setAttribute('stroke', activeColor);
                        arcPath.setAttribute('stroke-width', activeWidth);
                        if (activeDash !== 'none') arcPath.setAttribute('stroke-dasharray', activeDash);
                        gBranches.appendChild(arcPath);

                        // Invisible clickable arc overlay
                        const arcOverlay = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                        arcOverlay.setAttribute('d', pathData);
                        arcOverlay.setAttribute('fill', 'none');
                        arcOverlay.setAttribute('stroke', 'transparent');
                        arcOverlay.setAttribute('stroke-width', '12');
                        arcOverlay.setAttribute('cursor', 'pointer');
                        arcOverlay.addEventListener('click', (e) => {
                            e.stopPropagation();
                            self.selectNode(node);
                        });
                        gOverlays.appendChild(arcOverlay);
                    }
                }
            }

            // Draw Collapsed Clade Triangle
            if (self.collapsedNodes.has(node)) {
                // A solid triangle representing a collapsed clade
                const leaves = getLeaves(node);
                const info = computeDepthsAndHeights(node);
                const maxDepth = info.maxDepth;
                const maxDist = info.maxDist || 0.1;

                // Find depth boundaries in Cartesian or Polar
                let apexX = node.x;
                let apexY = node.y;
                let baseX, baseY1, baseY2;

                const style = self.nodeStyles.get(node.id) || { stroke: '#f16e53', strokeWidth: 2 };
                const triangleColor = isSelected ? '#f16e53' : style.stroke;
                const fillColor = isSelected ? 'rgba(78, 177, 177, 0.2)' : 'rgba(23, 74, 97, 0.08)';

                if (!self.layout.startsWith('circular')) {
                    // Rectangular collapsed clade triangle
                    const drawWidth = (width - 30 - (150 * self.labelSizeScale)) * self.widthScale;
                    const paddingLeft = 30;
                    
                    if (self.layout.endsWith('cladogram')) {
                        baseX = paddingLeft + (maxDepth / computeDepthsAndHeights(self.tree).maxDepth) * drawWidth;
                    } else {
                        baseX = paddingLeft + (node.maxDistBelow / computeDepthsAndHeights(self.tree).maxDist) * drawWidth;
                    }

                    // Height is proportion of drawing height spanned by collapsed items
                    const padding = 15;
                    baseY1 = node.y - padding - (leaves.length * 2.5);
                    baseY2 = node.y + padding + (leaves.length * 2.5);

                    const polygon = document.createElementNS('http://www.w3.org/2000/svg', 'polygon');
                    polygon.setAttribute('points', `${apexX},${apexY} ${baseX},${baseY1} ${baseX},${baseY2}`);
                    polygon.setAttribute('fill', fillColor);
                    polygon.setAttribute('stroke', triangleColor);
                    polygon.setAttribute('stroke-width', isSelected ? '3' : '1.5');
                    polygon.setAttribute('cursor', 'pointer');
                    polygon.addEventListener('click', (e) => {
                        e.stopPropagation();
                        self.selectNode(node);
                    });
                    gCollapsedTriangles.appendChild(polygon);

                    // Add collapsed clade label (number of sequences)
                    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                    label.setAttribute('x', baseX + 10);
                    label.setAttribute('y', node.y + 4);
                    label.setAttribute('fill', '#e2f1f8');
                    label.setAttribute('font-size', '10px');
                    label.setAttribute('font-family', 'Inter');
                    label.setAttribute('font-weight', '600');
                    label.textContent = `[${leaves.length} Táxons Colapsados]`;
                    gLabels.appendChild(label);

                } else {
                    // Circular collapsed clade triangle (drawn as pie slice or sector)
                    const centerX = width / 2;
                    const centerY = height / 2;
                    const maxRadius = (Math.min(width, height) / 2 - 100 * self.labelSizeScale) * Math.min(self.widthScale, self.heightScale);
                    const treeMaxDist = computeDepthsAndHeights(self.tree).maxDist || 0.1;
                    const treeMaxDepth = computeDepthsAndHeights(self.tree).maxDepth;

                    let baseRadius;
                    if (self.layout.endsWith('cladogram')) {
                        baseRadius = (maxDepth / treeMaxDepth) * maxRadius;
                    } else {
                        baseRadius = (node.maxDistBelow / treeMaxDist) * maxRadius;
                    }

                    const angleWidth = (leaves.length / computeDepthsAndHeights(self.tree).leavesCount) * Math.PI * 0.4 + 0.1;
                    const angle1 = node.angle - angleWidth / 2;
                    const angle2 = node.angle + angleWidth / 2;

                    const bx1 = centerX + baseRadius * Math.cos(angle1);
                    const by1 = centerY + baseRadius * Math.sin(angle1);
                    const bx2 = centerX + baseRadius * Math.cos(angle2);
                    const by2 = centerY + baseRadius * Math.sin(angle2);

                    const sectorPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');
                    const d = `M ${apexX} ${apexY} L ${bx1} ${by1} A ${baseRadius} ${baseRadius} 0 0 1 ${bx2} ${by2} Z`;
                    sectorPath.setAttribute('d', d);
                    sectorPath.setAttribute('fill', fillColor);
                    sectorPath.setAttribute('stroke', triangleColor);
                    sectorPath.setAttribute('stroke-width', isSelected ? '3' : '1.5');
                    sectorPath.setAttribute('cursor', 'pointer');
                    sectorPath.addEventListener('click', (e) => {
                        e.stopPropagation();
                        self.selectNode(node);
                    });
                    gCollapsedTriangles.appendChild(sectorPath);

                    // Label for circular collapsed
                    const labelX = centerX + (baseRadius + 15) * Math.cos(node.angle);
                    const labelY = centerY + (baseRadius + 15) * Math.sin(node.angle);
                    const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                    label.setAttribute('x', labelX);
                    label.setAttribute('y', labelY + 3);
                    label.setAttribute('fill', '#e2f1f8');
                    label.setAttribute('font-size', '9px');
                    label.setAttribute('font-family', 'Inter');
                    label.setAttribute('text-anchor', Math.cos(node.angle) > 0 ? 'start' : 'end');
                    label.textContent = `[${leaves.length} Táxons Colapsados]`;
                    gLabels.appendChild(label);
                }
            }

            // Draw Node Dot (clickable anchor point at splits)
            if (node.children.length > 0 && !self.collapsedNodes.has(node)) {
                const dot = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                dot.setAttribute('cx', node.x);
                dot.setAttribute('cy', node.y);
                dot.setAttribute('r', isSelected ? '5' : '3.5');
                dot.setAttribute('fill', isSelected ? '#f16e53' : '#ffffff');
                dot.setAttribute('stroke', isSelected ? '#f16e53' : '#e2f1f8');
                dot.setAttribute('stroke-width', '1.5');
                dot.setAttribute('cursor', 'pointer');
                dot.addEventListener('click', (e) => {
                    e.stopPropagation();
                    self.selectNode(node);
                });
                gNodes.appendChild(dot);
            }

            // ─── Render Text Labels (Leaf names, Bootstrap, Branch length) ───
            const isLeaf = (node.children.length === 0 || self.collapsedNodes.has(node));
            
            // 1. Leaf Labels (Taxa names)
            if (isLeaf && self.labelsConfig.showTaxa && node.name && !self.collapsedNodes.has(node)) {
                const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                
                if (!self.layout.startsWith('circular')) {
                    text.setAttribute('x', node.x + 10);
                    text.setAttribute('y', node.y + 4);
                    text.setAttribute('text-anchor', 'start');
                } else {
                    // Position labels outwards radially
                    const angleDeg = node.angle * (180 / Math.PI);
                    const isRight = Math.cos(node.angle) > 0;
                    
                    let rotateAngle = angleDeg;
                    if (!isRight) rotateAngle += 180; // keep text upright

                    text.setAttribute('x', node.x + (isRight ? 10 : -10));
                    text.setAttribute('y', node.y + 3);
                    text.setAttribute('text-anchor', isRight ? 'start' : 'end');
                    text.setAttribute('transform', `rotate(${rotateAngle}, ${node.x}, ${node.y})`);
                }

                text.setAttribute('fill', self.labelsConfig.taxaColor);
                text.setAttribute('font-size', `${self.labelsConfig.taxaSize * self.labelSizeScale}px`);
                text.setAttribute('font-family', self.labelsConfig.taxaFont);
                text.setAttribute('font-weight', self.labelsConfig.taxaWeight);
                text.textContent = node.name;
                gLabels.appendChild(text);
            }

            // 2. Node Labels (Bootstrap / Support values)
            if (self.labelsConfig.showBootstrap && !isLeaf && node.name) {
                const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                
                if (!self.layout.startsWith('circular')) {
                    text.setAttribute('x', node.x - 8);
                    text.setAttribute('y', node.y - 6);
                    text.setAttribute('text-anchor', 'end');
                } else {
                    text.setAttribute('x', node.x - 5);
                    text.setAttribute('y', node.y - 5);
                    text.setAttribute('text-anchor', 'middle');
                }

                text.setAttribute('fill', self.labelsConfig.bootstrapColor);
                text.setAttribute('font-size', `${self.labelsConfig.bootstrapSize * self.labelSizeScale}px`);
                text.setAttribute('font-family', self.labelsConfig.bootstrapFont);
                text.setAttribute('font-weight', self.labelsConfig.bootstrapWeight);
                text.textContent = node.name;
                gLabels.appendChild(text);
            }

            // 3. Branch Length Labels
            if (self.labelsConfig.showBranchLength && node.parent && node.branchLength !== null) {
                const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                const midX = (node.parent.x + node.x) / 2;
                const midY = (node.parent.y + node.y) / 2;

                if (!self.layout.startsWith('circular')) {
                    text.setAttribute('x', midX);
                    text.setAttribute('y', node.y - 5);
                    text.setAttribute('text-anchor', 'middle');
                } else {
                    const parentRad = node.parent.radius;
                    const childRad = node.radius;
                    const midRad = (parentRad + childRad) / 2;
                    const mx = width/2 + midRad * Math.cos(node.angle);
                    const my = height/2 + midRad * Math.sin(node.angle);
                    
                    text.setAttribute('x', mx);
                    text.setAttribute('y', my - 4);
                    text.setAttribute('text-anchor', 'middle');
                }

                text.setAttribute('fill', self.labelsConfig.branchLengthColor);
                text.setAttribute('font-size', `${self.labelsConfig.branchLengthSize * self.labelSizeScale}px`);
                text.setAttribute('font-family', self.labelsConfig.branchLengthFont);
                text.setAttribute('font-weight', self.labelsConfig.branchLengthWeight);
                text.textContent = node.branchLength.toFixed(4);
                gLabels.appendChild(text);
            }
        }

        draw(this.tree);

        // Event listener on SVG background to deselect when clicking outside nodes/branches
        svgElement.addEventListener('click', (e) => {
            if (e.target === svgElement) {
                self.selectNode(null);
            }
        });
    };

    // Helper: Describe circular arc path coordinates
    function describeArc(x, y, radius, startAngle, endAngle) {
        const startX = x + radius * Math.cos(startAngle);
        const startY = y + radius * Math.sin(startAngle);
        const endX = x + radius * Math.cos(endAngle);
        const endY = y + radius * Math.sin(endAngle);

        const largeArcFlag = endAngle - startAngle <= Math.PI ? '0' : '1';
        
        // Arc command: A rx ry x-axis-rotation large-arc-flag sweep-flag x y
        return `M ${startX} ${startY} A ${radius} ${radius} 0 ${largeArcFlag} 1 ${endX} ${endY}`;
    }

    // ─── Select/Highlight Node & Dispatch Event ───────────────────────────────
    MyTreesViz.selectNode = function(node) {
        this.selectedNode = node;
        
        // Redraw to show selection immediately
        const svg = document.getElementById('tree-svg');
        if (svg) this.render(svg);

        // Fire callback to notify main UI controller
        if (this.onNodeSelectedCallback) {
            this.onNodeSelectedCallback(node);
        }
    };

    // Update style properties for specific nodes
    MyTreesViz.updateNodeStyle = function(nodeId, properties) {
        let style = this.nodeStyles.get(nodeId);
        if (!style) {
            style = { stroke: '#e2f1f8', strokeWidth: 2, strokeDasharray: 'none' };
        }
        this.nodeStyles.set(nodeId, { ...style, ...properties });

        const svg = document.getElementById('tree-svg');
        if (svg) this.render(svg);
    };

    // Register globally
    window.MyTreesViz = MyTreesViz;

})(window);
