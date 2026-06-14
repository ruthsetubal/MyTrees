import os
import sys
import json
from flask import Flask, request, jsonify, send_from_directory

# Add parent directory to path so mytrees can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import mytrees

app = Flask(__name__)

# Resolve frontend absolute directory path
frontend_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'frontend'))

# Simple CORS decorator for Flask endpoints
def add_cors_headers(response):
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Methods'] = 'POST, GET, OPTIONS, PUT, DELETE'
    response.headers['Access-Control-Allow-Headers'] = 'Content-Type, X-Requested-With'
    return response

@app.after_request
def after_request_func(response):
    return add_cors_headers(response)

@app.route('/api/test', methods=['GET', 'OPTIONS'])
def test_connection():
    if request.method == 'OPTIONS':
        return jsonify({'status': 'ok'})
    return jsonify({
        'status': 'connected',
        'message': 'MyTrees Backend Server is Active!',
        'engine': 'Python v3 (Flask)',
        'capabilities': ['p-distance', 'JC69', 'K2P', 'F84', 'LogDet', 'UPGMA', 'WPGMA', 'NJ', 'Maximum Parsimony', 'Maximum Likelihood']
    })

@app.route('/api/build_tree', methods=['POST', 'OPTIONS'])
def build_tree():
    if request.method == 'OPTIONS':
        return jsonify({'status': 'ok'})
        
    data = request.json or {}
    
    # 1. Inputs
    input_type = data.get('inputType', 'alignment')  # 'alignment' or 'matrix'
    method = data.get('method', 'NJ')              # NJ, UPGMA, WPGMA, MP, ML
    model = data.get('model', 'JC69')              # JC69, K2P, F84, LogDet, p-distance
    
    # 2. Logic based on input type
    taxa = []
    matrix = []
    newick = ""
    stats = {}
    
    try:
        if input_type == 'alignment':
            # Receive sequence map
            seqs_raw = data.get('seqs', {})
            # Parse alignments if text is passed instead
            raw_text = data.get('alignmentText', '')
            
            if raw_text.strip():
                # Heuristic detect FASTA or PHYLIP Sequential/Interleaved
                raw_text_strip = raw_text.strip()
                if raw_text_strip.startswith('>'):
                    seqs = mytrees.parse_fasta(raw_text)
                else:
                    seqs = mytrees.parse_phylip_alignment(raw_text)
            else:
                seqs = seqs_raw
                
            if not seqs:
                return jsonify({'error': 'Nenhuma sequência biológica válida carregada.'}), 400
                
            # Clean seqs keys and values
            seqs = {k.strip(): v.replace(" ", "").upper() for k, v in seqs.items() if k.strip()}
            
            # Compute distance matrix
            taxa, matrix, matrix_stats = mytrees.compute_distance_matrix(seqs, model)
            
            # Build tree
            if method == 'NJ':
                tree = mytrees.neighbor_joining(taxa, matrix)
            elif method == 'UPGMA':
                tree = mytrees.upgma(taxa, matrix, weighted=False)
            elif method == 'WPGMA':
                tree = mytrees.upgma(taxa, matrix, weighted=True)
            elif method == 'MP':
                tree = mytrees.maximum_parsimony_tree(seqs)
            elif method == 'ML':
                tree = mytrees.maximum_likelihood_tree(seqs, model)
            else:
                tree = mytrees.neighbor_joining(taxa, matrix)
                
            newick = tree.to_newick()
            stats = mytrees.tree_stats(tree)
            # Add likelihood/parsimony metadata if exists
            if 'log_likelihood' in tree.metadata:
                stats['log_likelihood'] = tree.metadata['log_likelihood']
            if 'parsimony_score' in tree.metadata:
                stats['parsimony_score'] = tree.metadata['parsimony_score']
            
            # Merge matrix stats
            stats.update(matrix_stats)
            
        elif input_type == 'matrix':
            # Receive pre-computed matrix
            taxa = data.get('taxa', [])
            matrix = data.get('matrix', [])
            
            # Parse matrix from text if passed
            raw_matrix_text = data.get('matrixText', '')
            if raw_matrix_text.strip():
                taxa, matrix = mytrees.parse_phylip_matrix(raw_matrix_text)
                
            if not taxa or not matrix:
                return jsonify({'error': 'Matriz de distância ou Táxons vazios/inválidos.'}), 400
                
            # Build tree from matrix (MP & ML require alignments, fallback to NJ if selected)
            if method in ('MP', 'ML'):
                method = 'NJ'  # Override to Neighbor Joining
                
            if method == 'NJ':
                tree = mytrees.neighbor_joining(taxa, matrix)
            elif method == 'UPGMA':
                tree = mytrees.upgma(taxa, matrix, weighted=False)
            elif method == 'WPGMA':
                tree = mytrees.upgma(taxa, matrix, weighted=True)
            else:
                tree = mytrees.neighbor_joining(taxa, matrix)
                
            newick = tree.to_newick()
            stats = mytrees.tree_stats(tree)
            
            # Analyze matrix stats
            matrix_stats = mytrees.analyze_matrix_stats(taxa, matrix)
            stats.update(matrix_stats)
            
        else:
            return jsonify({'error': 'Tipo de entrada não reconhecido.'}), 400
            
        return jsonify({
            'newick': newick + ";",
            'taxa': taxa,
            'matrix': matrix,
            'stats': stats,
            'method': method,
            'model': model
        })
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': f'Falha no cálculo filogenético: {str(e)}'}), 500

# ─── Static Files Router ─────────────────────────────────────────────────────

@app.route('/')
def serve_index():
    return send_from_directory(frontend_dir, 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    if not os.path.exists(os.path.join(frontend_dir, path)):
        return send_from_directory(frontend_dir, 'index.html')
    return send_from_directory(frontend_dir, path)

if __name__ == '__main__':
    # Make sure frontend dir exists
    os.makedirs(frontend_dir, exist_ok=True)
    print(f"Iniciando MyTrees no endereço http://localhost:5000")
    print(f"Servindo arquivos estáticos de: {frontend_dir}")
    app.run(host='0.0.0.0', port=5000, debug=True)
