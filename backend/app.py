import os
import sys
import traceback

from flask import Flask, request, jsonify, send_from_directory
from flask_login import login_required

# Allow importing sibling modules (mytrees, models, auth) by absolute name
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)

# Load environment variables from a .env file at the repo root, if present
try:
    from dotenv import load_dotenv
    load_dotenv(os.path.join(os.path.dirname(BASE_DIR), '.env'))
except ImportError:
    pass

import mytrees
from models import db
from auth import auth_bp, init_auth

app = Flask(__name__)

# ─── Configuration ───────────────────────────────────────────────────────────
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-insecure-change-me')

# Database: default to a SQLite file inside backend/. A relative sqlite path is
# resolved against this directory so it works regardless of the launch CWD.
database_url = os.environ.get('DATABASE_URL', 'sqlite:///mytrees.db')
if database_url.startswith('sqlite:///'):
    rel = database_url[len('sqlite:///'):]
    if not os.path.isabs(rel):
        database_url = 'sqlite:///' + os.path.join(BASE_DIR, rel)
app.config['SQLALCHEMY_DATABASE_URI'] = database_url
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Session cookie hardening. SESSION_COOKIE_SECURE must stay off for local http
# development (a Secure cookie is not sent over http); enable it behind HTTPS.
app.config['SESSION_COOKIE_HTTPONLY'] = True
app.config['SESSION_COOKIE_SAMESITE'] = 'Lax'
app.config['SESSION_COOKIE_SECURE'] = os.environ.get('SESSION_COOKIE_SECURE', '0') == '1'

# ─── Extensions ──────────────────────────────────────────────────────────────
db.init_app(app)
init_auth(app)
app.register_blueprint(auth_bp)

with app.app_context():
    db.create_all()

# Resolve frontend absolute directory path
frontend_dir = os.path.abspath(os.path.join(BASE_DIR, '..', 'frontend'))

# Note: the frontend is served same-origin by this app, so no permissive CORS
# headers are needed (the previous Access-Control-Allow-Origin '*' was removed).

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
@login_required
def build_tree():
    if request.method == 'OPTIONS':
        return jsonify({'status': 'ok'})
        
    data = request.json or {}
    
    # 1. Inputs
    input_type = data.get('inputType', 'alignment')  # 'alignment' or 'matrix'
    method = data.get('method', 'NJ')              # NJ, UPGMA, WPGMA, FM, ME, MP, ML
    model = data.get('model', 'JC69')              # JC69, K2P, F84, LogDet, p-distance
    try:
        bootstrap = int(data.get('bootstrap', 0) or 0)
    except (TypeError, ValueError):
        bootstrap = 0

    # Distance-based methods that work on a distance matrix (and support bootstrap)
    DISTANCE_METHODS = ('NJ', 'UPGMA', 'WPGMA', 'FM', 'ME')

    def build_distance_tree(t, m):
        if method == 'UPGMA':
            return mytrees.upgma(t, m, weighted=False)
        if method == 'WPGMA':
            return mytrees.upgma(t, m, weighted=True)
        if method == 'FM':
            return mytrees.fitch_margoliash(t, m)
        if method == 'ME':
            return mytrees.minimum_evolution(t, m)
        return mytrees.neighbor_joining(t, m)
    
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
            if bootstrap > 0 and method in DISTANCE_METHODS:
                # Bootstrap resampling annotates internal nodes with support (%)
                tree = mytrees.bootstrap_support(seqs, bootstrap, method, model)
            elif method == 'MP':
                tree = mytrees.maximum_parsimony_tree(seqs)
            elif method == 'ML':
                tree = mytrees.maximum_likelihood_tree(seqs, model)
            else:
                tree = build_distance_tree(taxa, matrix)

            newick = tree.to_newick()
            stats = mytrees.tree_stats(tree)
            # Add likelihood/parsimony metadata if exists
            if 'log_likelihood' in tree.metadata:
                stats['log_likelihood'] = tree.metadata['log_likelihood']
            if 'parsimony_score' in tree.metadata:
                stats['parsimony_score'] = tree.metadata['parsimony_score']
            
            # Merge matrix stats
            stats.update(matrix_stats)
            # Empirical nucleotide frequencies (for the statistics view)
            try:
                stats['base_frequencies'] = mytrees.get_base_frequencies(list(seqs.values()))
            except Exception:
                pass

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
                
            # Build tree from matrix (MP & ML require alignments, fall back to NJ)
            if method in ('MP', 'ML'):
                method = 'NJ'  # Override to Neighbor Joining

            tree = build_distance_tree(taxa, matrix)

            newick = tree.to_newick()
            stats = mytrees.tree_stats(tree)
            
            # Analyze matrix stats
            matrix_stats = mytrees.analyze_matrix_stats(taxa, matrix)
            stats.update(matrix_stats)
            
        else:
            return jsonify({'error': 'Tipo de entrada não reconhecido.'}), 400
            
        # Surface bootstrap replicate count in stats when present
        if 'bootstrap_reps' in tree.metadata:
            stats['bootstrap_reps'] = tree.metadata['bootstrap_reps']

        return jsonify({
            # Structured tree preserves internal-node metadata (e.g. bootstrap
            # support) that a bare Newick string would lose.
            'tree': tree.to_dict(),
            'newick': newick + ";",
            'taxa': taxa,
            'matrix': matrix,
            'stats': stats,
            'method': method,
            'model': model
        })
        
    except Exception:
        # Log full details server-side; return a generic message to the client
        # so internal structure / stack traces are not leaked.
        traceback.print_exc()
        return jsonify({'error': 'Falha no cálculo filogenético. Verifique os dados de entrada.'}), 500

# ─── Static Files Router ─────────────────────────────────────────────────────

@app.route('/')
def serve_landing():
    # Public marketing landing page
    return send_from_directory(frontend_dir, 'index.html')

@app.route('/app')
def serve_app():
    # The phylogenetics workspace (SPA); gated by the login overlay
    return send_from_directory(frontend_dir, 'app.html')

@app.route('/<path:path>')
def serve_static(path):
    full = os.path.join(frontend_dir, path)
    if os.path.isfile(full):
        return send_from_directory(frontend_dir, path)
    # Unknown path → fall back to the landing page
    return send_from_directory(frontend_dir, 'index.html')

if __name__ == '__main__':
    # Make sure frontend dir exists
    os.makedirs(frontend_dir, exist_ok=True)
    debug = os.environ.get('FLASK_DEBUG', '0') == '1'
    port = int(os.environ.get('PORT', '5000'))
    print(f"Iniciando MyTrees em http://localhost:{port}  (debug={debug})")
    print(f"Servindo arquivos estáticos de: {frontend_dir}")
    app.run(host='0.0.0.0', port=port, debug=debug)
