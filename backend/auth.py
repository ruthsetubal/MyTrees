"""Rotas de autenticação e configuração do Flask-Login.

Expõe um Blueprint com /api/register, /api/login, /api/logout e /api/me.
Por ser uma SPA, respostas não-autenticadas retornam 401 em JSON em vez de
redirecionar para uma página de login.
"""

import re

from flask import Blueprint, request, jsonify
from flask_login import (
    LoginManager, login_user, logout_user, login_required, current_user
)

from models import db, User

auth_bp = Blueprint('auth', __name__)
login_manager = LoginManager()

EMAIL_RE = re.compile(r'^[^@\s]+@[^@\s]+\.[^@\s]+$')
USERNAME_RE = re.compile(r'^[A-Za-z0-9_.-]{3,32}$')
MIN_PASSWORD_LEN = 8


def init_auth(app):
    """Conecta o Flask-Login à aplicação."""
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        return db.session.get(User, int(user_id))

    @login_manager.unauthorized_handler
    def unauthorized():
        return jsonify({'error': 'Autenticação necessária.'}), 401


@auth_bp.route('/api/register', methods=['POST'])
def register():
    data = request.get_json(silent=True) or {}
    username = (data.get('username') or '').strip()
    email = (data.get('email') or '').strip().lower()
    password = data.get('password') or ''

    if not USERNAME_RE.match(username):
        return jsonify({'error': 'Usuário deve ter 3-32 caracteres (letras, números, . _ -).'}), 400
    if not EMAIL_RE.match(email):
        return jsonify({'error': 'E-mail inválido.'}), 400
    if len(password) < MIN_PASSWORD_LEN:
        return jsonify({'error': f'A senha deve ter pelo menos {MIN_PASSWORD_LEN} caracteres.'}), 400

    if User.query.filter_by(username=username).first():
        return jsonify({'error': 'Nome de usuário já cadastrado.'}), 409
    if User.query.filter_by(email=email).first():
        return jsonify({'error': 'E-mail já cadastrado.'}), 409

    user = User(username=username, email=email)
    user.set_password(password)
    db.session.add(user)
    db.session.commit()

    login_user(user)
    return jsonify({'user': user.to_dict()}), 201


@auth_bp.route('/api/login', methods=['POST'])
def login():
    data = request.get_json(silent=True) or {}
    identifier = (data.get('username') or data.get('email') or '').strip()
    password = data.get('password') or ''

    user = User.query.filter(
        (User.username == identifier) | (User.email == identifier.lower())
    ).first()

    if user is None or not user.check_password(password):
        # Mensagem genérica: não revela se o usuário existe.
        return jsonify({'error': 'Usuário ou senha inválidos.'}), 401

    login_user(user, remember=bool(data.get('remember', False)))
    return jsonify({'user': user.to_dict()})


@auth_bp.route('/api/logout', methods=['POST'])
@login_required
def logout():
    logout_user()
    return jsonify({'status': 'ok'})


@auth_bp.route('/api/me', methods=['GET'])
def me():
    if current_user.is_authenticated:
        return jsonify({'user': current_user.to_dict()})
    return jsonify({'user': None}), 401
