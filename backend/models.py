"""Modelos de banco de dados do MyTrees (SQLAlchemy)."""

from datetime import datetime, timezone

from flask_sqlalchemy import SQLAlchemy
from flask_login import UserMixin
from werkzeug.security import generate_password_hash, check_password_hash

db = SQLAlchemy()


class User(UserMixin, db.Model):
    """Usuário da aplicação. A senha nunca é guardada em texto puro."""

    __tablename__ = 'users'

    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), unique=True, nullable=False, index=True)
    email = db.Column(db.String(255), unique=True, nullable=False, index=True)
    password_hash = db.Column(db.String(255), nullable=False)
    created_at = db.Column(db.DateTime, default=lambda: datetime.now(timezone.utc))

    def set_password(self, password: str) -> None:
        self.password_hash = generate_password_hash(password)

    def check_password(self, password: str) -> bool:
        return check_password_hash(self.password_hash, password)

    def to_dict(self) -> dict:
        """Representação segura (sem o hash da senha) para enviar ao frontend."""
        return {
            'id': self.id,
            'username': self.username,
            'email': self.email,
            'created_at': self.created_at.isoformat() if self.created_at else None,
        }

    def __repr__(self) -> str:
        return f'<User {self.username}>'


# ─── Futuro: histórico de análises por usuário ───────────────────────────────
# Quando quiser que cada conta salve suas árvores, descomente este modelo,
# rode db.create_all() novamente (ou use Flask-Migrate) e crie as rotas.
#
# class Analysis(db.Model):
#     __tablename__ = 'analyses'
#     id = db.Column(db.Integer, primary_key=True)
#     user_id = db.Column(db.Integer, db.ForeignKey('users.id'), nullable=False, index=True)
#     name = db.Column(db.String(120))
#     method = db.Column(db.String(16))
#     model = db.Column(db.String(16))
#     newick = db.Column(db.Text)
#     created_at = db.Column(db.DateTime, default=lambda: datetime.now(timezone.utc))
#     user = db.relationship('User', backref=db.backref('analyses', lazy=True))
