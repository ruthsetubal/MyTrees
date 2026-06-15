# MyTrees

Ferramenta web para análise filogenética. Você carrega um alinhamento de sequências
(ou uma matriz de distâncias) e o MyTrees constrói e desenha a árvore evolutiva, com
vários métodos de reconstrução e modelos de distância.

O cálculo roda em um motor próprio em Python (Flask) no servidor, e também há uma
versão em JavaScript que roda direto no navegador como alternativa offline.

## Funcionalidades

- Entrada em FASTA, PHYLIP, CSV e Newick (o tipo é detectado automaticamente)
- Cinco modelos de distância: p-distância, JC69, K2P, F84 e LogDet
- Sete métodos de reconstrução: Neighbor Joining, UPGMA, WPGMA, Fitch-Margoliash,
  Evolução Mínima, Máxima Parcimônia e Máxima Verossimilhança
- Suporte estatístico por bootstrap (100/500/1000 réplicas)
- Visualização interativa: filograma e cladograma, retangular e circular; reenraizar,
  rotacionar, inverter, colapsar clados, ponto médio, zoom/escala
- Customização num painel flutuante: cor de fundo, rótulos, estilo por ramo e por clado,
  marcadores nas pontas e barra de escala; desfazer com Ctrl+Z e botão de reset
- Matriz de distâncias (heatmap) e de caracteres, calculadas a partir de sequências
- Estatísticas: frequências de bases, métricas da árvore e suporte dos ramos
- Exportação: árvore em SVG, PNG, Newick e NEXUS; matriz em CSV, PHYLIP e NEXUS
- Contas de usuário (login/registro) e documentação dentro do app

## Na aplicação

A aplicação tem três abas:

- **Árvore** — entrada de dados, construção e visualização. A customização fica num
  painel que abre pelo botão "Customização".
- **Matriz** — calcula a matriz de distâncias a partir de sequências (ou exibe uma matriz
  já pronta) e também a matriz de caracteres (o alinhamento).
- **Estatísticas** — métricas da árvore, frequências de bases e suporte dos ramos.

O botão "Documentação" no topo abre um guia rápido de uso.

## Stack

- Backend: Python 3 + Flask, Flask-Login, Flask-SQLAlchemy (SQLite)
- Frontend: HTML, CSS e JavaScript puro (sem framework)

## Como rodar

Requer Python 3.10+.

```bash
python -m venv venv
# Windows:
venv\Scripts\activate
# Linux/Mac:
# source venv/bin/activate

pip install -r backend/requirements.txt

# configure as variáveis de ambiente
cp .env.example .env
# edite o .env e defina uma SECRET_KEY (ex.: python -c "import secrets; print(secrets.token_hex(32))")

python backend/app.py
```

Depois abra http://localhost:5000.

- `/` — página inicial
- `/app` — aplicação (pede login; crie uma conta na primeira vez)

A porta pode ser trocada com a variável `PORT`. Em desenvolvimento, use `FLASK_DEBUG=1`.

## Estrutura

```
backend/
  app.py        # servidor Flask, rotas da API e arquivos estáticos
  mytrees.py    # motor filogenético (distâncias, métodos, bootstrap)
  models.py     # modelo de usuário (SQLAlchemy)
  auth.py       # rotas de login/registro
frontend/
  index.html    # landing page
  app.html      # aplicação
  mytrees.js    # motor filogenético em JS (modo offline)
  visualization.js  # renderização SVG da árvore
  styles.css
  assets/       # logo e imagens
```

## Formatos de entrada

FASTA:

```
>Humano
ATGCTAGGGTTCCTATGTTTGGTG
>Chimpanze
ATGCTAGGGTTCCTATGTTTAGTG
```

Matriz de distâncias em CSV (com ou sem cabeçalho):

```
,Humano,Chimpanze,Gorila
Humano,0,0.082,0.124
Chimpanze,0.082,0,0.118
Gorila,0.124,0.118,0
```

Também aceita matriz/alinhamento em PHYLIP e árvores em Newick.
