#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fig2 по трём сборкам: транспозазы (ожидание против наблюдения) и плазмиды (ожидание против наблюдения)."""
import sys, json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.family'] = 'DejaVu Sans'
import os
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..') + os.sep
ARG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'data', 'fig2_data.json')
D = json.loads(open(ARG, encoding='utf-8').read()) if os.path.exists(ARG) else json.loads(ARG)
asm = ['короткая', 'гибридная', 'Flye']
c = lambda x, f='%.1f': (f % x).replace('.', ',')
def pfmt(p):
    if p >= 0.01: return c(p, '%.2f')
    m, e = ('%.1e' % p).split('e'); return '$%s\\cdot10^{%d}$' % (m.replace('.', '{,}'), int(e))
fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4))
x = np.arange(3); w = 0.36
ax = axes[0]
e, o = D['tnp_exp'], D['tnp_obs']
b1 = ax.bar(x - w/2, e, w, color='#B0B0B0', edgecolor='black', label='ожидалось по фону своей сборки')
b2 = ax.bar(x + w/2, o, w, color='#4C72B0', edgecolor='black', label='наблюдалось')
for i in range(3):
    ax.text(x[i] - w/2, e[i] + 0.15, c(e[i]), ha='center', va='bottom', fontsize=9)
    ax.text(x[i] + w/2, o[i] + 0.15, str(o[i]), ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.text(x[i], -0.15, 'P(X≥набл.) = ' + pfmt(D['tnp_p'][i]), ha='center', va='top', fontsize=8, color='#444', transform=ax.get_xaxis_transform())
ax.set_xticks(x); ax.set_xticklabels(asm, fontsize=10)
ax.set_ylabel('Белков с доменами МГЭ в окнах ±10 генов', fontsize=9.5)
ax.set_title('Транспозазы: обогащения нет', fontsize=11)
ax.set_ylim(0, max(e + o) * 1.35); ax.legend(fontsize=8.3, loc='upper left', frameon=True)
ax.grid(axis='y', linestyle='--', alpha=0.45); ax.set_axisbelow(True)
ax = axes[1]
e, o, n = D['pl_exp'], D['pl_obs'], D['pl_n']
ax.bar(x - w/2, e, w, color='#B0B0B0', edgecolor='black', label='ожидалось по контролю той же длины')
ax.bar(x + w/2, o, w, color='#C44E52', edgecolor='black', label='наблюдалось')
for i in range(3):
    ax.text(x[i] - w/2, e[i] + 0.06, c(e[i], '%.2f'), ha='center', va='bottom', fontsize=9)
    ax.text(x[i] + w/2, o[i] + 0.06, '%d из %d' % (o[i], n[i]), ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.text(x[i], -0.15, 'p = ' + pfmt(D['pl_p'][i]), ha='center', va='top', fontsize=8, color='#444', transform=ax.get_xaxis_transform())
ax.set_xticks(x); ax.set_xticklabels(asm, fontsize=10)
ax.set_ylabel('Контигов с локусами, отнесённых к плазмидам', fontsize=9.5)
ax.set_title('Плазмиды: обогащение есть', fontsize=11)
ax.set_ylim(0, max(o) * 1.45 + 0.5); ax.legend(fontsize=8.3, loc='upper right', frameon=True)
ax.grid(axis='y', linestyle='--', alpha=0.45); ax.set_axisbelow(True)
fig.suptitle('Фторидные локусы Хойто-Гола: связь с плазмидами, но не с транспозазами', fontsize=12, y=1.02)
fig.text(0.5, -0.07, 'Окна вокруг генов CrcB и клады CLCᶠ; фон и контроль — из той же сборки с учётом длины контигов. '
         'Плазмиды — geNomad.', ha='center', fontsize=8.3, color='#444')
plt.tight_layout()
plt.savefig(OUT + 'fig2_transposase_enrichment.svg', bbox_inches='tight')
plt.savefig(OUT + 'fig2_transposase_enrichment.pdf', bbox_inches='tight')
print('fig2 готов')
