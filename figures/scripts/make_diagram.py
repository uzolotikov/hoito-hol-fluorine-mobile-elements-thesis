#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Схема локусов: гены стрелками, цвет по функции, рибопереключатель отдельным знаком."""
import json, io

import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', 'data', 'diagram_data.json')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, '..', 'fig3_locus_maps.svg')
D = json.load(io.open(DATA, encoding='utf-8'))

# функция -> (цвет, подпись)
def classify(dm):
    s = set(dm)
    if s & {'CRCB'}:                       return ('#C44E52', 'CrcB (экспортёр F⁻)')
    if s & {'Voltage_CLC'}:                return ('#DD8452', 'CLCᶠ (экспортёр F⁻)')
    if s & {'Enolase_N','Enolase_C'}:      return ('#4C72B0', 'енолаза')
    if s & {'Pyrophosphatase','PPase'}:    return ('#17becf', 'пирофосфатаза')
    if s & {'Rep_3','Rep_1','Rep_2','RepA_N','RepA_C','Rep_trans'}:
                                           return ('#55A868', 'репликация плазмиды')
    if s & {'MobA_MobL','MobC','Relaxase','MobL'}:
                                           return ('#8172B3', 'мобилизация')
    if s & {'T4SS-DNA_transf'}:            return ('#937860', 'конъюгация T4SS')
    if s & {'Resolvase','Recombinase','SpecificRecomb'}:
                                           return ('#DA8BC3', 'резольваза/рекомбиназа')
    if s & {'ParBc','ParB'}:               return ('#64B5CD', 'разделение (ParB)')
    if s & {'DUF190'}:                     return ('#CCB974', 'DUF190')
    if s:                                  return ('#BBBBBB', 'прочее')
    return ('#E8E8E8', 'без домена')

ROW_H, GENE_H, PAD_L, PAD_R, W_TRACK = 108, 22, 300, 40, 900
Y0 = 100
H = 86 + ROW_H * len(D) + 150
W = PAD_L + W_TRACK + PAD_R

def esc(s): return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')

svg = []
svg.append('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" font-family="DejaVu Sans, Arial, sans-serif">' % (W,H,W,H))
svg.append('<rect width="100%" height="100%" fill="#ffffff"/>')
svg.append('<text x="%d" y="34" font-size="19" font-weight="bold" fill="#222">Фторидные локусы Хойто-Гола: генное окружение</text>' % (24,))
svg.append('<text x="24" y="54" font-size="12.5" fill="#666">Гены — стрелки по направлению транскрипции; ▼ — фторидный рибопереключатель RF01734; числа под стрелками — межгенные промежутки, н.</text>')
svg.append('<text x="24" y="70" font-size="11.5" fill="#666">Первые две панели — один и тот же элемент в двух вариантах сборки; он собран в противоположной ориентации, поэтому схемы зеркальны.</text>')

legend_used = {}
y = 100
for panel in D:
    a, b = panel['view']; span = max(1, b - a)
    sc = lambda x: PAD_L + (x - a) * W_TRACK / span

    svg.append('<text x="24" y="%d" font-size="13" font-weight="bold" fill="#222">%s</text>' % (y+4, esc(panel['title'])))
    svg.append('<text x="24" y="%d" font-size="11" fill="#777">%s · %s н.%s</text>'
               % (y+20, esc(panel['contig'][:34]), format(panel['length'],',').replace(',',' '),
                  '' if panel['view'][0]==1 and panel['view'][1]==panel['length'] else '  (показан участок %s–%s)' % (format(a,',').replace(',',' '), format(b,',').replace(',',' '))))

    ymid = y + 44
    svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#999" stroke-width="1.6"/>' % (sc(a), ymid, sc(b), ymid))

    genes = sorted(panel['genes'], key=lambda g: g['start'])
    for i, g in enumerate(genes):
        col, lab = classify(g['dom']); legend_used[lab] = col
        x1, x2 = sc(g['start']), sc(g['stop'])
        if x2 - x1 < 3: x2 = x1 + 3
        tip = min(11.0, max(3.0, (x2-x1)*0.32))
        top, bot = ymid - GENE_H/2, ymid + GENE_H/2
        if g['strand'] == '1':
            pts = "%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f" % (
                x1,top, x2-tip,top, x2,ymid, x2-tip,bot, x1,bot)
        else:
            pts = "%.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f %.1f,%.1f" % (
                x2,top, x1+tip,top, x1,ymid, x1+tip,bot, x2,bot)
        svg.append('<polygon points="%s" fill="%s" stroke="#333" stroke-width="0.8"/>' % (pts, col))
        if x2 - x1 > 26:
            svg.append('<text x="%.1f" y="%.1f" font-size="9.5" text-anchor="middle" fill="#fff">%d</text>'
                       % ((x1+x2)/2, ymid+3.5, g['n']))
        # межгенный промежуток
        if i+1 < len(genes):
            nxt = genes[i+1]; gap = nxt['start'] - g['stop'] - 1
            xm = (sc(g['stop']) + sc(nxt['start'])) / 2
            svg.append('<text x="%.1f" y="%.1f" font-size="8.5" text-anchor="middle" fill="#888">%d</text>' % (xm, bot+11, gap))

    for r in panel['riboswitches']:
        xr = sc((r['lo']+r['hi'])/2.0)
        svg.append('<polygon points="%.1f,%.1f %.1f,%.1f %.1f,%.1f" fill="#222"/>' % (xr-6, ymid-24, xr+6, ymid-24, xr, ymid-13))
        svg.append('<text x="%.1f" y="%.1f" font-size="9" text-anchor="middle" fill="#444">RF01734</text>' % (xr, ymid-28))
    y += ROW_H

# легенда
ly = y + 6
svg.append('<text x="24" y="%d" font-size="12" font-weight="bold" fill="#333">Обозначения</text>' % ly)
ly += 16
cx, n = 24, 0
order = ['CrcB (экспортёр F⁻)','CLCᶠ (экспортёр F⁻)','енолаза','пирофосфатаза',
         'репликация плазмиды','мобилизация','конъюгация T4SS','резольваза/рекомбиназа',
         'разделение (ParB)','DUF190','прочее','без домена']
for lab in order:
    if lab not in legend_used: continue
    svg.append('<rect x="%d" y="%d" width="13" height="13" fill="%s" stroke="#333" stroke-width="0.7"/>' % (cx, ly-10, legend_used[lab]))
    svg.append('<text x="%d" y="%d" font-size="11" fill="#333">%s</text>' % (cx+18, ly, esc(lab)))
    cx += 210; n += 1
    if n % 4 == 0: cx = 24; ly += 20
svg.append('</svg>')

io.open(OUT,'w',encoding='utf-8').write('\n'.join(svg))
print("готово: fig3_locus_maps.svg  (%d панелей)" % len(D))
