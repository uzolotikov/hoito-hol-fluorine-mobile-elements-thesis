#!/usr/bin/env python3
"""Плазмиды среди контигов с локусами в одной сборке.

Контиг считается плазмидным, если он попал в plasmid_summary geNomad (пороги по умолчанию).
Ожидаемое число: для каждого контига с локусом берётся доля плазмидных среди его контрольных
контигов той же сборки и длины +-25 %; сумма этих долей и есть ожидание, а число плазмидных
контигов сравнивается с пуассон-биномиальным распределением с этими вероятностями.

Запуск: workflow/04_plasmids.py <каталог сборки в OUT>
"""
import io, os, re, sys
from collections import defaultdict

D = sys.argv[1]
NAME = os.path.basename(os.path.normpath(D))


def path(*f):
    return os.path.join(D, *f)


def table(fn):
    rows = list(io.open(fn, encoding="utf-8"))
    hdr = rows[0].rstrip("\n").split("\t")
    return [dict(zip(hdr, r.rstrip("\n").split("\t"))) for r in rows[1:] if r.strip()]


summary = {r["seq_name"]: r for r in table(path("genomad_out", "genomad_input_summary", "genomad_input_plasmid_summary.tsv"))}
scores = {r["seq_name"]: r for r in table(path("genomad_out", "genomad_input_aggregated_classification",
                                               "genomad_input_aggregated_classification.tsv"))}
clen = dict(l.rstrip("\n").split("\t") for l in io.open(path("contig_lengths.tsv"), encoding="utf-8"))
controls = {}
for l in io.open(path("controls.tsv"), encoding="utf-8"):
    a, b = l.rstrip("\n").split("\t")
    controls[a] = [x for x in b.split(",") if x]
markers = defaultdict(set)
for l in io.open(path("plasmid_markers.tbl"), encoding="utf-8"):
    if l[0] != "#" and l.strip():
        p = l.split()
        markers[re.sub(r"_\d+$", "", p[0])].add(p[2])
loci = table(path("loci.tsv"))


def poisson_binomial_ge(k, probs):
    dp = [1.0] + [0.0] * len(probs)
    for p in probs:
        for j in range(len(dp) - 1, 0, -1):
            dp[j] = dp[j] * (1 - p) + dp[j - 1] * p
        dp[0] *= 1 - p
    return sum(dp[k:])


contigs = sorted({r["contig"] for r in loci})
probs = [sum(y in summary for y in controls[c]) / len(controls[c]) if controls.get(c) else 0.0 for c in contigs]
k = sum(c in summary for c in contigs)
pval = poisson_binomial_ge(k, probs)
all_ctrl = sorted({y for v in controls.values() for y in v})
n_ctrl_pl = sum(y in summary for y in all_ctrl)
with io.open(path("plasmid_test.tsv"), "w", encoding="utf-8") as f:
    f.write("assembly\tlocus_contigs\tplasmid_observed\tplasmid_expected\tP_ge\tcontrol_contigs\tcontrol_plasmid\n")
    f.write("%s\t%d\t%d\t%.2f\t%.3g\t%d\t%d\n" % (NAME, len(contigs), k, sum(probs), pval, len(all_ctrl), n_ctrl_pl))
print("%s: контигов с локусами %d, плазмидных %d, ожидалось %.2f, P(X>=%d) = %.3g; в контроле %d из %d" % (
    NAME, len(contigs), k, sum(probs), k, pval, n_ctrl_pl, len(all_ctrl)))

with io.open(path("locus_contigs.tsv"), "w", encoding="utf-8") as f:
    f.write("contig\tcontig_len\tloci\tgenomad_chromosome\tgenomad_plasmid\tgenomad_virus\tgenomad_plasmid_call"
            "\tconjugation_genes\tplasmid_marker_domains\n")
    for c in contigs:
        s = scores.get(c, {})
        conj = summary.get(c, {}).get("conjugation_genes", "")
        f.write("\t".join([c, clen[c], ",".join(r["exporter"] for r in loci if r["contig"] == c),
                           s.get("chromosome_score", "NA"), s.get("plasmid_score", "NA"), s.get("virus_score", "NA"),
                           "yes" if c in summary else "no", conj if conj not in ("", "NA") else "",
                           ",".join(sorted(markers.get(c, [])))]) + "\n")
        if c in summary:
            print("  плазмидный: %s, длина %s, маркеры: %s" % (c, clen[c], ",".join(sorted(markers.get(c, []))) or "нет"))
