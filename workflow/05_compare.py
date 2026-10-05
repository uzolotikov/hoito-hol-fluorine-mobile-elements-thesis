#!/usr/bin/env python3
"""Итоговые таблицы по всем сборкам.

Участок локуса (рибопереключатель и ген экспортёра) с флангами по 300 п. н. ищется blastn в другой
сборке. Локусы двух сборок считаются одним и тем же, если выравнивание с идентичностью >= 99 %
покрывает >= 90 % самого локуса без флангов и перекрывает локус другой сборки. Фланги в критерий не
входят: у краёв контигов сборки расходятся. Связанные так локусы объединяются в группы; группа и есть
«различный локус». Если такое выравнивание есть, но локуса на этом месте нет, в таблице стоит
«участок без локуса»; если выравнивание покрывает от 50 до 90 % локуса, стоит «частично».

Запуск: workflow/05_compare.py <OUT> <имя сборки> [<имя сборки> ...]
"""
import io, os, statistics, sys
from collections import defaultdict

OUT, NAMES = sys.argv[1], sys.argv[2:]
PAD = 300  # фланги участка, как в 02_loci.py
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(OUT, "results")
os.makedirs(RES, exist_ok=True)


def table(fn):
    rows = list(io.open(fn, encoding="utf-8"))
    hdr = rows[0].rstrip("\n").split("\t")
    return [dict(zip(hdr, r.rstrip("\n").split("\t"))) for r in rows[1:] if r.strip()]


loci = {}
for a in NAMES:
    for r in table(os.path.join(OUT, a, "loci.tsv")):
        lid = "%s:%s" % (a, r["orf"])
        loci[lid] = dict(id=lid, asm=a, contig=r["contig"], L=int(r["contig_len"]), kind=r["exporter"],
                         lo=min(int(r["ribo_start"]), int(r["orf_start"])), hi=max(int(r["ribo_end"]), int(r["orf_end"])))
contig_info = {}
for a in NAMES:
    for r in table(os.path.join(OUT, a, "locus_contigs.tsv")):
        contig_info[(a, r["contig"])] = r
taxon = {}
tx = os.path.join(ROOT, "resources", "locus_taxonomy.tsv")
if os.path.exists(tx):
    for r in table(tx):
        taxon[(r["assembly"], r["contig"])] = r["taxon"]

edges = defaultdict(set)
status = {}  # (локус, сборка-мишень) -> статус, если связи нет
RANK = {"участок без локуса": 3, "частично": 2, "нет": 1}
for a in NAMES:
    for b in NAMES:
        if a == b:
            continue
        hits = defaultdict(list)
        for l in io.open(os.path.join(OUT, "compare", "%s_vs_%s.tsv" % (a, b)), encoding="utf-8"):
            p = l.rstrip("\n").split("\t")
            s1, s2 = int(p[6]), int(p[7])
            hits[p[0]].append(dict(s=p[1], pid=float(p[2]), qs=int(p[4]), qe=int(p[5]),
                                   ss=min(s1, s2), se=max(s1, s2), qlen=int(p[8])))
        targets = [x for x in loci.values() if x["asm"] == b]
        for q, hs in hits.items():
            lid = "%s:%s" % (a, q.split("@")[2])
            lo0, hi0 = map(int, q.split("@")[3].split("-"))
            cs = lo0 - max(1, lo0 - PAD) + 1          # начало локуса внутри участка
            ce = cs + hi0 - lo0                        # конец локуса внутри участка

            def core_cov(h):
                return max(0, min(h["qe"], ce) - max(h["qs"], cs) + 1) / (ce - cs + 1)

            good = [h for h in hs if h["pid"] >= 99.0 and core_cov(h) >= 0.9]
            linked = False
            for h in good:
                for t in targets:
                    if t["contig"] == h["s"] and t["lo"] <= h["se"] and h["ss"] <= t["hi"]:
                        edges[lid].add(t["id"]); edges[t["id"]].add(lid); linked = True
            if not linked:
                st = "участок без локуса" if good else (
                    "частично" if any(core_cov(h) >= 0.5 for h in hs) else "нет")
                status[(lid, b)] = st

seen, groups = set(), []
for lid in sorted(loci):
    if lid in seen:
        continue
    stack, comp = [lid], []
    while stack:
        v = stack.pop()
        if v in seen:
            continue
        seen.add(v); comp.append(v); stack += list(edges[v])
    groups.append(sorted(comp))


def group_key(g):
    for a in NAMES:
        m = [loci[x] for x in g if loci[x]["asm"] == a]
        if m:
            return (NAMES.index(a), m[0]["contig"], m[0]["lo"])


groups.sort(key=group_key)
pattern = defaultdict(int)
pairs = defaultdict(list)
with io.open(os.path.join(RES, "loci_across_assemblies.tsv"), "w", encoding="utf-8") as f:
    hdr = ["locus", "exporter", "taxon"]
    for a in NAMES:
        hdr += ["%s_contig" % a, "%s_contig_len" % a, "%s_plasmid_call" % a, "%s_plasmid_score" % a, "%s_plasmid_markers" % a]
    f.write("\t".join(hdr) + "\n")
    for i, g in enumerate(groups, 1):
        members = [loci[x] for x in g]
        kinds = sorted({m["kind"] for m in members})
        # таксон берётся из MAG; при расхождении сборок приоритет у гибридной
        ordered = sorted(members, key=lambda m: (m["asm"] != "hybrid", m["asm"]))
        tax = next((taxon[(m["asm"], m["contig"])] for m in ordered if (m["asm"], m["contig"]) in taxon), "NA")
        row = [str(i), "/".join(kinds), tax]
        present = ""
        for a in NAMES:
            m = [x for x in members if x["asm"] == a]
            if m:
                present += a[0].upper()
                ci = contig_info.get((a, m[0]["contig"]), {})
                row += [m[0]["contig"], str(m[0]["L"]), ci.get("genomad_plasmid_call", "NA"),
                        ci.get("genomad_plasmid", "NA"), ci.get("plasmid_marker_domains", "")]
            else:
                st = max((status.get((x["id"], a), "нет") for x in members), key=lambda s: RANK[s])
                row += [st, "NA", "NA", "NA", ""]
        pattern[present] += 1
        for a in NAMES:
            for b in NAMES:
                ma = [x for x in members if x["asm"] == a]; mb = [x for x in members if x["asm"] == b]
                if a < b and ma and mb:
                    pairs[(a, b)].append((ma[0]["L"], mb[0]["L"]))
        f.write("\t".join(row) + "\n")

print("различных локусов: %d" % len(groups))
print("присутствие по сборкам (первые буквы имён):", dict(sorted(pattern.items())))
for (a, b), v in sorted(pairs.items()):
    print("общих локусов %s и %s: %d; медиана длины контига %d и %d п. н." % (
        a, b, len(v), statistics.median(x for x, _ in v), statistics.median(y for _, y in v)))

with io.open(os.path.join(RES, "enrichment_tests.tsv"), "w", encoding="utf-8") as f:
    f.write("assembly\tgenes_CrcB_CLCF\twindow_orfs\tmge_observed\tmge_expected\tmge_P_ge\tmge_P_le"
            "\tlocus_contigs\tplasmid_observed\tplasmid_expected\tplasmid_P_ge\n")
    for a in NAMES:
        t = next(r for r in table(os.path.join(OUT, a, "transposase_test.tsv")) if r["set"] == "CrcB+CLC^F")
        p = table(os.path.join(OUT, a, "plasmid_test.tsv"))[0]
        f.write("\t".join([a, t["genes"], t["window_orfs"], t["mge_observed"], t["mge_expected"], t["P_ge"], t["P_le"],
                           p["locus_contigs"], p["plasmid_observed"], p["plasmid_expected"], p["P_ge"]]) + "\n")

with io.open(os.path.join(RES, "assembly_stats.tsv"), "w", encoding="utf-8") as f:
    f.write("assembly\tcontigs_500\ttotal_bp_500\tN50\tmax_len\tGC\tcontigs_1000\tproteins\tRF01734\tCrcB\tCLC\tCLC_F\tloci\n")
    for a in NAMES:
        s = table(os.path.join(OUT, a, "stats500.tsv"))[0]
        ex = table(os.path.join(OUT, a, "exporters.tsv"))
        n1k = sum(1 for _ in io.open(os.path.join(OUT, a, "contig_lengths.tsv")))
        nprot = sum(1 for l in io.open(os.path.join(OUT, a, "proteins.faa")) if l[0] == ">")
        nribo = sum(1 for l in io.open(os.path.join(OUT, a, "ribo.tbl")) if l.strip() and l[0] != "#")
        f.write("\t".join(map(str, [a, s["num_seqs"], s["sum_len"], s["N50"], s["max_len"], s["GC(%)"], n1k, nprot, nribo,
                                    sum(r["type"] == "CrcB" for r in ex), sum(r["type"] != "CrcB" for r in ex),
                                    sum(r["type"] == "CLC^F" for r in ex), sum(1 for x in loci.values() if x["asm"] == a)])) + "\n")
print("таблицы записаны в", RES)
