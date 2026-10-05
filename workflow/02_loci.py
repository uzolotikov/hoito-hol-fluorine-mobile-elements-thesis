#!/usr/bin/env python3
"""Локусы и транспозазы в одной сборке.

1. Клада CLC^F: на дереве белков CLC ищется разрез, лучше всего отделяющий референсы CLC^F
   от референсов ClcA/ClcB; находки PF00654 по ту же сторону, что CLC^F, считаются CLC^F.
2. Локус: рибопереключатель RF01734 и ген экспортёра (CrcB или CLC^F любой клады) на одной цепи,
   ген ниже рибопереключателя, промежуток не более 200 п. н.
3. Транспозазы: гены с доменами из 146 профилей ферментов МГЭ в окне +-10 генов вокруг экспортёров.
   Ожидаемое число считается по доле таких генов во всей сборке отдельно для четырёх классов длины
   контига (< 5, 5-20, 20-100, > 100 тыс. п. н.); значимость по биномиальному распределению.
4. Для geNomad отбираются все контиги с находками и до 60 контрольных контигов длиной +-25 %
   на каждый контиг с локусом.

Запуск: workflow/02_loci.py <каталог сборки в OUT>
"""
import bisect, io, math, os, random, re, sys
from collections import defaultdict

D = sys.argv[1]
NAME = os.path.basename(os.path.normpath(D))
MAXGAP, FLANK, PAD = 200, 10, 300
LB = [0, 5000, 20000, 100000, 10**12]
LB_NAMES = ["<5k", "5-20k", "20-100k", ">100k"]


def path(f):
    return os.path.join(D, f)


clen = {}
for l in io.open(path("contig_lengths.tsv"), encoding="utf-8"):
    p = l.rstrip("\n").split("\t")
    clen[p[0]] = int(p[1])
coord, norf = defaultdict(dict), defaultdict(int)
for l in io.open(path("proteins.faa"), encoding="utf-8"):
    if l[0] != ">":
        continue
    m = re.match(r">(\S+)_(\d+) # (\d+) # (\d+) # (-?1) #", l)
    c, n = m.group(1), int(m.group(2))
    coord[c][n] = (int(m.group(3)), int(m.group(4)), "+" if m.group(5) == "1" else "-")
    norf[c] = max(norf[c], n)
TOTAL = sum(norf.values())
print("%s: белков %d на %d контигах (контигов >= 1000 п. н.: %d)" % (NAME, TOTAL, len(norf), len(clen)))


def tbl_ids(f):
    return list(dict.fromkeys(l.split()[0] for l in io.open(path(f), encoding="utf-8") if l.strip() and l[0] != "#"))


mge = set(tbl_ids("mge.tbl"))


# ---------- клада CLC^F ----------
class Node:
    __slots__ = ("name", "children", "parent")

    def __init__(self):
        self.name, self.children, self.parent = None, [], None


def parse_newick(t):
    i, root = 0, Node()
    cur = root
    while i < len(t):
        ch = t[i]
        if ch == "(":
            c = Node(); c.parent = cur; cur.children.append(c); cur = c; i += 1
        elif ch == ",":
            c = Node(); c.parent = cur.parent; cur.parent.children.append(c); cur = c; i += 1
        elif ch == ")":
            cur = cur.parent; i += 1
        elif ch == ";":
            break
        else:
            j = i
            while j < len(t) and t[j] not in "(),;":
                j += 1
            lab = t[i:j].split(":")[0]
            if lab and not re.match(r"^[\d.]+$", lab):
                cur.name = lab
            i = j
    return root


def leaves(n):
    if not n.children:
        return [n.name] if n.name else []
    return [x for c in n.children for x in leaves(c)]


def nodes(n):
    yield n
    for c in n.children:
        yield from nodes(c)


root = parse_newick(io.open(path("clc_all.nwk")).read().strip())
ALL = set(leaves(root))
RF = {x for x in ALL if x.startswith("REF_CLCF_")}
RA = {x for x in ALL if x.startswith(("REF_ClcA_", "REF_ClcB_"))}
best = None
for n in nodes(root):
    if not n.children:
        continue
    side = set(leaves(n))
    if not side or side == ALL:
        continue
    for s in (side, ALL - side):
        score = len(RF & s) + len(RA - s)
        if best is None or score > best[0]:
            best = (score, s)
clcf_side = best[1]
print("дерево CLC: референсы разделены %d из %d" % (best[0], len(RF) + len(RA)))
clade = {x: ("CLC^F" if x in clcf_side else "other_CLC") for x in ALL - RF - RA}

exporters = [(o, "CrcB") for o in tbl_ids("crcb.tbl")] + [(o, clade.get(o, "other_CLC")) for o in tbl_ids("clc.tbl")]
print("экспортёры: CrcB %d, CLC %d, из них CLC^F %d" % (
    sum(k == "CrcB" for _, k in exporters), sum(k != "CrcB" for _, k in exporters),
    sum(k == "CLC^F" for _, k in exporters)))


def split_orf(o):
    c, n = re.match(r"(.+)_(\d+)$", o).groups()
    return c, int(n)


by_contig = defaultdict(list)
for o, k in exporters:
    c, n = split_orf(o)
    s, e, st = coord[c][n]
    by_contig[c].append(dict(orf=o, c=c, n=n, gs=min(s, e), ge=max(s, e), strand=st, kind=k))
with io.open(path("exporters.tsv"), "w", encoding="utf-8") as f:
    f.write("orf\tcontig\tstart\tend\tstrand\ttype\n")
    for gs in by_contig.values():
        for g in gs:
            f.write("%s\t%s\t%d\t%d\t%s\t%s\n" % (g["orf"], g["c"], g["gs"], g["ge"], g["strand"], g["kind"]))

ribo = []
for l in io.open(path("ribo.tbl"), encoding="utf-8"):
    if l[0] == "#" or not l.strip():
        continue
    p = l.split()
    a, b = int(p[7]), int(p[8])
    ribo.append(dict(c=p[0], lo=min(a, b), hi=max(a, b), strand=p[9], score=float(p[14])))
print("рибопереключателей RF01734: %d" % len(ribo))


# ---------- локусы ----------
def gap(r, g):
    if r["strand"] != g["strand"]:
        return None
    if g["strand"] == "+":
        return g["gs"] - r["hi"] - 1 if r["hi"] < g["gs"] else None
    return r["lo"] - g["ge"] - 1 if r["lo"] > g["ge"] else None


loci = []
for r in ribo:
    best_g = None
    for g in by_contig.get(r["c"], []):
        gp = gap(r, g)
        if gp is not None and 0 <= gp <= MAXGAP and (best_g is None or gp < best_g[1]):
            best_g = (g, gp)
    if best_g:
        loci.append((r, best_g[0], best_g[1]))


# ---------- транспозазы ----------
def window(c, n, f=FLANK):
    return ["%s_%d" % (c, i) for i in range(max(1, n - f), min(norf.get(c, 0), n + f) + 1) if i != n]


def lbin(c):
    return bisect.bisect_right(LB, clen.get(c, 0)) - 1


bn, bk = defaultdict(int), defaultdict(int)
for c, n in norf.items():
    bn[lbin(c)] += n
for o in mge:
    bk[lbin(split_orf(o)[0])] += 1
rate = {b: (bk[b] / bn[b] if bn[b] else 0.0) for b in range(4)}
print("доля генов ферментов МГЭ: %.3f %%; по классам длины %s" % (
    100 * len(mge) / TOTAL, ", ".join("%s %.3f %%" % (LB_NAMES[b], 100 * rate[b]) for b in range(4))))


def binom_tail(k, n, p, upper=True):
    if n == 0 or p <= 0:
        return 1.0
    rng = range(k, n + 1) if upper else range(0, k + 1)
    return min(1.0, sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1)
                                 + i * math.log(p) + (n - i) * math.log1p(-p)) for i in rng))


allhits = [g for gs in by_contig.values() for g in gs]
sets = [("CrcB+CLC^F", [g for g in allhits if g["kind"] in ("CrcB", "CLC^F")]),
        ("all_exporters", allhits),
        ("other_CLC", [g for g in allhits if g["kind"] == "other_CLC"]),
        ("loci", [g for _, g, _ in loci])]
with io.open(path("transposase_test.tsv"), "w", encoding="utf-8") as f:
    f.write("assembly\tset\tgenes\twindow_orfs\tmge_observed\tmge_expected\tP_ge\tP_le\n")
    for name, sub in sets:
        orfs = list(dict.fromkeys(x for g in sub for x in window(g["c"], g["n"])))
        k, n = sum(1 for x in orfs if x in mge), len(orfs)
        e = sum(rate[lbin(split_orf(x)[0])] for x in orfs)
        p = e / n if n else 0
        row = (NAME, name, len(sub), n, k, e, binom_tail(k, n, p), binom_tail(k, n, p, False))
        f.write("%s\t%s\t%d\t%d\t%d\t%.2f\t%.3g\t%.3g\n" % row)
        print("транспозазы, %-13s генов %3d, ORF в окнах %4d, найдено %d, ожидалось %.1f, P(>=) %.3g, P(<=) %.3g" % row[1:])

with io.open(path("loci.tsv"), "w", encoding="utf-8") as f:
    f.write("contig\tcontig_len\tribo_start\tribo_end\tstrand\tribo_bitscore\texporter\torf\torf_start\torf_end"
            "\tgap\twindow_orfs\tedge_dist\tmge_in_window\n")
    for r, g, gp in sorted(loci, key=lambda x: -x[0]["score"]):
        L = clen[g["c"]]
        lo, hi = min(r["lo"], g["gs"]), max(r["hi"], g["ge"])
        w = window(g["c"], g["n"])
        f.write("\t".join(map(str, [g["c"], L, r["lo"], r["hi"], r["strand"], r["score"], g["kind"], g["orf"],
                                    g["gs"], g["ge"], gp, len(w), min(lo, L - hi), sum(x in mge for x in w)])) + "\n")
print("локусов: %d на %d контигах" % (len(loci), len({g["c"] for _, g, _ in loci})))

# ---------- набор для geNomad ----------
fluo = sorted(set(by_contig) | {r["c"] for r in ribo})
locus_contigs = sorted({g["c"] for _, g, _ in loci})
pool = sorted((c for c in clen if c not in set(fluo)), key=lambda c: clen[c])
pool_len = [clen[c] for c in pool]
random.seed(1)
controls = {}
for c in locus_contigs:
    L, fr = clen[c], 0.25
    while True:
        a = bisect.bisect_left(pool_len, int(L * (1 - fr)))
        b = bisect.bisect_right(pool_len, int(L * (1 + fr)))
        if b - a >= 20 or fr > 3:
            break
        fr *= 1.5
    cand = pool[a:b]
    controls[c] = random.sample(cand, min(60, len(cand)))
with io.open(path("controls.tsv"), "w", encoding="utf-8") as f:
    for c in locus_contigs:
        f.write(c + "\t" + ",".join(controls[c]) + "\n")
ids = sorted(set(fluo) | {x for v in controls.values() for x in v})
io.open(path("genomad_ids.txt"), "w").write("\n".join(ids) + "\n")
io.open(path("locus_contigs.txt"), "w").write("\n".join(locus_contigs) + "\n")
print("для geNomad: контигов с находками %d, контрольных %d" % (len(fluo), len(ids) - len(fluo)))


# ---------- участки локусов для сверки сборок ----------
def read_fasta(fn, want):
    seqs, name, buf = {}, None, []
    for l in io.open(fn, encoding="utf-8"):
        if l[0] == ">":
            if name in want:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
        elif name in want:
            buf.append(l.strip())
    if name in want:
        seqs[name] = "".join(buf)
    return seqs


seq = read_fasta(path("contigs1k.fna"), set(locus_contigs))
with io.open(path("regions.fna"), "w") as f:
    for r, g, gp in loci:
        s = seq[g["c"]]
        lo0, hi0 = min(r["lo"], g["gs"]), max(r["hi"], g["ge"])
        f.write(">%s@%s@%s@%d-%d\n%s\n" % (NAME, g["c"], g["orf"], lo0, hi0,
                                           s[max(0, lo0 - 1 - PAD):min(len(s), hi0 + PAD)]))
