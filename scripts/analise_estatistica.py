"""Análise estatística pré-registrada (PRE_REGISTRO.md) para CTSC / Papillon-Lefèvre.

Executa:
  1. Comparação contínua (Mann-Whitney exato/assintótico via scipy.stats + AUC + IC bootstrap 95%).
  2. Enriquecimento por mecanismo estrutural (Fisher exato bicaudal).
  3. Validação cruzada leave-one-out (LOO) aninhada no conjunto de calibração (P + P-lit<=2015 vs. controle AC>=2).
  4. Avaliação no conjunto de teste externo independente (P-lit>2015, n=31) e no controle homozigoto (hom>=1).
  5. Análise de sensibilidade aos cortes e decomposição das VUS em subgrupos disjuntos.

Uso: py -I -X utf8 analise_estatistica.py <pasta_papillon_lefevre>
"""
import csv
import math
import random
import sys
from pathlib import Path
from statistics import median
from scipy import stats


def wilson_ci(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def auc_val(pos, neg):
    u, _ = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (len(pos) * len(neg))


def bootstrap_auc(pos, neg, n_boot=2000, seed=42):
    rng = random.Random(seed)
    vals = []
    for _ in range(n_boot):
        p_s = [rng.choice(pos) for _ in range(len(pos))]
        n_s = [rng.choice(neg) for _ in range(len(neg))]
        u, _ = stats.mannwhitneyu(p_s, n_s, alternative="two-sided")
        vals.append(u / (len(p_s) * len(n_s)))
    vals.sort()
    return vals[int(0.025 * n_boot)], vals[int(0.975 * n_boot)]


def bootstrap_lr(tp, fn, tn, fp, n_boot=2000, seed=42):
    rng = random.Random(seed)
    p_list = [1] * tp + [0] * fn
    n_list = [1] * fp + [0] * tn
    lr_pos, lr_neg = [], []
    for _ in range(n_boot):
        ps = [rng.choice(p_list) for _ in range(len(p_list))]
        ns = [rng.choice(n_list) for _ in range(len(n_list))]
        sens = sum(ps) / len(ps)
        fpr = sum(ns) / len(ns)
        spec = 1 - fpr
        if fpr > 0:
            lr_pos.append(sens / fpr)
        if spec > 0:
            lr_neg.append((1 - sens) / spec)
    lr_pos.sort()
    lr_neg.sort()
    ci_p = (lr_pos[int(0.025 * len(lr_pos))], lr_pos[int(0.975 * len(lr_pos))]) if lr_pos else (0, 0)
    ci_n = (lr_neg[int(0.025 * len(lr_neg))], lr_neg[int(0.975 * len(lr_neg))]) if lr_neg else (0, 0)
    return ci_p, ci_n


def regra_primaria(r, d_act=4.0, d_tet=4.0, cg=1.5):
    da = float(r["dist_sitio_ativo"])
    dt = float(r["dist_interface_tetramero"])
    ddg = float(r["ddg_1k3b"])
    return da < d_act or dt < d_tet or r["dissulfeto"] == "sim" or ddg >= cg


def regra_minimalista(r, d_act=4.0, cg=1.5):
    da = float(r["dist_sitio_ativo"])
    ddg = float(r["ddg_1k3b"])
    return da < d_act or ddg >= cg


def main():
    base = Path(sys.argv[1])
    rows = [r for r in csv.DictReader(open(base / "variantes_CTSC.csv", encoding="utf-8-sig")) if r["na_1k3b"] == "sim"]

    P_cal = [r for r in rows if r["coorte_patogenica"] == "calibracao"]
    P_ext = [r for r in rows if r["coorte_patogenica"] == "teste_externo"]
    P_all = P_cal + P_ext
    C_hom = [r for r in rows if r["rotulo"] not in ("P", "P-lit") and int(r["gnomad_hom"] or 0) >= 1]
    C_ac2 = [r for r in rows if r["rotulo"] not in ("P", "P-lit") and int(r["gnomad_ac"] or 0) >= 2]

    print(f"=== CONJUNTOS NA ESTRUTURA 1K3B ===")
    print(f"Patogênicas calibração (P + P-lit<=2015): {len(P_cal)}")
    print(f"Patogênicas teste externo (P-lit>2015)  : {len(P_ext)}")
    print(f"Patogênicas totais                      : {len(P_all)}")
    print(f"Controle primário (gnomAD hom>=1)       : {len(C_hom)}")
    print(f"Controle sensibilidade (gnomAD AC>=2)   : {len(C_ac2)}")

    # 1. Variáveis contínuas (P_all vs C_ac2, e P_all vs C_hom)
    print("\n=== 1. COMPARAÇÕES CONTÍNUAS (Patogênicas totais n=51 vs. Controle AC>=2 n=234) ===")
    for nome, col in (("ΔΔG ThermoMPNN (1K3B)", "ddg_1k3b"),
                      ("ΔΔG ThermoMPNN (AF-P53634)", "ddg_af"),
                      ("AlphaMissense", "alphamissense"),
                      ("RSA no tetrâmero", "rsa_tetramero")):
        a = [float(r[col]) for r in P_all if r[col] != ""]
        b = [float(r[col]) for r in C_ac2 if r[col] != ""]
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        auc = u / (len(a) * len(b))
        ci_lo, ci_hi = bootstrap_auc(a, b)
        print(f"  {nome:28s}: med P={median(a):.3f} vs C={median(b):.3f} | AUC={auc:.3f} (IC95% {ci_lo:.3f}-{ci_hi:.3f}) | p={p:.3e}")

    # 2. Mecanismos individuais (Fisher exato)
    print("\n=== 2. MECANISMOS ESTRUTURAIS (Patogênicas totais n=51 vs. Controle AC>=2 n=234) ===")
    mecs = [
        ("Sítio ativo < 4,0 Å", lambda r: float(r["dist_sitio_ativo"]) < 4.0),
        ("Interface tetrâmero < 4,0 Å", lambda r: float(r["dist_interface_tetramero"]) < 4.0),
        ("Íon cloreto < 5,0 Å", lambda r: float(r["dist_cloreto"]) < 5.0),
        ("Ponte dissulfeto", lambda r: r["dissulfeto"] == "sim"),
        ("N-glicosilação", lambda r: r["glicosilacao"] == "sim"),
        ("ΔΔG (1K3B) >= 1,5 kcal/mol", lambda r: float(r["ddg_1k3b"]) >= 1.5),
    ]
    for nome, fn in mecs:
        tp = sum(fn(r) for r in P_all)
        fp = sum(fn(r) for r in C_ac2)
        _, p = stats.fisher_exact([[tp, len(P_all) - tp], [fp, len(C_ac2) - fp]], alternative="two-sided")
        print(f"  {nome:30s}: P={tp}/{len(P_all)} ({tp/len(P_all):.1%}) vs C={fp}/{len(C_ac2)} ({fp/len(C_ac2):.1%}) | Fisher p={p:.3e}")

    # 3. LOO aninhada no conjunto de calibração (P_cal vs C_ac2)
    print("\n=== 3. VALIDAÇÃO CRUZADA LOO ANINHADA (Calibração n=20 vs. Controle AC>=2 n=234) ===")
    cortes = [i * 0.25 for i in range(13)]
    for nome_r, fn_r in (("Regra Primária (4 termos)", regra_primaria), ("Regra Minimalista (2 termos)", regra_minimalista)):
        todos = [(r, 1) for r in P_cal] + [(r, 0) for r in C_ac2]
        tp = fn = tn = fp = 0
        escolhidos = []
        for i, (r, y) in enumerate(todos):
            treino = [t for j, t in enumerate(todos) if j != i]
            def youden(c):
                pp = [fn_r(t, cg=c) for t, yy in treino if yy == 1]
                cc = [not fn_r(t, cg=c) for t, yy in treino if yy == 0]
                return sum(pp) / len(pp) + sum(cc) / len(cc) - 1
            c_opt = max(cortes, key=youden)
            escolhidos.append(c_opt)
            pred = fn_r(r, cg=c_opt)
            if y == 1:
                tp += pred; fn += not pred
            else:
                tn += not pred; fp += pred
        sens = tp / (tp + fn)
        esp = tn / (tn + fp)
        s_lo, s_hi = wilson_ci(tp, tp + fn)
        e_lo, e_hi = wilson_ci(tn, tn + fp)
        lr_p_ci, lr_n_ci = bootstrap_lr(tp, fn, tn, fp)
        moda = max(set(escolhidos), key=escolhidos.count)
        print(f"  {nome_r}:")
        print(f"    LOO cortes escolhidos: {min(escolhidos):.2f}-{max(escolhidos):.2f} (moda={moda:.2f})")
        print(f"    Sensibilidade LOO: {tp}/{tp+fn} = {sens:.3f} (Wilson IC95% {s_lo:.3f}-{s_hi:.3f})")
        print(f"    Especificidade LOO: {tn}/{tn+fp} = {esp:.3f} (Wilson IC95% {e_lo:.3f}-{e_hi:.3f})")
        print(f"    LR+={sens/(1-esp):.2f} (IC95% {lr_p_ci[0]:.2f}-{lr_p_ci[1]:.2f}), LR-={(1-sens)/esp:.2f} (IC95% {lr_n_ci[0]:.2f}-{lr_n_ci[1]:.2f})")

        # 4. Teste externo independente (P_ext n=31) com a moda da calibração
        tp_ext = sum(fn_r(r, cg=moda) for r in P_ext)
        ext_lo, ext_hi = wilson_ci(tp_ext, len(P_ext))
        tn_hom = sum(not fn_r(r, cg=moda) for r in C_hom)
        hom_lo, hom_hi = wilson_ci(tn_hom, len(C_hom))
        print(f"    TESTE EXTERNO (P-lit>2015, corte={moda:.2f}): {tp_ext}/{len(P_ext)} = {tp_ext/len(P_ext):.3f} (Wilson IC95% {ext_lo:.3f}-{ext_hi:.3f})")
        print(f"    CONTROLE HOMOZIGOTO (hom>=1, corte={moda:.2f}): {tn_hom}/{len(C_hom)} = {tn_hom/len(C_hom):.3f} (Wilson IC95% {hom_lo:.3f}-{hom_hi:.3f})")

    # 5. Sensibilidade aos cortes pré-registrados (todos P vs C_ac2)
    print("\n=== 4. SENSIBILIDADE AOS CORTES (Regra Primária: P_all n=51 vs C_ac2 n=234) ===")
    for cg in (0.5, 1.0, 1.5, 2.0):
        tp = sum(regra_primaria(r, cg=cg) for r in P_all)
        tn = sum(not regra_primaria(r, cg=cg) for r in C_ac2)
        _, p = stats.fisher_exact([[tp, len(P_all) - tp], [len(C_ac2) - tn, tn]])
        print(f"  ΔΔG >= {cg:.1f}: Sens={tp}/{len(P_all)} ({tp/len(P_all):.3f}), Espec={tn}/{len(C_ac2)} ({tn/len(C_ac2):.3f}), Fisher p={p:.3e}")

    # 6. Decomposição das VUS em subgrupos disjuntos
    print("\n=== 5. DECOMPOSIÇÃO DAS VUS (U) EM SUBGRUPOS DISJUNTOS (Regra Primária, c=1,5) ===")
    vus = [r for r in rows if r["rotulo"] == "U"]
    subgrupos = [
        ("Ausente do gnomAD (AC=0)", [r for r in vus if int(r["gnomad_ac"] or 0) == 0]),
        ("Rara no gnomAD (AC=1, fora calibração)", [r for r in vus if int(r["gnomad_ac"] or 0) == 1]),
        ("Recorrente no gnomAD (AC>=2, sobreposta ao controle)", [r for r in vus if int(r["gnomad_ac"] or 0) >= 2]),
    ]
    for nome_s, grp in subgrupos:
        pos = sum(regra_primaria(r, cg=1.5) for r in grp)
        lo, hi = wilson_ci(pos, len(grp))
        print(f"  {nome_s:52s}: {pos}/{len(grp)} ({pos/len(grp):.1%}, IC95% {lo:.1%}-{hi:.1%})")


if __name__ == "__main__":
    main()
