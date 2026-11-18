"""Analisa distância ao resíduo Asp25 (bolso S2 do domínio de exclusão) e à interface
entre o domínio de exclusão (cadeia A) e o domínio catalítico (cadeias B+C) na 1K3B.

Uso: py -I -X utf8 analise_subsitio_exclusao.py <pasta_papillon_lefevre>
"""
import csv
import math
import sys
from pathlib import Path
from scipy import stats


def parse_1k3b(cif_path):
    cols, res_coords = [], {}
    for ln in cif_path.read_text(encoding="utf-8").splitlines():
        if ln.startswith("_atom_site."):
            cols.append(ln.split(".", 1)[1].strip())
        elif ln.startswith("ATOM") and cols:
            r = dict(zip(cols, ln.split()))
            if r.get("type_symbol") == "H":
                continue
            ch = r["auth_asym_id"]
            if ch not in ("A", "B", "C"):
                continue
            up_pos = int(r["auth_seq_id"]) + 24
            xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            res_coords.setdefault((ch, up_pos), []).append(xyz)
    return res_coords


def min_dist(c1, c2):
    return min(math.dist(a, b) for a in c1 for b in c2)


def main():
    base = Path(sys.argv[1])
    rc = parse_1k3b(base / "dados_brutos" / "1K3B.cif")
    asp25_coords = rc[("A", 25)]
    excl_coords = [xyz for (ch, _), pts in rc.items() if ch == "A" for xyz in pts]
    cat_coords = [xyz for (ch, _), pts in rc.items() if ch in ("B", "C") for xyz in pts]

    pos_coords = {pos: pts for (_, pos), pts in rc.items()}
    d_asp25 = {p: min_dist(c, asp25_coords) for p, c in pos_coords.items()}
    d_excl_cat = {
        p: (min_dist(c, cat_coords) if p <= 143 else min_dist(c, excl_coords))
        for p, c in pos_coords.items()
    }

    rows = [r for r in csv.DictReader(open(base / "variantes_CTSC.csv", encoding="utf-8-sig")) if r["na_1k3b"] == "sim"]
    P = [r for r in rows if r["coorte_patogenica"] in ("calibracao", "teste_externo")]
    C = [r for r in rows if r["rotulo"] not in ("P", "P-lit") and int(r["gnomad_ac"] or 0) >= 2]

    for nome, dmap, corte in (
        ("Distância a Asp25 (bolso S2) < 6,0 Å", d_asp25, 6.0),
        ("Interface exclusão-catalítica < 4,0 Å", d_excl_cat, 4.0),
    ):
        tp = sum(dmap[int(r["posicao"])] < corte for r in P)
        fp = sum(dmap[int(r["posicao"])] < corte for r in C)
        _, p_f = stats.fisher_exact([[tp, len(P) - tp], [fp, len(C) - fp]])
        vp = [dmap[int(r["posicao"])] for r in P]
        vc = [dmap[int(r["posicao"])] for r in C]
        _, p_mw = stats.mannwhitneyu(vp, vc, alternative="two-sided")
        print(f"{nome:42s}: P={tp}/{len(P)} ({tp/len(P):.1%}) vs C={fp}/{len(C)} ({fp/len(C):.1%}) | Fisher p={p_f:.3e} | MW p={p_mw:.3e}")


if __name__ == "__main__":
    main()
