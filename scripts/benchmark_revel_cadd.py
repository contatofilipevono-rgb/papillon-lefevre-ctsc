"""Busca escores CADD e REVEL via myvariant.info (dbNSFP) para variantes missense do CTSC com coordenadas GRCh38 no gnomAD,
e compara AUC contra ThermoMPNN ΔΔG e AlphaMissense.

Uso: py -I -X utf8 benchmark_revel_cadd.py <pasta_papillon_lefevre>
"""
import csv
import json
import sys
import urllib.request
from pathlib import Path
from scipy import stats


def main():
    base = Path(sys.argv[1])
    aa1to3 = {
        "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
        "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
        "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
        "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    }
    gn_map = {}
    for r in csv.DictReader(open(base / "dados_brutos" / "gnomad_CTSC.csv", encoding="utf-8")):
        if r["hgvsp"] and r["variant_id"].count("-") == 3:
            ch, pos, ref, alt = r["variant_id"].split("-")
            if len(ref) == 1 and len(alt) == 1:
                gn_map[r["hgvsp"]] = f"chr{ch}:g.{pos}{ref}>{alt}"

    var_rows = [r for r in csv.DictReader(open(base / "variantes_CTSC.csv", encoding="utf-8-sig")) if r["na_1k3b"] == "sim"]
    P_all = [r for r in var_rows if r["coorte_patogenica"] in ("calibracao", "teste_externo")]
    C_ac2 = [r for r in var_rows if r["rotulo"] not in ("P", "P-lit") and int(r["gnomad_ac"] or 0) >= 2]

    ids_to_query = {}
    for r in P_all + C_ac2:
        v = r["variante"]
        p3 = f"p.{aa1to3[v[0]]}{v[1:-1]}{aa1to3[v[-1]]}"
        hgvs = gn_map.get(p3)
        if hgvs:
            ids_to_query[hgvs] = r

    print(f"Consultando myvariant.info para {len(ids_to_query)} variantes (P={sum(1 for r in ids_to_query.values() if r in P_all)}, C={sum(1 for r in ids_to_query.values() if r in C_ac2)})...")
    keys = list(ids_to_query.keys())
    scores = {}
    for i in range(0, len(keys), 100):
        chunk = keys[i:i + 100]
        body = json.dumps({"ids": chunk}).encode()
        req = urllib.request.Request(
            "https://myvariant.info/v1/variant?assembly=hg38&fields=dbnsfp.revel.score,dbnsfp.cadd.phred",
            data=body,
            headers={"Content-Type": "application/json", "User-Agent": "academic/0.1"},
        )
        res = json.loads(urllib.request.urlopen(req, timeout=30).read())
        for item in res:
            q = item.get("query")
            db = item.get("dbnsfp") or {}
            rev = (db.get("revel") or {}).get("score")
            cad = (db.get("cadd") or {}).get("phred")
            if isinstance(rev, list):
                rev = rev[0]
            if isinstance(cad, list):
                cad = cad[0]
            scores[q] = (rev, cad)

    for label, idx in (("REVEL", 0), ("CADD phred", 1)):
        p_vals = [float(scores[k][idx]) for k, r in ids_to_query.items() if r in P_all and k in scores and scores[k][idx] is not None]
        c_vals = [float(scores[k][idx]) for k, r in ids_to_query.items() if r in C_ac2 and k in scores and scores[k][idx] is not None]
        if p_vals and c_vals:
            u, p = stats.mannwhitneyu(p_vals, c_vals, alternative="two-sided")
            auc = u / (len(p_vals) * len(c_vals))
            print(f"  {label:12s}: AUC={auc:.3f} (n={len(p_vals)} P vs {len(c_vals)} C), p={p:.3e}")


if __name__ == "__main__":
    main()
