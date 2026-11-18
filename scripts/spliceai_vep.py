"""Busca escores SpliceAI (pré-calculados) pelo Ensembl VEP REST e compara com o AlphaGenome v2.

Uso: py -I -X utf8 spliceai_vep.py <alphagenome_nao_codificantes.csv> <saida.csv>
"""
import csv
import json
import sys
import time
import urllib.request
from scipy import stats


URL = "https://rest.ensembl.org/vep/human/region?SpliceAI=1"


def vep(lines):
    out_all = []
    for i in range(0, len(lines), 100):
        chunk = lines[i:i + 100]
        body = json.dumps({"variants": chunk}).encode()
        req = urllib.request.Request(URL, data=body, headers={
            "Content-Type": "application/json", "Accept": "application/json",
            "User-Agent": "academic-variant-curation/0.1"})
        for t in range(4):
            try:
                res = json.loads(urllib.request.urlopen(req, timeout=180).read())
                out_all.extend(res)
                break
            except Exception:  # noqa: BLE001
                if t == 3:
                    raise
                time.sleep(5 * (t + 1))
    return out_all


def auc(pos, neg):
    u, _ = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    return u / (len(pos) * len(neg))


def main():
    entrada, saida = sys.argv[1:3]
    with open(entrada, encoding="utf-8-sig") as f:
        rows = [r for r in csv.DictReader(f) if r.get("max_splice_sites_q", "") not in ("", "nan")]
    lines, ids = [], {}
    for r in rows:
        chrom = r["cromossomo"].replace("chr", "")
        vid = f"v{len(lines)}"
        lines.append(f"{chrom} {r['posicao_1based']} {vid} {r['ref']} {r['alt']} . . .")
        ids[vid] = r
    out = vep(lines)
    for res in out:
        r = ids.get(res.get("id"))
        if not r:
            continue
        ds, cons = [], set()
        for tc in res.get("transcript_consequences", []):
            if tc.get("gene_symbol") != r["gene"]:
                continue
            cons.update(tc.get("consequence_terms", []))
            s = tc.get("spliceai")
            if s:
                ds.append(max(s.get(k, 0) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")))
        r["spliceai_max"] = max(ds) if ds else ""
        r["vep_consequencias"] = ";".join(sorted(cons))
    campos = list(rows[0].keys()) + [c for c in ("spliceai_max", "vep_consequencias") if c not in rows[0]]
    with open(saida, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    ok = [r for r in rows if r.get("spliceai_max") not in ("", None) and "start_lost" not in r.get("vep_consequencias", "")]
    for nome, key in (("AlphaGenome v2 splicing (K562)", "k562_mieloide_splice_sites_q"),
                      ("AlphaGenome v2 splicing (max)", "max_splice_sites_q"),
                      ("SpliceAI", "spliceai_max")):
        P = [float(r[key]) for r in ok if r["grupo"] == "P" and r.get(key, "") != ""]
        B = [float(r[key]) for r in ok if r["grupo"] == "B" and r.get(key, "") != ""]
        if P and B:
            print(f"{nome:32s}: AUC={auc(P, B):.3f} ({len(P)} P, {len(B)} B); max B={max(B):.3f}; min P={min(P):.3f}")
    print("\nIncertas (U) ordenadas pelo AlphaGenome v2 (max_splice_sites_q):")
    for r in sorted([r for r in ok if r["grupo"] == "U"], key=lambda r: -float(r["max_splice_sites_q"]))[:10]:
        print(f"  {r['title'][:52]:52s} AG={float(r['max_splice_sites_q']):.3f} SpliceAI={float(r['spliceai_max']):.2f}")


if __name__ == "__main__":
    main()
