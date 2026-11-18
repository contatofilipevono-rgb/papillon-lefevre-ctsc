"""Baixa do ClinVar as variantes NÃO codificantes (UTR, íntron, splice, regulatórias) de CTSC,
com coordenadas GRCh38 (SPDI), para pontuar no AlphaGenome e SpliceAI.

Uso: py -I -X utf8 clinvar_nao_codificantes.py <saida.csv>
"""
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter

UA = {"User-Agent": "academic-variant-curation/0.1"}
CROMO = {"NC_000011.10": "chr11"}  # CTSC em 11q14.2 (GRCh38)
EXCLUIR = ("missense", "nonsense", "frameshift", "synonymous", "inframe", "stop")


def get(url):
    for i in range(4):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120).read()
        except Exception:  # noqa: BLE001
            if i == 3:
                raise
            time.sleep(2 * (i + 1))


def main():
    saida = sys.argv[1]
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    linhas = []
    gene = "CTSC"
    q = urllib.parse.urlencode({"db": "clinvar", "term": f"{gene}[gene] AND single_gene[prop]",
                                "retmax": 5000, "retmode": "json"})
    ids = json.loads(get(base + "esearch.fcgi?" + q))["esearchresult"]["idlist"]
    for i in range(0, len(ids), 200):
        q = urllib.parse.urlencode({"db": "clinvar", "id": ",".join(ids[i:i + 200]), "retmode": "json"})
        res = json.loads(get(base + "esummary.fcgi?" + q))["result"]
        for uid in res.get("uids", []):
            d = res[uid]
            cons = ";".join(d.get("molecular_consequence_list", []) or [])
            if any(x in cons.lower() for x in EXCLUIR):
                continue
            vs = (d.get("variation_set") or [{}])[0]
            spdi = vs.get("canonical_spdi", "")
            partes = spdi.split(":")
            if len(partes) != 4 or partes[0] not in CROMO:
                continue
            acc, pos0, ref, alt = partes
            if len(ref) > 50 or len(alt) > 50:
                continue
            germ = d.get("germline_classification") or {}
            linhas.append({
                "gene": gene, "accession": d.get("accession", ""), "title": d.get("title", ""),
                "consequencia": cons or vs.get("variant_type", ""),
                "classificacao": germ.get("description", ""),
                "revisao": germ.get("review_status", ""),
                "cromossomo": CROMO[acc], "posicao_1based": int(pos0) + 1, "ref": ref, "alt": alt,
                "spdi": spdi,
            })
        time.sleep(0.4)
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()) if linhas else ["gene"])
        w.writeheader()
        w.writerows(linhas)
    print(len(linhas))
    for k, n in Counter((r["gene"], r["classificacao"]) for r in linhas).most_common():
        print(n, k)


if __name__ == "__main__":
    main()
