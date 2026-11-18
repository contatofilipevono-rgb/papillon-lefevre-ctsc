"""Busca e extrai variantes missense de CTSC e dados funcionais da literatura (Europe PMC / PubMed).

Fontes principais:
  - Nagy et al. 2014 (PMC4049362): revisão e banco de variantes PLS/HMS/AP
  - Chitsamankhun et al. 2024 (PMC11337848): revisão recente de CTSC
  - Relatos de caso e estudos funcionais no Europe PMC (ex.: Pham 2004, Hamon et al., estudos de urina 2022)

Uso: py -I -X utf8 extrair_literatura.py <pasta_dados_brutos>
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

UA = {"User-Agent": "academic-variant-curation/0.1"}
PMCIDS = [
    "PMC4049362",   # Nagy 2014
    "PMC11337848",  # Chitsamankhun 2024
    "PMC8404484",   # 2021 Iranian PLS
    "PMC8372118",   # 2021 Chinese PLS
    "PMC6801622",   # 2019 Processing & Maturation of CTSC Zymogen
    "PMC13447139",  # 2026 Front Immunol PLS
    "PMC8273319",   # 2021 Ann Dermatol Chinese PLS
    "PMC4992098",   # 2016 Saudi PLS
    "PMC4617240",   # 2010 Brazilian PLS
    "PMC2921986",   # 2010 Pakistani PLS
    "PMC1734492",   # Hart 2000 J Med Genet
]


def get(url, timeout=25, tries=4):
    req = urllib.request.Request(url, headers=UA)
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def main():
    out = Path(sys.argv[1]) / "literatura"
    out.mkdir(parents=True, exist_ok=True)

    # 1. Busca ampla no Europe PMC por artigos sobre CTSC / Papillon-Lefevre com variantes
    q = urllib.parse.urlencode({
        "query": '("CTSC" OR "cathepsin C") AND ("Papillon-Lefevre" OR "Haim-Munk" OR "aggressive periodontitis") AND ("missense" OR "mutation" OR "variant")',
        "format": "json",
        "resultType": "lite",
        "pageSize": 100,
    })
    results = []
    try:
        res = json.loads(get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{q}", timeout=30))
        results = res.get("resultList", {}).get("result", [])
        (out / "europepmc_ctsc_search.json").write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"Europe PMC search: {len(results)} artigos encontrados (hitCount={res.get('hitCount')})")
    except Exception as e:  # noqa: BLE001
        print(f"Aviso: busca Europe PMC falhou ({e}); usando lista de PMCIDs.")

    # Coleta PMCIDs open access relevantes publicados após 2014
    extra_pmcids = list(PMCIDS)
    for r in results:
        pmcid = r.get("pmcid")
        yr = int(r.get("pubYear") or 0)
        title = (r.get("title") or "").lower()
        if pmcid and r.get("isOpenAccess") == "Y" and yr >= 2014 and pmcid not in extra_pmcids:
            if any(k in title for k in ("papillon", "haim-munk", "ctsc", "cathepsin c")):
                extra_pmcids.append(pmcid)

    print(f"Baixando XML de {len(extra_pmcids)} artigos Open Access...", flush=True)
    ok = 0
    for pmcid in extra_pmcids:
        dest = out / f"{pmcid}.xml"
        if dest.exists() and dest.stat().st_size > 1000:
            ok += 1
            continue
        try:
            xml_bytes = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML", timeout=15, tries=2)
            dest.write_bytes(xml_bytes)
            ok += 1
            print(f"  OK {pmcid} ({len(xml_bytes)} bytes)", flush=True)
            time.sleep(0.2)
        except Exception as e:  # noqa: BLE001
            print(f"  Falhou {pmcid}: {e}", flush=True)
    print(f"XMLs salvos: {ok}/{len(extra_pmcids)}", flush=True)


if __name__ == "__main__":
    main()
