"""Busca de ineditismo no Europe PMC (título + resumo) para análise estrutural de variantes do CTSC."""
import json
import time
import urllib.parse
import urllib.request

GENE = '(TITLE:"cathepsin C" OR ABSTRACT:"cathepsin C" OR TITLE:CTSC OR ABSTRACT:CTSC OR ABSTRACT:"dipeptidyl peptidase I")'
Q = {
    "estrutura_mutacoes": GENE + ' AND (ABSTRACT:mutation* OR ABSTRACT:variant*) AND (ABSTRACT:structur* OR ABSTRACT:modelling OR ABSTRACT:modeling OR ABSTRACT:"molecular dynamics" OR ABSTRACT:stability)',
    "in_silico": GENE + ' AND (ABSTRACT:"in silico" OR ABSTRACT:computational OR ABSTRACT:nsSNP* OR ABSTRACT:bioinformatic*)',
    "ferramentas_ia": GENE + ' AND (ABSTRACT:AlphaFold OR ABSTRACT:AlphaMissense OR ABSTRACT:FoldX OR ABSTRACT:"free energy" OR ABSTRACT:DDG OR ABSTRACT:ThermoMPNN OR ABSTRACT:AlphaGenome OR ABSTRACT:SpliceAI)',
    "espectro_mutacoes": GENE + ' AND (ABSTRACT:"Papillon-Lefevre" OR ABSTRACT:"Papillon-Lefèvre") AND (ABSTRACT:review OR ABSTRACT:"genotype-phenotype" OR ABSTRACT:spectrum OR ABSTRACT:"systematic")',
}


def busca(q):
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "pageSize": 60, "resultType": "lite", "sort": "P_PDATE_D desc"})
    for tentativa in range(5):
        try:
            return json.load(urllib.request.urlopen(url, timeout=60))
        except Exception as e:
            print("  (erro", e, "— nova tentativa)")
            time.sleep(5 * (tentativa + 1))
    raise SystemExit("Europe PMC indisponível")


for nome, q in Q.items():
    r = busca(q)
    print(f"\n### {nome}: {r['hitCount']} resultados")
    for a in r["resultList"]["result"]:
        print(f"- {a.get('pubYear')} | {a.get('title', '')[:160]} | {a.get('journalTitle', '')} | pmid:{a.get('pmid', '')} | pmcid:{a.get('pmcid', '')}")
    time.sleep(2)
