"""Pontua variantes não codificantes de CTSC com o AlphaGenome (API v2).

Uso: <alphagenome-venv-python> alphagenome_score.py <clinvar_nao_codificantes.csv> <.env> <saida.csv>
Condições selecionadas por expressão prevista de CTSC (expressao_por_celula.csv):
  - EFO:0002067 (K562, mieloide, principal; 162.5x ALB, observado)
  - EFO:0002793 (HL-60, promielocítica; 159.7x ALB)
  - CL:0011022  (fibroblasto de pele; 39.5x ALB, observado)
"""
import csv
import sys
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client
from alphagenome.models.v2 import variant_scorers as vs

TIPOS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "RNA_SEQ"]
CONDICOES = {
    "EFO:0002067": "k562_mieloide",
    "EFO:0002793": "hl60_promielocito",
    "CL:0011022": "fibroblasto_pele",
}


def chave(env):
    for ln in open(env, encoding="utf-8"):
        if ln.startswith("ALPHAGENOME_API_KEY="):
            return ln.split("=", 1)[1].strip()
    raise SystemExit("chave não encontrada")


def main():
    entrada, env, saida = sys.argv[1:4]
    with open(entrada, encoding="utf-8") as f:
        linhas = [r for r in csv.DictReader(f) if "p." not in r["title"]]
    for r in linhas:
        c = r["classificacao"]
        r["grupo"] = "P" if "athogenic" in c and "onflict" not in c else ("B" if "enign" in c else "U")
    linhas = [r for r in linhas if r["ref"] and r["alt"]]
    print("variantes não codificantes:", len(linhas), flush=True)
    model = dna_client.create(chave(env))
    scorers = [vs.RECOMMENDED_VARIANT_SCORERS[t] for t in TIPOS]
    variantes = [genome.Variant(chromosome=r["cromossomo"], position=int(r["posicao_1based"]),
                                reference_bases=r["ref"], alternate_bases=r["alt"]) for r in linhas]
    intervalos = [v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB) for v in variantes]
    resumo = [{k: r[k] for k in ("gene", "accession", "title", "consequencia", "classificacao", "grupo",
                                 "cromossomo", "posicao_1based", "ref", "alt")} for r in linhas]
    for cond, nome in CONDICOES.items():
        res = model.score_variants(intervalos, variantes, conditions=[cond] * len(variantes),
                                   variant_scorers=scorers, max_workers=2)
        brutos = []
        for r, linha, scores in zip(linhas, resumo, res):
            df = vs.tidy_scores(scores)
            df["accession"] = r["accession"]
            df["condicao"] = nome
            brutos.append(df)
            if "gene_name" in df:
                df = df[df["gene_name"].isna() | (df["gene_name"] == r["gene"])]
            col = "raw_score" if "raw_score" in df else "quantile_score"
            col_q = "quantile_score" if "quantile_score" in df else col
            for t in TIPOS:
                sub = df[df["output_type"].astype(str).str.upper().str.endswith(t)]
                q = float(sub[col_q].abs().max()) if len(sub) else float("nan")
                raw = float(sub[col].abs().max()) if len(sub) else float("nan")
                linha[f"{nome}_{t.lower()}_q"] = round(q, 4) if q == q else ""
                linha[f"{nome}_{t.lower()}_raw"] = round(raw, 4) if raw == raw else ""
        bruto = pd.concat(brutos)
        bruto.to_csv(saida.replace(".csv", f"_bruto_{nome}.csv"), index=False, encoding="utf-8-sig")
        print(f"  {nome}: {len(bruto)} linhas brutas", flush=True)
    for linha in resumo:
        for t in TIPOS:
            vals = [float(linha[f"{n}_{t.lower()}_q"]) for n in CONDICOES.values() if linha[f"{n}_{t.lower()}_q"] != ""]
            linha[f"max_{t.lower()}_q"] = max(vals) if vals else ""
            vals_raw = [float(linha[f"{n}_{t.lower()}_raw"]) for n in CONDICOES.values() if linha[f"{n}_{t.lower()}_raw"] != ""]
            linha[f"max_{t.lower()}_raw"] = max(vals_raw) if vals_raw else ""
    pd.DataFrame(resumo).to_csv(saida, index=False, encoding="utf-8-sig")
    print(f"Salvo {saida}", flush=True)


if __name__ == "__main__":
    main()
