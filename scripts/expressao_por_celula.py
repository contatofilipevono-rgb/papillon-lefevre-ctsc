"""Checa se CTSC é expresso (previsão de RNA-seq do AlphaGenome 2) nas células candidatas.

Candidatos (PROMPT_ANTIGRAVITY.md): neutrófilo/mieloide, queratinócito, fibroblasto.
Controles: GAPDH (positivo) e ALB (fígado, negativo).
Uso: <alphagenome-venv-python> expressao_por_celula.py <.env> <saida.csv>
"""
import sys
import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

GENES = {  # GRCh38
    "CTSC": ("chr11", 88291570, 88337763),
    "GAPDH": ("chr12", 6534012, 6538930),
    "ALB": ("chr4", 73397114, 73426607),
}
CELULAS = {
    "EFO:0002067": "K562_mieloide",
    "EFO:0002793": "HL60_promielocito",
    "CL:0011022": "fibroblasto_pele",
    "UBERON:0036149": "pele_suprapubica",
    "UBERON:0013777": "pele_palma",
}


def chave(env):
    for ln in open(env, encoding="utf-8"):
        if ln.startswith("ALPHAGENOME_API_KEY="):
            return ln.split("=", 1)[1].strip()
    raise SystemExit("ALPHAGENOME_API_KEY não encontrada")


def main():
    env, saida = sys.argv[1:3]
    model = dna_client.create(chave(env))
    linhas = []
    for cond, nome in CELULAS.items():
        for g, (ch, ini, fim) in GENES.items():
            meio = (ini + fim) // 2
            iv = genome.Interval(ch, meio - dna_client.SEQUENCE_LENGTH_100KB // 2,
                                 meio + dna_client.SEQUENCE_LENGTH_100KB // 2)
            try:
                out = model.predict_interval(iv, requested_outputs=[dna_client.OutputType.RNA_SEQ], condition=cond)
                td = out.rna_seq
                vals = np.asarray(td.values)
                res = vals.shape[0] / iv.width
                a = max(0, int((ini - iv.start) * res))
                b = min(vals.shape[0], int((fim - iv.start) * res))
                media = float(vals[a:b].mean(axis=0).max()) if b > a else float("nan")
                pico = float(np.percentile(vals[a:b], 99, axis=0).max()) if b > a else float("nan")
                nt = vals.shape[1]
            except Exception as e:  # noqa: BLE001
                print(f"  Aviso {nome} ({cond}) {g}: {e}", flush=True)
                media, pico, nt = float("nan"), float("nan"), 0
            linhas.append({
                "celula": nome, "condicao": cond, "gene": g,
                "rna_media_no_gene": round(media, 4) if media == media else "",
                "rna_p99_no_gene": round(pico, 4) if pico == pico else "",
                "n_trilhas": nt,
            })
            print(nome, g, round(media, 4) if media == media else "NA", round(pico, 4) if pico == pico else "NA", flush=True)
    df = pd.DataFrame(linhas)
    df_num = df.copy()
    df_num["rna_p99_no_gene"] = pd.to_numeric(df_num["rna_p99_no_gene"], errors="coerce")
    gap = df_num[df_num.gene == "GAPDH"].set_index("celula")["rna_p99_no_gene"]
    alb = df_num[df_num.gene == "ALB"].set_index("celula")["rna_p99_no_gene"]
    df["p99_relativo_gapdh"] = (df_num["rna_p99_no_gene"] / df_num["celula"].map(gap)).round(4)
    df["p99_vezes_alb"] = (df_num["rna_p99_no_gene"] / df_num["celula"].map(alb)).round(1)
    df.to_csv(saida, index=False, encoding="utf-8-sig")
    print(df.pivot(index="gene", columns="celula", values="rna_p99_no_gene").to_string())


if __name__ == "__main__":
    main()
