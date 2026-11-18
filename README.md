# Structural Mechanisms of Pathogenic *CTSC* Variants in Papillon-Lefèvre Syndrome

Computational structural biology and regulatory variant analysis of human cathepsin C (*CTSC*, UniProt `P53634`, PDB `1K3B`) in Papillon-Lefèvre syndrome, Haim-Munk syndrome, and aggressive periodontitis.

## Terms of Use & Disclaimer
- **Academic research only:** All data and scripts in this repository are provided strictly for non-commercial, basic academic research.
- **No clinical or diagnostic use:** Outputs from ThermoMPNN, AlphaMissense, and AlphaGenome v2 must **not** be used for clinical decision-making, diagnosis, or commercial applications, in accordance with the AlphaGenome and AlphaFold terms of service.

## Repository Structure
- `PRE_REGISTRO.md`: Pre-registered thresholds, rules, and train/external-test split (committed before analysis).
- `PASSAGEM_METODOLOGICA.md`: Methodological audit trail (`[EXECUTADO]`, `[FALHOU]`, limitations, dataset independence).
- `LEIAME.md`: Chronological research log.
- `manuscrito.md`: Full manuscript draft (English).
- `revisao_literatura.md`: Literature review with verified PMIDs/PMCIDs.
- `figuras/`: Figures 1–5 (`PNG` 300 dpi + `PDF`) and `LEGENDAS.md`.
- `scripts/`: Reproducible Python 3.12 pipeline (`baixar_dados.py`, `preparar_estruturas.py`, `consolidar.py`, `analise_estatistica.py`, `validacao_funcional_fenotipo.py`, `expressao_por_celula.py`, `alphagenome_score.py`, `figuras.py`).
- `estabilidade/`: Sanitized ThermoMPNN saturation mutagenesis outputs (`1K3B` and `AF-P53634-F1`).
- `variantes_CTSC.csv`, `descritores_estruturais.csv`, `validacao_funcional.csv`, `fenotipos_chitsamankhun2024.csv`, `alphagenome_nao_codificantes.csv`, `expressao_por_celula.csv`.
