# Passagem metodológica — *CTSC* / Síndrome de Papillon-Lefèvre

*Redigido em 18/11/2026 com base exclusivamente nos scripts e saídas executados em `papillon-lefevre/`.*

## Status dos itens

- **[EXECUTADO] Verificação de ineditismo:** Europe PMC (`busca_literatura_saida.txt`), bioRxiv, medRxiv e busca web. Nenhuma análise sistemática combinando estrutura experimental do CTSC, $\Delta\Delta G$ e controle populacional do gnomAD foi encontrada.
- **[EXECUTADO] Pré-registro commitado antes da análise:** `PRE_REGISTRO.md` commitado em `574901d` antes de cruzar patogênicas com controles.
- **[EXECUTADO] Dados públicos (`scripts/baixar_dados.py`, `scripts/clinvar_nao_codificantes.py`):**
  - ClinVar: 522 registros totais; 153 variantes não codificantes GRCh38 (`chr11`).
  - gnomAD v4.1 (`dataset=gnomad_r4`, `ENST00000227266`, max AN = 1.614.196): 3.456 variantes, 79 com `homozygote_count >= 1` (4 missense dentro da `1K3B`, 1 no propeptídeo `I153T`).
  - AlphaMissense: 8.797 substituições para `P53634`.
  - Consolidação (`scripts/consolidar.py` → `variantes_CTSC.csv`): 832 missense sem `M1*` (`18 P`, `35 P-lit`, `7 B`, `772 U`; na `1K3B`: `16 P`, `35 P-lit`, `5 B`, `569 U`).
- **[EXECUTADO] Estrutura experimental (`scripts/preparar_estruturas.py`):**
  - PDB `1K3B` (2,15 Å) + tetrâmero biológico `1K3B-assembly1.cif` (4 íons `CL`, 352 resíduos/monômero mapeados 100% em `uniprot_pos = auth_seq_id + 24`).
- **[EXECUTADO] Estabilidade (`ThermoMPNN`):**
  - `1K3B_mature_chainA.pdb` (7.040 predições, 0 discordâncias com `P53634` em `position + 25`) e `AF-P53634-F1.pdb` (9.260 predições). Caminhos locais sanitizados.
  - Convenção de sinal confirmada em Dieckhaus et al. (2024, *PNAS*): $\Delta\Delta G > 0$ = desestabilizante.
- **[EXECUTADO] Análise estatística pré-registrada (`scripts/analise_estatistica.py`):**
  - `scipy.stats.mannwhitneyu` (com correção de empates), `scipy.stats.fisher_exact`, IC 95% de Wilson e bootstrap (2.000 réplicas).
  - Desestabilização ($\Delta\Delta G$) domina: AUC = 0,740 (IC95% 0,657–0,813; `p = 8,33e-08`), `60,8%` vs. `26,5%` em `>= 1,5 kcal/mol` (`p = 5,37e-06`). Sítio ativo `< 4,0 Å`: `15,7%` vs. `3,8%` (`p = 4,14e-03`). Dissulfetos: `9,8%` vs. `2,1%` (`p = 0,0188`).
  - **Resultados negativos:** interface do tetrâmero `< 4,0 Å` (`11,8%` vs. `11,5%`, `p = 1,00`), íon cloreto `< 5,0 Å` (`3,9%` vs. `3,0%`, `p = 0,665`) e N-glicosilação (`0%` vs. `0,4%`, `p = 1,00`) não discriminam.
  - **Calibração LOO vs. Teste externo independente:** LOO na calibração (`n = 20` vs. `234`) escolheu corte `2,25 kcal/mol` (sens. `75,0%`, espec. `74,4%`). No teste externo independente (`P-lit > 2015`, `n = 31`), a sensibilidade caiu para `58,1%` (IC95% `40,8–73,6%`). Nos homozigotos (`n = 4`), especificidade `75,0%` (regra 4 termos) e `100,0%` (regra 2 termos).
- **[EXECUTADO] Validação funcional e fenótipo (`scripts/validacao_funcional_fenotipo.py`):**
  - 4/4 mutantes com perda de atividade DPPI em leucócitos na `1K3B` (Pham 2004: `Q252L`, `R272P`, `G301S`, `Y347C`) foram previstos como positivos.
  - `A253T` (`c.757G>A`, Abdel-Hamid 2024) foi negativa na regra proteica (`ΔΔG = -0,43`), mas positiva no AlphaGenome v2 (`SPLICE_SITES raw_score = 0,8589`).
  - Fenótipos extraídos de Chitsamankhun 2024 (`PMC11337848`): `61` PLS e `1` periodontite agressiva isolada na `1K3B`.
- **[EXECUTADO] Variantes não codificantes (`scripts/expressao_por_celula.py`, `scripts/alphagenome_score.py`):**
  - K562 (`EFO:0002067`) selecionada como condição mieloide principal (`CTSC` `162,5× ALB`).
  - AlphaGenome v2 (`n = 152`): `SPLICE_SITES` AUC = `1,000` (`12 P` vs. `105 B`, `p = 1,55e-08`), `RNA_SEQ` AUC = `0,953` (`p = 2,96e-07`).
- **[FALHOU] SpliceAI via Ensembl VEP REST (`scripts/spliceai_vep.py`):** `rest.ensembl.org` deu timeout nesta rede (`WinError 10060`). Documentado honestamente; não foi usado.
- **[FALHOU] LOVD:** bloqueado por IP (conforme registrado no protocolo inicial).
- **[EXECUTADO] Figuras 1–5 (`scripts/figuras.py`):** PNG 300 dpi + PDF em `figuras/`.
