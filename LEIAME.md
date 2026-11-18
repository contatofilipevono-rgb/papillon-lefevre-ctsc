# Diário de pesquisa — *CTSC* / Síndrome de Papillon-Lefèvre

**Autores:** Filipe Vono, Andre Naufel, Valeria Santos (Odontologia, UniBH).  
**Modelo metodológico:** `agenesia/` (PAX9/MSX1), corrigindo independência de conjuntos, pré-registro de cortes, estrutura experimental e controle recessivo.

---

## 18/11/2026 — Etapa 1: Ineditismo, dados, estruturas e pré-registro

1. **Ineditismo:**
   - Europe PMC (`busca_literatura_saida.txt`): relatos de caso com modelagem por homologia de 1–2 mutações (Moghaddasian 2014, *Gene*, PMID 24374475), revisão (Chitsamankhun 2024, *J Transl Med*, PMID 39164687), bancos de variantes (Nagy 2014, *Mol Genet Genomic Med*; Hart 2000; Toomes 1999) e ensaios isolados (Abdel-Hamid 2024).
   - Confirmado em busca web (bioRxiv, medRxiv, Google Scholar): não há estudo sistemático de variantes *missense* do *CTSC* com estrutura experimental + $\Delta\Delta G$ + controle populacional do gnomAD.
2. **Dados públicos (`scripts/baixar_dados.py`, `scripts/clinvar_nao_codificantes.py`):**
   - ClinVar: 522 registros totais; 153 variantes não codificantes com coordenadas GRCh38 (`chr11`).
   - gnomAD v4.1 (`dataset=gnomad_r4`, transcrito canônico `ENST00000227266`, max AN = 1.614.196): 3.456 variantes, das quais 79 têm `homozygote_count >= 1`.
   - AlphaMissense: 8.797 substituições para `P53634` (463 × 19).
   - UniProt `P53634`: peptídeo sinal 1–24; domínio de exclusão 25–134; propeptídeo 135–230; cadeia pesada 231–394; cadeia leve 395–463; tríade catalítica Cys258, His405, Asn427; dissulfetos 30–118, 54–136, 255–298, 291–331, 321–337; N-glicosilação 29, 53, 119, 276.
   - Estrutura experimental: PDB `1K3B` (2,15 Å) + tetrâmero biológico `1K3B-assembly1.cif` (4 íons `CL`, 352 resíduos por monômero mapeados 100% em `uniprot_pos = auth_seq_id + 24`).
3. **Estabilidade (`ThermoMPNN`):**
   - `1K3B_mature_chainA.pdb` (352 resíduos, 7.040 predições, 0 discordâncias com `P53634` usando `position + 25`).
   - `AF-P53634-F1.pdb` (463 resíduos, 9.260 predições).
   - Caminhos locais sanitizados da coluna `pdb`.
   - Convenção de sinal confirmada em Dieckhaus et al. (2024, *PNAS*): $\Delta\Delta G > 0$ = desestabilizante.
4. **Expressão prevista no AlphaGenome v2 (`scripts/expressao_por_celula.py`):**
   - Condição principal selecionada: **K562 (`EFO:0002067`)**, trilha observada, `CTSC` p99 = 9,7511 vs. `ALB` 0,0600 (**162,5× ALB**), `GAPDH` = 743,87.
   - Sensibilidade: **HL-60 (`EFO:0002793`)** (159,7× ALB) e **fibroblasto de pele (`CL:0011022`)** (39,5× ALB, observado).
5. **Pré-registro:** `PRE_REGISTRO.md` commitado antes de qualquer comparação entre patogênicas e controles.
