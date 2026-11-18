# Pré-registro metodológico — *CTSC* / Síndrome de Papillon-Lefèvre

**Data:** 2026-11-18  
**Status:** Pré-registrado e commitado **antes** de qualquer comparação estatística entre variantes patogênicas e controles populacionais.

---

## 1. Pergunta e desenho

Avaliar os mecanismos estruturais pelos quais variantes *missense* no gene *CTSC* (catepsina C / dipeptidil peptidase I, UniProt `P53634`, 463 aa) causam síndrome de Papillon-Lefèvre (PLS), síndrome de Haim-Munk (HMS) ou periodontite agressiva, e testar se uma regra estrutural simples discrimina variantes patogênicas de variantes populacionais em uma doença autossômica recessiva.

## 2. Estrutura e estabilidade

- **Estrutura experimental principal:** PDB `1K3B` (resolução 2,15 Å; tetrâmero biológico `1K3B-assembly1.cif`).
  - Cadeia A: domínio de exclusão (`1–119` PDB = `25–143` UniProt).
  - Cadeia B: cadeia pesada (`207–370` PDB = `231–394` UniProt).
  - Cadeia C: cadeia leve (`371–439` PDB = `395–463` UniProt).
  - Mapeamento verificado: `uniprot_pos = auth_seq_id + 24` (352 resíduos, 100% de identidade com `P53634`).
- **Estrutura complementar:** AlphaFold DB (`AF-P53634-F1`), apenas para cobertura do propeptídeo (`144–230`) e análise de sensibilidade.
- **Estabilidade ($\Delta\Delta G$):** ThermoMPNN (`thermoMPNN_default.pt`) sobre `1K3B_mature_chainA.pdb`. Convenção confirmada em Dieckhaus et al. (2024, *PNAS*): $\Delta\Delta G > 0$ indica mutação **desestabilizante**.

## 3. Definição dos grupos e independência dos conjuntos

### 3.1 Variantes patogênicas
- **ClinVar (`P`):** *Pathogenic* ou *Likely pathogenic* sem classificação conflitante e sem termo *Benign*.
- **Literatura (`P-lit`):** variantes *missense* descritas em pacientes com PLS, HMS ou periodontite agressiva em artigos revisados por pares, sem conflito no ClinVar e com `gnomAD AF < 1e-4`.
- **Partição calibração vs. teste externo:**
  - **Calibração:** variantes `P` do ClinVar + `P-lit` publicadas até 2015 (`coorte_temporal <= 2015`).
  - **Teste externo independente:** variantes `P-lit` publicadas **após 2015** (`coorte_temporal > 2015`) que não estejam no conjunto de calibração. Não serão usadas na seleção de cortes.

### 3.2 Controles populacionais (doença recessiva)
- Exclusão obrigatória de qualquer variante `P` ou `P-lit`.
- **Controle primário (forte):** variantes *missense* no gnomAD v4.1 (na região estruturada da `1K3B`) observadas em **homozigose (`homozygote_count >= 1`)**, pois homozigotos na população geral são quase certamente benignos em uma doença recessiva grave de início precoce.
- **Controle de sensibilidade:** variantes com **`AC >= 2`** no gnomAD v4.1.

### 3.3 Variantes incertas (VUS / `U`)
- As VUS serão analisadas de forma **descritiva** e estratificadas em subgrupos disjuntos:
  1. Ausentes do gnomAD (`AC = 0`);
  2. Raras no gnomAD (`AC = 1`, `hom = 0`, fora da calibração);
  3. Recorrentes/homozigotas no gnomAD (`AC >= 2` ou `hom >= 1`, sobrepostas ao controle).

## 4. Descritores, cortes e regra pré-registrada

### 4.1 Cortes individuais pré-registrados
- **Sítio ativo (`dist_sitio_ativo`):** distância mínima de átomos pesados à tríade catalítica (`Cys258`, `His405`, `Asn427`) `< 4,0 Å` (sensibilidade: `3,5 Å` e `5,0 Å`).
- **Interface do tetrâmero (`dist_interface_tetramero`):** distância mínima a outro monômero no tetrâmero biológico `< 4,0 Å` (sensibilidade: `5,0 Å`).
- **Íon cloreto (`dist_cloreto`):** distância mínima ao íon `CL` `< 5,0 Å`.
- **Sítios covalentes:** participação em ponte dissulfeto (`30, 54, 118, 136, 255, 291, 298, 321, 331, 337`) ou sítio de N-glicosilação (`29, 53, 119, 276`).
- **Desestabilização ($\Delta\Delta G$):** grade de cortes $\{0{,}5;\; 1{,}0;\; 1{,}5;\; 2{,}0\}$ kcal/mol; seleção do corte ótimo por índice de Youden dentro de **validação cruzada leave-one-out (LOO) aninhada** na grade $0{,}00;\; 0{,}25;\; \dots;\; 3{,}00$ kcal/mol.
- **Splicing em *missense*:** `SpliceAI max DS >= 0,2` (moderado) e `>= 0,5` (forte).

### 4.2 Forma da regra
- **Regra estrutural completa (primária):** positiva se `dist_sitio_ativo < 4,0 Å` **OU** `dist_interface_tetramero < 4,0 Å` **OU** `dissulfeto == "sim"` **OU** `ddG >= c`.
- **Regra minimalista (secundária, comparável à agenesia):** positiva se `dist_sitio_ativo < 4,0 Å` **OU** `ddG >= c`.

## 5. Análises estatísticas e validação

- **Comparações contínuas:** `scipy.stats.mannwhitneyu` (bicaudal, com correção de empates) e AUC.
- **Comparações categóricas:** `scipy.stats.fisher_exact` (bicaudal).
- **Incerteza:** intervalos de confiança de 95% de **Wilson** para proporções (sensibilidade, especificidade) e **bootstrap** (2.000 réplicas) para LR+, LR− e AUC.
- **Validação funcional:** comparação da regra com atividade enzimática residual de DPPI medida na literatura (Pham et al. 2004; Hamon et al.; estudos de urina de 2022) e análise descritiva por fenótipo (PLS vs. HMS vs. periodontite agressiva).
- **Variantes não codificantes:** AlphaGenome v2 (condição primária `K562 EFO:0002067`, selecionada por expressão de `CTSC` 162,5× `ALB`; sensibilidade `HL-60` e `fibroblasto de pele`) e SpliceAI via Ensembl VEP REST.
