# Structural Mechanisms of Pathogenic *CTSC* Variants in Papillon-Lefèvre Syndrome: A Pre-Registered Computational Study

**Authors:** Filipe Vono, Andre Naufel, Valeria Santos  
**Affiliation:** School of Dentistry, Centro Universitário de Belo Horizonte (UniBH), Belo Horizonte, MG, Brazil  
**Target Journal:** *Journal of Periodontal Research* / *Journal of Clinical Periodontology*  

> **Research-Only & Terms Notice:** This study uses exclusively public datasets and open computational tools for basic academic research. Predictions from ThermoMPNN and AlphaGenome v2 are not intended for clinical, diagnostic, or commercial use, in compliance with the AlphaGenome and AlphaFold terms of service.

---

## Abstract

**Background:** Loss-of-function variants in *CTSC* (encoding cathepsin C / dipeptidyl peptidase I, DPPI) cause Papillon-Lefèvre syndrome (PLS), Haim-Munk syndrome, and nonsyndromic aggressive periodontitis. While isolated case reports have modeled individual mutations, no systematic study has evaluated structural mechanisms against recessive population controls with pre-registered thresholds and independent external validation.

**Methods:** We consolidated 832 *CTSC* missense variants from ClinVar, gnomAD v4.1 (including homozygote counts), AlphaMissense, and curated literature (`n = 51` pathogenic variants in the mature enzyme; `20` pre-2016 calibration, `31` post-2015 external test). Using the experimental human DPPI tetramer crystal structure (PDB `1K3B`, 2.15 Å) and ThermoMPNN saturation mutagenesis (`ΔΔG`), we pre-registered distance and stability cutoffs before comparing pathogenic variants against gnomAD `AC ≥ 2` (`n = 234`) and homozygous (`homozygote_count ≥ 1`, `n = 4`) controls. Non-coding ClinVar variants (`n = 152`) were scored with AlphaGenome v2 in myeloid K562 cells (`EFO:0002067`).

**Results:** Pathogenic *CTSC* missense variants were strongly destabilizing compared to `AC ≥ 2` controls (median `ΔΔG` `1.796` vs. `0.762 kcal/mol`; `AUC = 0.740`, 95% CI `0.657–0.813`, `p = 8.33 × 10⁻⁸`) and enriched near the catalytic triad (`< 4.0 Å`: `15.7%` vs. `3.8%`, `p = 4.14 × 10⁻³`) and disulfide bonds (`9.8%` vs. `2.1%`, `p = 0.0188`). Conversely, tetramer interface proximity (`p = 1.00`), chloride ion proximity (`p = 0.665`), and N-glycosylation sites (`p = 1.00`) did not discriminate groups. AlphaMissense achieved `AUC = 0.819` (`p = 8.94 × 10⁻¹³`). In nested leave-one-out calibration (`c = 2.25 kcal/mol`), the four-term structural rule achieved `75.0%` sensitivity (`15/20`) and `74.4%` specificity (`174/234`), dropping to `58.1%` (`18/31`, 95% CI `40.8–73.6%`) on the independent post-2015 external test set and `75.0%` (`3/4`) on homozygous controls (`4/4 = 100%` for the two-term catalytic/stability rule). Functional validation confirmed `4/4` leukocyte DPPI-deficient mutants (`Q252L`, `R272P`, `G301S`, `Y347C`), whereas the exonic splicing mutant `A253T` (`ΔΔG = -0.43 kcal/mol`) was negative by protein structure yet strongly flagged by AlphaGenome v2 (`SPLICE_SITES raw_score = 0.8589`). Across `152` non-coding variants, AlphaGenome v2 separated pathogenic from benign variants with `AUC = 1.000` for splicing (`p = 1.55 × 10⁻⁸`) and `0.953` for RNA-seq (`p = 2.96 × 10⁻⁷`).

**Conclusions:** Thermodynamic destabilization and catalytic triad disruption account for most pathogenic missense *CTSC* variants, while exonic splicing variants require complementary nucleotide-level models. Pre-registration and independent post-2015 testing reveal realistic generalization bounds (`75.0%` vs. `58.1%`) that prevent overoptimistic classification claims.

---

## 1. Introduction

Papillon-Lefèvre syndrome (PLS) is a rare autosomal recessive disorder characterized by severe early-onset periodontitis leading to premature loss of primary and permanent dentition, accompanied by palmoplantar keratoderma (Toomes et al., 1999; Hart et al., 2000). PLS, Haim-Munk syndrome, and a subset of nonsyndromic aggressive periodontitis cases are caused by biallelic loss-of-function variants in *CTSC*, which encodes lysosomal cathepsin C (dipeptidyl peptidase I, DPPI; UniProt `P53634`, 463 amino acids). Mature DPPI functions as a homotetramer assembled from an exclusion domain (`25–134`), a heavy chain (`231–394`) and a light chain (`395–463`) following cleavage of the signal peptide (`1–24`) and propeptide (`135–230`), with a catalytic triad formed by `Cys258`, `His405`, and `Asn427` (Turk et al., 2001). Because DPPI activates neutrophil serine proteases (neutrophil elastase, cathepsin G, proteinase 3), loss of DPPI activity impairs periodontal host defense (Pham et al., 2004).

Although variant databases (Nagy et al., 2014; Chitsamankhun et al., 2024) and isolated homology-modeling case reports (Moghaddasian et al., 2014) have cataloged dozens of *CTSC* missense mutations, no systematic study has evaluated structural mechanisms across all known variants using an experimental crystal structure, thermodynamic stability predictions (`ΔΔG`), and recessive population controls. Moreover, computational variant studies often suffer from circularity when cutoffs are chosen after inspecting the data or when variants of uncertain significance (VUS) overlap with control sets. Here, we apply a pre-registered, non-circular framework to classify *CTSC* missense and non-coding variants.

---

## 2. Materials and Methods

### 2.1 Pre-registration and dataset independence
All distance thresholds, stability cutoffs, rule forms, and calibration/test partitions were written and committed to version control (`PRE_REGISTRO.md`, commit `574901d`) prior to any statistical comparison between pathogenic and control groups. To ensure an independent external test set, ClinVar pathogenic/likely pathogenic (`P`) variants and literature-curated (`P-lit`, `gnomAD AF < 10⁻⁴`, no ClinVar conflict) variants published up to 2015 were assigned to calibration (`n = 20` in `1K3B`), whereas `P-lit` variants published after 2015 (`n = 31` in `1K3B`) were held out strictly for external testing. Because PLS is autosomal recessive, heterozygous carriers of pathogenic alleles appear in population databases; we therefore excluded all `P`/`P-lit` variants from controls and defined a stringent homozygous control set (`gnomAD homozygote_count ≥ 1`, `n = 4` in `1K3B`) alongside a sensitivity control set (`AC ≥ 2`, `n = 234` in `1K3B`). VUS (`U`, `n = 569` in `1K3B`) were stratified into mutually disjoint subgroups (`AC = 0`, `AC = 1`, `AC ≥ 2`).

### 2.2 Public data, experimental structure, and stability (`ΔΔG`)
ClinVar (`522` total records; `153` non-coding GRCh38 variants), gnomAD v4.1 (`dataset=gnomad_r4`, canonical transcript `ENST00000227266`, max `AN = 1,614,196`, `3,456` variants), and AlphaMissense (`8,797` substitutions for `P53634`) were retrieved via public APIs. Structural descriptors were computed on the experimental human DPPI biological tetramer (PDB `1K3B`, `2.15 Å`, `1K3B-assembly1.cif`, `352` resolved residues per monomer mapped 100% to UniProt `P53634` via `uniprot_pos = auth_seq_id + 24`): minimum heavy-atom distance to the catalytic triad (`Cys258`, `His405`, `Asn427`), tetramer interface, chloride ions (`CL`), relative solvent accessibility (Shrake-Rupley RSA), disulfide bonds (`30–118`, `54–136`, `255–298`, `291–331`, `321–337`), and N-glycosylation sites (`29`, `53`, `119`, `276`). Saturation mutagenesis `ΔΔG` was computed locally using ThermoMPNN (`thermoMPNN_default.pt`) on the unified mature monomer (`1K3B_mature_chainA.pdb`, `7,040` predictions) and full-length AlphaFold DB model (`AF-P53634-F1.pdb`, `9,260` predictions), where `ΔΔG > 0` denotes destabilization (Dieckhaus et al., 2024).

### 2.3 Non-coding variants and statistical analysis
AlphaGenome v2 was first benchmarked for *CTSC* RNA-seq expression against `GAPDH` (positive) and `ALB` (negative) across five supported cell/tissue conditions; myeloid `K562` (`EFO:0002067`, observed tracks, *CTSC* p99 = `9.7511` vs. `ALB = 0.0600`, `162.5× ALB`) was selected as primary, with `HL-60` (`159.7× ALB`) and skin fibroblast (`CL:0011022`, `39.5× ALB`) as sensitivity conditions. Non-coding ClinVar variants (`n = 152`) and the exonic splicing variant `A253T` (`c.757G>A`) were scored with AlphaGenome v2 (`SPLICE_SITES`, `SPLICE_SITE_USAGE`, `SPLICE_JUNCTIONS`, `RNA_SEQ`). Statistical comparisons used two-sided `scipy.stats.mannwhitneyu`, two-sided Fisher's exact test (`scipy.stats.fisher_exact`), Wilson 95% confidence intervals for proportions, and 2,000-replicate bootstrap 95% CIs for AUC and likelihood ratios.

---

## 3. Results

### 3.1 Thermodynamic destabilization and catalytic proximity drive *CTSC* pathogenicity
Across all `51` pathogenic missense variants in `1K3B` versus `234` `gnomAD AC ≥ 2` controls, predicted destabilization (`ΔΔG`) was markedly elevated (median `1.796` vs. `0.762 kcal/mol`; `AUC = 0.740`, 95% CI `0.657–0.813`, `p = 8.33 × 10⁻⁸`; `AF-P53634` `AUC = 0.741`, `p = 6.83 × 10⁻⁸`). Pathogenic variants were significantly buried (`RSA` median `0.000` vs. `0.006`, `p = 5.92 × 10⁻³`) and enriched for `ΔΔG ≥ 1.5 kcal/mol` (`31/51 = 60.8%` vs. `62/234 = 26.5%`, Fisher `p = 5.37 × 10⁻⁶`), catalytic triad proximity `< 4.0 Å` (`8/51 = 15.7%` vs. `9/234 = 3.8%`, `p = 4.14 × 10⁻³`), and disulfide bond disruption (`5/51 = 9.8%` vs. `5/234 = 2.1%`, `p = 0.0188`). By contrast, tetramer interface proximity `< 4.0 Å` (`11.8%` vs. `11.5%`, `p = 1.00`), chloride ion proximity `< 5.0 Å` (`3.9%` vs. `3.0%`, `p = 0.665`), and N-glycosylation (`0.0%` vs. `0.4%`, `p = 1.00`) showed no enrichment. AlphaMissense strongly separated groups (`0.853` vs. `0.208`; `AUC = 0.819`, 95% CI `0.753–0.877`, `p = 8.94 × 10⁻¹³`).

### 3.2 Nested cross-validation, independent external test, and homozygous controls
In nested leave-one-out calibration (`n = 20` vs. `234`), Youden's index selected `ΔΔG ≥ 2.25 kcal/mol` (`range 2.25–2.75`), yielding `75.0%` sensitivity (`15/20`, 95% CI `53.1–88.8%`) and `74.4%` specificity (`174/234`, 95% CI `68.4–79.5%`; `LR+ = 2.93`, `LR− = 0.34`) for the four-term rule. On the independent post-2015 external test set (`n = 31`), sensitivity was `58.1%` (`18/31`, 95% CI `40.8–73.6%`). On homozygous gnomAD controls (`L71M`, `T233A`, `M420L`, `I453V`), `ΔΔG` ranged from `-0.138` to `0.289 kcal/mol`, giving `75.0%` (`3/4`) specificity for the four-term rule and `100.0%` (`4/4`) for the two-term catalytic/stability rule (`L71M` lies `3.57 Å` from the tetramer interface). Across disjoint VUS subgroups (`ΔΔG ≥ 1.5 kcal/mol`), `44.5%` (`49/110`) of gnomAD-absent VUS, `42.6%` (`98/230`) of `AC = 1` VUS, and `41.5%` (`95/229`) of `AC ≥ 2` VUS were rule-positive.

### 3.3 Functional DPPI validation and splicing rescue of `A253T`
All four leukocyte DPPI-deficient missense mutants from Pham et al. (2004) resolved in `1K3B` (`Q252L`, `R272P`, `G301S`, `Y347C`) were correctly predicted as pathogenic (`4/4`), whereas the experimentally proven exonic splicing mutant `A253T` (`c.757G>A`; Abdel-Hamid et al., 2024) was predicted benign at the protein level (`ΔΔG = -0.43 kcal/mol`, `7.52 Å` from catalytic triad) but strongly positive by AlphaGenome v2 (`SPLICE_SITES raw_score = 0.8589`, `SPLICE_SITE_USAGE = 0.8645`).

### 3.4 Non-coding ClinVar variants
Across `152` non-coding *CTSC* variants scored in `K562`, AlphaGenome v2 achieved `AUC = 1.000` for `SPLICE_SITES` (`12 P` vs. `105 B`, median `0.953` vs. `0.006`, `p = 1.55 × 10⁻⁸`) and `AUC = 0.953` for `RNA_SEQ` (`p = 2.96 × 10⁻⁷`), prioritizing canonical splice VUS (`c.173-2A>T/G`, `c.641+1G>T`, `c.889+2T>A`, `c.757+1G>A`, `c.890-2A>G`) and deep intronic candidates (`c.318+6826G>T`, `0.8932`).

---

## 4. Discussion and Limitations

Our pre-registered analysis demonstrates that protein destabilization (`ΔΔG`) and catalytic triad disruption are the principal structural mechanisms of *CTSC* missense pathogenicity, while tetramer interface proximity alone is non-discriminative—an informative negative result that prevents overfitting. Importantly, evaluating a held-out post-2015 external cohort revealed a realistic sensitivity drop (`75.0%` to `58.1%`), underscoring that leave-one-out cross-validation alone overestimates generalization. Finally, the mechanistic dissociation of `A253T` highlights the necessity of pairing 3D stability models with nucleotide-level splicing predictors such as AlphaGenome v2.

Limitations include: (1) small sample size of homozygous missense controls in `1K3B` (`n = 4`); (2) reliance on a single `ΔΔG` predictor (ThermoMPNN); (3) timeout of the Ensembl VEP REST endpoint preventing direct SpliceAI comparison in this run; and (4) lack of laboratory validation for prioritized non-coding VUS.
