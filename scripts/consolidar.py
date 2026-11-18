"""Consolida variantes missense do CTSC: ClinVar + gnomAD v4.1 (com homozigotos) + AlphaMissense + Literatura + 1K3B + ThermoMPNN.

Uso: py -I -X utf8 consolidar.py <pasta_papillon_lefevre>
Rótulos:
  P     = Pathogenic / Likely pathogenic no ClinVar (sem conflito e sem Benign)
  P-lit = descrita na literatura curada (sem conflito no ClinVar e gnomAD AF < 1e-4)
  B     = Benign / Likely benign no ClinVar, ou AF >= 1e-3, ou homozygote_count >= 1 (sem P/P-lit)
  U     = resto (VUS, conflitante, ou rara sem classificação)
"""
import csv
import re
import sys
from pathlib import Path

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
MIS = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def ler(p):
    with open(p, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def main():
    base = Path(sys.argv[1])
    brutos = base / "dados_brutos"
    est = base / "estabilidade"

    am = {r["protein_variant"]: r for r in ler(brutos / "alphamissense_CTSC.csv")}
    desc = {int(r["uniprot_pos"]): r for r in ler(base / "descritores_estruturais.csv")}

    # ThermoMPNN 1K3B: uniprot_pos = position + 25
    ddg_1k3b = {}
    for r in ler(est / "ThermoMPNN_inference_1K3B_mature_chainA.csv"):
        if r["wildtype"] != r["mutation"]:
            up = int(r["position"]) + 25
            ddg_1k3b[(up, r["mutation"])] = float(r["ddG_pred"])

    # ThermoMPNN AF-P53634-F1: uniprot_pos = position + 1
    ddg_af = {}
    for r in ler(est / "ThermoMPNN_inference_AF-P53634-F1.csv"):
        if r["wildtype"] != r["mutation"]:
            up = int(r["position"]) + 1
            ddg_af[(up, r["mutation"])] = float(r["ddG_pred"])

    var = {}
    for r in ler(brutos / "clinvar_CTSC.csv"):
        if "missense" not in r["molecular_consequence"]:
            continue
        for pc in r["protein_change"].split(","):
            pc = pc.strip()
            if MIS.match(pc) and pc in am:
                v = var.setdefault(pc, {})
                v.setdefault("clinvar", []).append(r)
                break

    for r in ler(brutos / "gnomad_CTSC.csv"):
        if r["consequence"] != "missense_variant":
            continue
        m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", r["hgvsp"])
        if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
            continue
        pc = AA3[m.group(1)] + m.group(2) + AA3[m.group(3)]
        if pc not in am:
            continue
        v = var.setdefault(pc, {})
        v["gnomad_af"] = max(float(r["af"] or 0), v.get("gnomad_af", 0))
        v["gnomad_ac"] = v.get("gnomad_ac", 0) + int(r["ac"] or 0)
        v["gnomad_hom"] = v.get("gnomad_hom", 0) + int(r["homozygote_count"] or 0)

    lit_path = base / "literatura_missense.csv"
    lit = {r["variante"]: r for r in ler(lit_path)} if lit_path.exists() else {}
    for pc in lit:
        if pc in am:
            var.setdefault(pc, {})

    linhas = []
    for pc, v in var.items():
        pos = int(MIS.match(pc).group(2))
        mut_aa = pc[-1]
        if pos == 1:
            continue  # exclui perdas de códon de início
        cls = [c["classification"] for c in v.get("clinvar", [])]
        rev = [c["review_status"] for c in v.get("clinvar", [])]
        cond = sorted({x for c in v.get("clinvar", []) for x in c["conditions"].split(";") if x})
        joined = " | ".join(cls)
        af = v.get("gnomad_af", 0.0)
        ac = v.get("gnomad_ac", 0)
        hom = v.get("gnomad_hom", 0)

        if re.search(r"(?<!Likely )Pathogenic|Likely pathogenic", joined) and "onflicting" not in joined and "enign" not in joined:
            rot = "P"
        elif pc in lit and "onflicting" not in joined and af < 1e-4:
            rot = "P-lit"
        elif ("enign" in joined and "athogenic" not in joined) or af >= 1e-3 or hom >= 1:
            rot = "B"
        else:
            rot = "U"

        coorte = ""
        if rot == "P":
            coorte = "calibracao"
        elif rot == "P-lit":
            coorte = "calibracao" if lit[pc]["coorte_temporal"] == "<=2015" else "teste_externo"

        d = desc.get(pos, {})
        d1k = ddg_1k3b.get((pos, mut_aa))
        daf = ddg_af.get((pos, mut_aa))

        linhas.append({
            "gene": "CTSC",
            "variante": pc,
            "posicao": pos,
            "dominio": d.get("dominio", ""),
            "na_1k3b": d.get("na_1k3b", "nao"),
            "clinvar": joined,
            "clinvar_revisao": " | ".join(rev),
            "condicoes": "; ".join(cond),
            "gnomad_af": f"{af:.2e}" if af else "",
            "gnomad_ac": ac if ac else "",
            "gnomad_hom": hom if hom else "",
            "alphamissense": am[pc]["am_pathogenicity"],
            "am_classe": am[pc]["am_class"],
            "rotulo": rot,
            "coorte_patogenica": coorte,
            "literatura_pmcids": lit[pc]["pmcids"] if pc in lit else "",
            "literatura_ano": lit[pc]["primeiro_ano_observado"] if pc in lit else "",
            "dist_cys258": d.get("dist_cys258", ""),
            "dist_sitio_ativo": d.get("dist_sitio_ativo", ""),
            "dist_cloreto": d.get("dist_cloreto", ""),
            "dist_interface_tetramero": d.get("dist_interface_tetramero", ""),
            "rsa_tetramero": d.get("rsa_tetramero", ""),
            "dissulfeto": d.get("dissulfeto", "nao"),
            "glicosilacao": d.get("glicosilacao", "nao"),
            "ddg_1k3b": round(d1k, 3) if d1k is not None else "",
            "ddg_af": round(daf, 3) if daf is not None else "",
            "clinvar_ids": " ".join(c["accession"] for c in v.get("clinvar", [])),
        })

    linhas.sort(key=lambda r: (r["posicao"], r["variante"]))
    out_csv = base / "variantes_CTSC.csv"
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)

    print(f"Total variantes missense consolidadas (sem M1*): {len(linhas)}")
    for rot in ("P", "P-lit", "B", "U"):
        sub = [r for r in linhas if r["rotulo"] == rot]
        in_1k = [r for r in sub if r["na_1k3b"] == "sim"]
        print(f"  {rot:5s}: {len(sub):3d} total, {len(in_1k):3d} na estrutura 1K3B")
    cal = [r for r in linhas if r["coorte_patogenica"] == "calibracao" and r["na_1k3b"] == "sim"]
    ext = [r for r in linhas if r["coorte_patogenica"] == "teste_externo" and r["na_1k3b"] == "sim"]
    hom_c = [r for r in linhas if r["rotulo"] not in ("P", "P-lit") and int(r["gnomad_hom"] or 0) >= 1 and r["na_1k3b"] == "sim"]
    ac2_c = [r for r in linhas if r["rotulo"] not in ("P", "P-lit") and int(r["gnomad_ac"] or 0) >= 2 and r["na_1k3b"] == "sim"]
    print(f"Na 1K3B — Calibração (P + P-lit<=2015): {len(cal)}; Teste externo (P-lit>2015): {len(ext)}")
    print(f"Na 1K3B — Controle homozigoto (hom>=1): {len(hom_c)}; Controle AC>=2: {len(ac2_c)}")


if __name__ == "__main__":
    main()
