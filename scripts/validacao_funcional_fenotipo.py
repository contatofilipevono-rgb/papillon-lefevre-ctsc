"""Extrai fenótipos (PLS vs. periodontite agressiva isolada vs. HMS) da Tabela 1 de Chitsamankhun 2024 (PMC11337848)
e cruza com a regra estrutural e com variantes de validação funcional da literatura.

Uso: py -I -X utf8 validacao_funcional_fenotipo.py <pasta_papillon_lefevre>
"""
import csv
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
P3 = re.compile(r"p\.([A-Z][a-z]{2})(\d{1,3})([A-Z][a-z]{2})")


def regra(r, cg=1.5):
    return (
        float(r["dist_sitio_ativo"]) < 4.0
        or float(r["dist_interface_tetramero"]) < 4.0
        or r["dissulfeto"] == "sim"
        or float(r["ddg_1k3b"]) >= cg
    )


def main():
    base = Path(sys.argv[1])
    var_map = {r["variante"]: r for r in csv.DictReader(open(base / "variantes_CTSC.csv", encoding="utf-8-sig"))}

    # 1. Fenótipo a partir de Chitsamankhun 2024 (PMC11337848)
    root = ET.fromstring((base / "dados_brutos" / "literatura" / "PMC11337848.xml").read_bytes())
    tw = root.find(".//table-wrap")
    fen_rows = []
    for tr in tw.findall(".//tr"):
        cells = [re.sub(r"\s+", " ", " ".join(td.itertext())).strip() for td in tr]
        if len(cells) < 11:
            continue
        m = P3.search(cells[1])
        if not m:
            continue
        wt3, pos_s, mut3 = m.groups()
        if wt3 not in AA3 or mut3 not in AA3:
            continue
        var = f"{AA3[wt3]}{pos_s}{AA3[mut3]}"
        if var not in var_map or var_map[var]["na_1k3b"] != "sim":
            continue
        r = var_map[var]
        pls_flag = cells[6] if len(cells) > 6 else ""
        perio = cells[10] if len(cells) > 10 else ""
        other = cells[11] if len(cells) > 11 else ""
        if "(-)" in pls_flag and "aggressive" in perio.lower():
            fen = "Periodontite_agressiva"
        elif "haim" in other.lower() or "arachnodactyly" in other.lower() or "acro-osteolysis" in other.lower():
            fen = "Haim_Munk"
        else:
            fen = "PLS"
        fen_rows.append({
            "variante": var,
            "fenotipo": fen,
            "heranca": cells[2],
            "ddg_1k3b": r["ddg_1k3b"],
            "dist_sitio_ativo": r["dist_sitio_ativo"],
            "dissulfeto": r["dissulfeto"],
            "regra_c1_5": "sim" if regra(r, 1.5) else "nao",
            "regra_c2_25": "sim" if regra(r, 2.25) else "nao",
            "alphamissense": r["alphamissense"],
        })

    out_fen = base / "fenotipos_chitsamankhun2024.csv"
    with open(out_fen, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(fen_rows[0].keys()))
        w.writeheader()
        w.writerows(fen_rows)

    print(f"Fenótipos extraídos de PMC11337848 (na 1K3B): {len(fen_rows)} registros")
    for fen in sorted({x["fenotipo"] for x in fen_rows}):
        sub = [x for x in fen_rows if x["fenotipo"] == fen]
        n15 = sum(x["regra_c1_5"] == "sim" for x in sub)
        n225 = sum(x["regra_c2_25"] == "sim" for x in sub)
        print(f"  {fen:24s}: n={len(sub):2d} | regra(1.5)={n15}/{len(sub)} ({n15/len(sub):.0%}) | regra(2.25)={n225}/{len(sub)} ({n225/len(sub):.0%})")

    # 2. Variantes com atividade enzimática / ensaio funcional na literatura
    func_lit = [
        ("L172P", "Hamon et al. 2019 (PMC6801622)", "expressao recombinante", "perda de processamento e atividade DPPI", "sim"),
        ("Q252L", "Pham et al. 2004 (PMID 15585850)", "leucocitos / PMNs", "ausencia de atividade DPPI", "sim"),
        ("A253T", "Abdel-Hamid et al. 2024 (PMID 38104461)", "splicing (minigene/RNA)", "altera splicing (exon skipping)", "sim"),
        ("R272P", "Pham et al. 2004 (PMID 15585850)", "leucocitos / PMNs", "ausencia de atividade DPPI", "sim"),
        ("G301S", "Pham et al. 2004 (PMID 15585850)", "leucocitos / PMNs", "ausencia de atividade DPPI", "sim"),
        ("Y347C", "Pham et al. 2004 (PMID 15585850)", "leucocitos / PMNs", "ausencia de atividade DPPI", "sim"),
    ]
    func_rows = []
    print("\n=== VALIDAÇÃO FUNCIONAL NA LITERATURA ===")
    for var, ref, ensaio, desfecho, perda in func_lit:
        r = var_map.get(var, {})
        na_1k = r.get("na_1k3b", "nao")
        p15 = ("sim" if regra(r, 1.5) else "nao") if na_1k == "sim" else "fora_1k3b"
        p225 = ("sim" if regra(r, 2.25) else "nao") if na_1k == "sim" else "fora_1k3b"
        func_rows.append({
            "variante": var, "referencia": ref, "ensaio": ensaio, "desfecho": desfecho,
            "perda_funcao": perda, "na_1k3b": na_1k,
            "ddg_1k3b": r.get("ddg_1k3b", ""), "ddg_af": r.get("ddg_af", ""),
            "dist_sitio_ativo": r.get("dist_sitio_ativo", ""),
            "regra_c1_5": p15, "regra_c2_25": p225,
            "alphamissense": r.get("alphamissense", ""),
        })
        print(f"  {var:6s} | {ensaio:24s} | ddG_1K3B={r.get('ddg_1k3b','NA'):6s} | ddG_AF={r.get('ddg_af','NA'):6s} | d_act={r.get('dist_sitio_ativo','NA'):5s} | regra1.5={p15} | regra2.25={p225}")

    with open(base / "validacao_funcional.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(func_rows[0].keys()))
        w.writeheader()
        w.writerows(func_rows)


if __name__ == "__main__":
    main()
