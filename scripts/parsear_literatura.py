"""Extrai variantes missense do CTSC a partir dos XMLs do PMC (Nagy 2014, Chitsamankhun 2024 e relatos)
e valida cada aminoácido de referência contra a sequência UniProt P53634.

Uso: py -I -X utf8 parsear_literatura.py <pasta_papillon_lefevre>
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
P3 = re.compile(r"\bp\.([A-Z][a-z]{2})(\d{1,3})([A-Z][a-z]{2})\b")
P1 = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY])(\d{2,3})([ACDEFGHIKLMNPQRSTVWY])\b")


def parse_fasta(path):
    return "".join(l.strip() for l in path.read_text(encoding="utf-8").splitlines() if l.strip() and not l.startswith(">"))


def xml_year(root):
    for y in root.findall(".//pub-date/year"):
        if y.text and y.text.strip().isdigit():
            return int(y.text.strip())
    return 0


def main():
    base = Path(sys.argv[1])
    seq = parse_fasta(base / "dados_brutos" / "uniprot_P53634.fasta")
    lit_dir = base / "dados_brutos" / "literatura"

    # Mapeia variante -> lista de (pmcid, year)
    encontradas = {}
    for xml_path in sorted(lit_dir.glob("*.xml")):
        pmcid = xml_path.stem
        try:
            root = ET.fromstring(xml_path.read_bytes())
        except ET.ParseError:
            continue
        yr = xml_year(root)
        txt = " ".join(root.itertext())
        for wt3, pos_s, mut3 in P3.findall(txt):
            pos = int(pos_s)
            if wt3 in AA3 and mut3 in AA3 and 1 <= pos <= len(seq):
                wt1, mut1 = AA3[wt3], AA3[mut3]
                if wt1 != mut1 and seq[pos - 1] == wt1:
                    var = f"{wt1}{pos}{mut1}"
                    encontradas.setdefault(var, []).append((pmcid, yr))
        for wt1, pos_s, mut1 in P1.findall(txt):
            pos = int(pos_s)
            if 1 <= pos <= len(seq) and wt1 != mut1 and seq[pos - 1] == wt1:
                var = f"{wt1}{pos}{mut1}"
                encontradas.setdefault(var, []).append((pmcid, yr))

    rows = []
    for var, refs in sorted(encontradas.items(), key=lambda x: (int(x[0][1:-1]), x[0])):
        pos = int(var[1:-1])
        anos = sorted({yr for _, yr in refs if yr > 0})
        pmcids = sorted({p for p, _ in refs})
        min_ano = min(anos) if anos else 0
        # Se aparece em Nagy 2014 (PMC4049362) ou artigo <= 2015, é <=2015; caso contrário >2015
        coorte = "<=2015" if (min_ano and min_ano <= 2015) else ">2015"
        rows.append({
            "gene": "CTSC",
            "variante": var,
            "posicao": pos,
            "primeiro_ano_observado": min_ano or "",
            "coorte_temporal": coorte,
            "pmcids": ";".join(pmcids),
            "n_artigos": len(pmcids),
        })

    out_csv = base / "literatura_missense.csv"
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["gene"])
        w.writeheader()
        w.writerows(rows)
    n_pre = sum(1 for r in rows if r["coorte_temporal"] == "<=2015")
    n_post = sum(1 for r in rows if r["coorte_temporal"] == ">2015")
    print(f"Total variantes missense validadas na literatura OA: {len(rows)} (<=2015: {n_pre}, >2015: {n_post})")


if __name__ == "__main__":
    main()
