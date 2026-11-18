"""Prepara estruturas do CTSC (PDB 1K3B + tetrâmero biológico 1K3B-assembly1 + AlphaFold DB)
e calcula descritores estruturais por resíduo (numeração UniProt P53634).

Saídas:
  - estruturas/1K3B_mature_chainA.pdb (cadeias A+B+C unificadas para ThermoMPNN)
  - estruturas/1K3B_residuo_mapa.csv (índice 0-based -> posição UniProt -> AA)
  - descritores_estruturais.csv (distância ao sítio ativo, cloreto, interface do tetrâmero, RSA, dissulfeto, glicosilação, domínio)

Uso: <thermompnn_python> preparar_estruturas.py <pasta_papillon_lefevre>
"""
import csv
import math
import sys
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Tien et al. 2013 theoretical max ASA (Å^2)
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

ACTIVE_SITE = {258, 405, 427}
DISULFIDE = {30, 54, 118, 136, 255, 291, 298, 321, 331, 337}
GLYCOSYLATION = {29, 53, 119, 276}


def dominio_uniprot(pos):
    if 1 <= pos <= 24:
        return "peptideo_sinal"
    if 25 <= pos <= 134:
        return "exclusao"
    if 135 <= pos <= 230:
        return "propeptideo"
    if 231 <= pos <= 394:
        return "cadeia_pesada"
    if 395 <= pos <= 463:
        return "cadeia_leve"
    return "desconhecido"


def parse_fasta(path):
    lines = [l.strip() for l in path.read_text(encoding="utf-8").splitlines() if l.strip() and not l.startswith(">")]
    return "".join(lines)


def heavy_coords(res):
    return [tuple(a.get_coord()) for a in res.get_atoms() if a.element != "H"]


def min_dist_coords(c1, c2):
    if not c1 or not c2:
        return None
    return min(math.dist(a, b) for a in c1 for b in c2)


def main():
    base = Path(sys.argv[1])
    brutos = base / "dados_brutos"
    est_dir = base / "estruturas"
    est_dir.mkdir(parents=True, exist_ok=True)

    seq = parse_fasta(brutos / "uniprot_P53634.fasta")
    assert len(seq) == 463, f"Comprimento inesperado: {len(seq)}"

    parser = MMCIFParser(QUIET=True)
    asym_struct = parser.get_structure("1K3B", str(brutos / "1K3B.cif"))[0]
    asm_struct = parser.get_structure("1K3B_asm", str(brutos / "1K3B-assembly1.cif"))[0]

    # 1. Verifica mapeamento 1K3B (auth_seq_id + 24 = UniProt) e gera 1K3B_mature_chainA.pdb
    mapa = []
    pdb_lines = []
    atom_serial = 1
    idx0 = 0
    for ch_id in ("A", "B", "C"):
        ch = asym_struct[ch_id]
        for res in ch:
            if res.id[0] != " ":
                continue
            auth_num = res.id[1]
            up_pos = auth_num + 24
            aa1 = AA3TO1.get(res.resname, "X")
            assert seq[up_pos - 1] == aa1, f"Mismatch {ch_id} {auth_num} -> UP {up_pos}: {aa1} != {seq[up_pos - 1]}"
            mapa.append({
                "idx_0based": idx0,
                "uniprot_pos": up_pos,
                "wt_aa": aa1,
                "pdb_chain": ch_id,
                "pdb_auth_resnum": auth_num,
                "dominio": dominio_uniprot(up_pos),
            })
            idx0 += 1
            for atom in res:
                if atom.element == "H":
                    continue
                x, y, z = atom.get_coord()
                aname = atom.get_name()
                aname_fmt = f" {aname:<3s}" if len(aname) < 4 else aname[:4]
                pdb_lines.append(
                    f"ATOM  {atom_serial:5d} {aname_fmt} {res.resname:>3s} A{up_pos:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{atom.get_bfactor():6.2f}          {atom.element:>2s}  \n"
                )
                atom_serial += 1
    pdb_lines.append("END\n")
    (est_dir / "1K3B_mature_chainA.pdb").write_text("".join(pdb_lines), encoding="utf-8")
    with open(est_dir / "1K3B_residuo_mapa.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(mapa[0].keys()))
        w.writeheader()
        w.writerows(mapa)
    print(f"1K3B_mature_chainA.pdb: {len(mapa)} resíduos mapeados 100% com P53634")

    # 2. Analisa o tetrâmero biológico (1K3B-assembly1.cif)
    # Identifica monômeros: cada monômero tem cadeias equivalentes a A (1..119), B (207..370), C (371..439)
    print("Cadeias no assembly1:", [c.id for c in asm_struct])
    # Calcula SASA no tetrâmero completo
    sr = ShrakeRupley()
    sr.compute(asm_struct, level="R")

    # Agrupa resíduos por posição UniProt do primeiro monômero (cadeias A, B, C) e demais cadeias
    mono1_chains = {"A", "B", "C"}
    mono1_res = {}
    other_coords = []
    cl_coords = []
    cys258_coords = []
    active_coords = []

    for ch in asm_struct:
        for res in ch:
            if res.resname == "CL":
                cl_coords.extend(heavy_coords(res))
                continue
            if res.id[0] != " ":
                continue
            up_pos = res.id[1] + 24
            coords = heavy_coords(res)
            if ch.id in mono1_chains:
                mono1_res[up_pos] = (res, coords)
                if up_pos == 258:
                    cys258_coords.extend(coords)
                if up_pos in ACTIVE_SITE:
                    active_coords.extend(coords)
            else:
                other_coords.extend(coords)

    print(f"Assembly1: {len(mono1_res)} resíduos no monômero 1, {len(other_coords)} átomos nos outros monômeros, {len(cl_coords)} íons CL")

    rows = []
    for up_pos in range(1, 464):
        aa = seq[up_pos - 1]
        dom = dominio_uniprot(up_pos)
        info = mono1_res.get(up_pos)
        if info:
            res, coords = info
            d_cys = min_dist_coords(coords, cys258_coords)
            d_act = min_dist_coords(coords, active_coords)
            d_cl = min_dist_coords(coords, cl_coords)
            d_tet = min_dist_coords(coords, other_coords)
            rsa = res.sasa / MAX_ASA[aa] if aa in MAX_ASA else None
        else:
            d_cys = d_act = d_cl = d_tet = rsa = None

        rows.append({
            "uniprot_pos": up_pos,
            "wt_aa": aa,
            "dominio": dom,
            "na_1k3b": "sim" if info else "nao",
            "dist_cys258": round(d_cys, 2) if d_cys is not None else "",
            "dist_sitio_ativo": round(d_act, 2) if d_act is not None else "",
            "dist_cloreto": round(d_cl, 2) if d_cl is not None else "",
            "dist_interface_tetramero": round(d_tet, 2) if d_tet is not None else "",
            "rsa_tetramero": round(min(rsa, 1.0), 3) if rsa is not None else "",
            "dissulfeto": "sim" if up_pos in DISULFIDE else "nao",
            "glicosilacao": "sim" if up_pos in GLYCOSYLATION else "nao",
        })

    out_csv = base / "descritores_estruturais.csv"
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Salvo {out_csv} ({len(rows)} resíduos)")


if __name__ == "__main__":
    main()
