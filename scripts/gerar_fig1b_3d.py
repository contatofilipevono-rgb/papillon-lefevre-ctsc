"""Renderiza a Figura 1B em 3D (WebGL via 3Dmol.js + Microsoft Edge headless em 4400x2240 px):
  - Painel Esquerdo: Estrutura experimental em cristal da catepsina C madura (PDB 1K3B, 2,15 Å) com os sítios patogênicos
    coloridos pelo mecanismo (azul #2a78d6 = tríade catalítica <4 Å; violeta #6250d6 = ponte dissulfeto;
    laranja #eb6834 = desestabilizante ΔΔG >= 2,25 kcal/mol; cinza #8a8984 = outros / splicing A253T).
  - Painel Direito: Modelo completo do zimogênio no AlphaFold DB (AF-P53634-F1, resíduos 25-463) colorido por confiança
    pLDDT (azul escuro >90, azul claro 70-90, amarelo 50-70, laranja <50), destacando o propeptídeo (135-230) e suas
    mutações patogênicas (G139R, Y168C, L172P, L196P; AUC = 0,946).

Uso: py -I -X utf8 gerar_fig1b_3d.py <pasta_papillon_lefevre>
"""
import csv
import json
import os
import subprocess
import sys
from pathlib import Path
from PIL import Image

EDGE_PATHS = [
    r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
    r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
]


def find_edge():
    for p in EDGE_PATHS:
        if os.path.exists(p):
            return p
    raise RuntimeError("Microsoft Edge não encontrado")


def main():
    base = Path(sys.argv[1])
    pdb_1k3b = (base / "estruturas" / "1K3B_mature_chainA.pdb").read_text(encoding="utf-8")
    pdb_af = (base / "dados_brutos" / "AF-P53634-F1.pdb").read_text(encoding="utf-8")

    var_rows = list(csv.DictReader(open(base / "variantes_CTSC.csv", encoding="utf-8-sig")))
    P_1k = [r for r in var_rows if r["coorte_patogenica"] in ("calibracao", "teste_externo") and r["na_1k3b"] == "sim"]

    sitios = {}
    for r in P_1k:
        pos = int(r["posicao"])
        da = float(r["dist_sitio_ativo"])
        ddg = float(r["ddg_1k3b"])
        if da < 4.0:
            m = "ativo"
        elif r["dissulfeto"] == "sim":
            m = "dissulfeto"
        elif ddg >= 2.25:
            m = "desestab"
        else:
            m = "nenhum"
        # Prioriza ativo > dissulfeto > desestab > nenhum se houver mais de uma variante no mesmo resíduo
        prio = {"ativo": 4, "dissulfeto": 3, "desestab": 2, "nenhum": 1}
        if pos not in sitios or prio[m] > prio[sitios[pos]["mec"]]:
            sitios[pos] = {"pos": pos, "mec": m, "var": r["variante"], "ddg": ddg, "da": da}

    # Rótulos representativos no painel 1K3B (espalhados sem sobreposição)
    labels_1k3b = [
        {"pos": 258, "txt": "C258 (triad)", "cor": "#2a78d6", "dx": -95, "dy": -55},
        {"pos": 405, "txt": "H405R (triad)", "cor": "#2a78d6", "dx": -105, "dy": 20},
        {"pos": 252, "txt": "Q252L*", "cor": "#2a78d6", "dx": -90, "dy": -95},
        {"pos": 255, "txt": "C255R/S", "cor": "#6250d6", "dx": -110, "dy": 75},
        {"pos": 291, "txt": "C291Y", "cor": "#6250d6", "dx": 85, "dy": -85},
        {"pos": 337, "txt": "C337S", "cor": "#6250d6", "dx": 95, "dy": -25},
        {"pos": 272, "txt": "R272P/H", "cor": "#eb6834", "dx": 95, "dy": 45},
        {"pos": 301, "txt": "G301S", "cor": "#eb6834", "dx": -100, "dy": 115},
        {"pos": 347, "txt": "Y347C", "cor": "#eb6834", "dx": 95, "dy": 95},
        {"pos": 386, "txt": "G386R", "cor": "#eb6834", "dx": 80, "dy": 135},
        {"pos": 127, "txt": "H127P", "cor": "#eb6834", "dx": -110, "dy": -15},
        {"pos": 253, "txt": "A253T (splicing)", "cor": "#8a8984", "dx": 10, "dy": -115},
    ]

    # Rótulos no painel AlphaFold DB (zimogênio com propeptídeo 135-230)
    labels_af = [
        {"pos": 139, "txt": "G139R (ΔΔG 3.12)", "cor": "#eb6834", "dx": -115, "dy": -65},
        {"pos": 168, "txt": "Y168C (ΔΔG 2.35)", "cor": "#eb6834", "dx": 95, "dy": -75},
        {"pos": 172, "txt": "L172P (zymogen)", "cor": "#eb6834", "dx": 105, "dy": 10},
        {"pos": 196, "txt": "L196P (ΔΔG 2.07)", "cor": "#eb6834", "dx": 95, "dy": 85},
        {"pos": 258, "txt": "Cys258 (active site)", "cor": "#2a78d6", "dx": -115, "dy": 55},
    ]

    html_content = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ width: 2200px; height: 1120px; background: #ffffff; font-family: 'Segoe UI', Arial, sans-serif; color: #111; overflow: hidden; }}
  .header {{ height: 54px; padding: 10px 24px; display: flex; align-items: center; justify-content: space-between; border-bottom: 1px solid #e6e5e1; white-space: nowrap; }}
  .title {{ font-size: 17.5px; font-weight: 700; color: #0b0b0b; }}
  .legend {{ display: flex; gap: 16px; align-items: center; font-size: 13px; color: #333; }}
  .dot {{ width: 13px; height: 13px; border-radius: 50%; display: inline-block; margin-right: 6px; vertical-align: middle; border: 1px solid rgba(0,0,0,0.25); }}
  .rect {{ width: 16px; height: 11px; display: inline-block; margin-right: 6px; vertical-align: middle; border: 1px solid rgba(0,0,0,0.2); }}
  .grid {{ display: flex; width: 2200px; height: 1014px; padding: 14px 22px; gap: 22px; }}
  .card {{ flex: 1; border: 1px solid #e2e1dd; border-radius: 10px; position: relative; background: #ffffff; overflow: hidden; }}
  .card-hdr {{ position: absolute; top: 16px; left: 22px; z-index: 10; pointer-events: none; }}
  .card-ttl {{ font-size: 18px; font-weight: 700; color: #0b0b0b; }}
  .card-sub {{ font-size: 13.5px; color: #52514e; margin-top: 3px; }}
  .card-ftr {{ position: absolute; bottom: 14px; left: 22px; right: 22px; display: flex; justify-content: space-between; font-size: 13px; font-weight: 600; color: #52514e; z-index: 10; pointer-events: none; }}
  .viewer {{ width: 100%; height: 100%; }}
  .footer {{ height: 48px; padding: 0 28px; display: flex; align-items: center; justify-content: space-between; font-size: 12.5px; color: #52514e; border-top: 1px solid #e6e5e1; }}
</style>
</head>
<body>
  <div class="header">
    <div class="title">B. Three-dimensional mapping of pathogenic CTSC mechanisms (PDB 1K3B) and AlphaFold full-length zymogen (AF-P53634-F1)</div>
    <div class="legend">
      <span><i class="dot" style="background:#2a78d6"></i>Catalytic triad (&lt; 4.0 Å)</span>
      <span><i class="dot" style="background:#6250d6"></i>Disulfide bond</span>
      <span><i class="dot" style="background:#eb6834"></i>Destabilizing (ΔΔG ≥ 2.25 kcal/mol)</span>
      <span><i class="dot" style="background:#8a8984"></i>Other / splicing</span>
      <span><i class="rect" style="background:#0053d6"></i>AF pLDDT &gt; 90</span>
      <span><i class="rect" style="background:#65cbf3"></i>70–90</span>
      <span><i class="rect" style="background:#ffdb13"></i>&lt; 70</span>
    </div>
  </div>
  <div class="grid">
    <div class="card" id="card1">
      <div class="card-hdr">
        <div class="card-ttl">Experimental mature human cathepsin C crystal structure (PDB 1K3B, 2.15 Å)</div>
        <div class="card-sub">352 resolved residues (exclusion domain 25–143, heavy chain 231–394, light chain 395–463) | 51 pathogenic variants</div>
      </div>
      <div id="v1" class="viewer"></div>
      <div class="card-ftr">
        <span>Exclusion domain (β-barrel, aa 25–134)</span>
        <span>Catalytic heavy &amp; light chains (papain-like fold, aa 231–463)</span>
      </div>
    </div>
    <div class="card" id="card2">
      <div class="card-hdr">
        <div class="card-ttl">AlphaFold full-length pro-cathepsin C zymogen model (AF-P53634-F1, pLDDT confidence)</div>
        <div class="card-sub">Covers unobserved activation propeptide (aa 135–230, orange spheres; ThermoMPNN propeptide AUC = 0.946, p = 0.0031)</div>
      </div>
      <div id="v2" class="viewer"></div>
      <div class="card-ftr">
        <span>Mature enzyme scaffold (pLDDT &gt; 90, dark blue)</span>
        <span>Activation propeptide cap (aa 135–230) &amp; signal peptide (aa 1–24)</span>
      </div>
    </div>
  </div>
  <div class="footer">
    <span>*Q252L satisfies both catalytic triad contact (3.13 Å) and destabilization (ΔΔG = 2.21 kcal/mol). A253T (c.757G&gt;A) acts via exonic splicing (AlphaGenome v2 = 0.8589).</span>
    <span>Crystal coordinates from RCSB PDB 1K3B (Turk et al. 2001); zymogen from AlphaFold DB (AF-P53634-F1, non-commercial research use only).</span>
  </div>

<script>
const pdb1k3b = {json.dumps(pdb_1k3b)};
const pdbAf = {json.dumps(pdb_af)};
const sitios = {json.dumps(list(sitios.values()))};
const labels1k3b = {json.dumps(labels_1k3b)};
const labelsAf = {json.dumps(labels_af)};
const COR = {{ ativo: "#2a78d6", dissulfeto: "#6250d6", desestab: "#eb6834", nenhum: "#8a8984" }};

window.addEventListener("DOMContentLoaded", () => {{
  // Viewer 1: 1K3B
  const v1 = $3Dmol.createViewer("v1", {{ backgroundColor: "white", antialias: true }});
  v1.addModel(pdb1k3b, "pdb");
  v1.setStyle({{}}, {{ cartoon: {{ color: "#d5dde6", opacity: 0.88, thickness: 0.32 }} }});
  // Destaca domínio de exclusão (25-134) num tom levemente diferente do esqueleto
  v1.setStyle({{ resi: "25-134" }}, {{ cartoon: {{ color: "#c2cedb", opacity: 0.88, thickness: 0.32 }} }});

  for (const s of sitios) {{
    const c = COR[s.mec];
    v1.addStyle({{ resi: s.pos }}, {{ stick: {{ radius: 0.18, color: c }} }});
    v1.addStyle({{ resi: s.pos, atom: "CA" }}, {{ sphere: {{ radius: 1.35, color: c }} }});
  }}
  for (const lb of labels1k3b) {{
    v1.addLabel(lb.txt, {{
      fontSize: 12, fontColor: "#111111", backgroundColor: "#ffffff", backgroundOpacity: 0.88,
      borderColor: lb.cor, borderThickness: 1.5, inFront: true, showBackground: true
    }}, {{ resi: lb.pos, atom: "CA" }});
  }}
  v1.zoomTo();
  v1.rotate(18, "y");
  v1.rotate(-10, "x");
  v1.zoom(1.08);
  v1.render();

  // Viewer 2: AlphaFold DB AF-P53634-F1 colorido por pLDDT (b-factor)
  const v2 = $3Dmol.createViewer("v2", {{ backgroundColor: "white", antialias: true }});
  v2.addModel(pdbAf, "pdb");
  const plddtColor = function(atom) {{
    const b = atom.b;
    if (b >= 90) return "#0053d6";
    if (b >= 70) return "#65cbf3";
    if (b >= 50) return "#ffdb13";
    return "#ff7d45";
  }};
  // Oculta peptídeo sinal flexível 1-24 para focar no zimogênio 25-463
  v2.setStyle({{ resi: "1-24" }}, {{ cartoon: {{ color: "#ff7d45", opacity: 0.35, thickness: 0.2 }} }});
  v2.setStyle({{ resi: "25-463" }}, {{ cartoon: {{ colorfunc: plddtColor, opacity: 0.90, thickness: 0.32 }} }});
  for (const lb of labelsAf) {{
    v2.addStyle({{ resi: lb.pos }}, {{ stick: {{ radius: 0.22, color: lb.cor }} }});
    v2.addStyle({{ resi: lb.pos, atom: "CA" }}, {{ sphere: {{ radius: 1.55, color: lb.cor }} }});
    v2.addLabel(lb.txt, {{
      fontSize: 12.5, fontColor: "#111111", backgroundColor: "#ffffff", backgroundOpacity: 0.90,
      borderColor: lb.cor, borderThickness: 1.5, inFront: true, showBackground: true
    }}, {{ resi: lb.pos, atom: "CA" }});
  }}
  v2.zoomTo({{ resi: "25-463" }});
  v2.rotate(18, "y");
  v2.rotate(-10, "x");
  v2.zoom(1.08);
  v2.render();
}});
</script>
</body>
</html>"""

    tmp_html = base / "figuras" / "_tmp_fig1b_3d.html"
    out_png = base / "figuras" / "Fig1B_3D_1K3B_AF_mecanismos.png"
    out_pdf = base / "figuras" / "Fig1B_3D_1K3B_AF_mecanismos.pdf"
    sub_png = base / "submissao" / "figuras" / "Figure_1B_3D_1K3B_AF_mechanisms.png"
    sub_tif = base / "submissao" / "figuras" / "Figure_1B_3D_1K3B_AF_mechanisms.tif"
    sub_pdf = base / "submissao" / "figuras" / "Figure_1B_3D_1K3B_AF_mechanisms.pdf"

    tmp_html.write_text(html_content, encoding="utf-8")
    edge = find_edge()
    try:
        subprocess.run(
            [
                edge, "--headless", "--disable-gpu", "--hide-scrollbars",
                "--window-size=2200,1120", "--force-device-scale-factor=2",
                "--virtual-time-budget=5000",
                f"--screenshot={out_png}",
                tmp_html.resolve().as_uri(),
            ],
            check=True,
            capture_output=True,
            timeout=120,
        )
    finally:
        if tmp_html.exists():
            tmp_html.unlink()

    img = Image.open(out_png).convert("RGB")
    img.save(out_pdf, "PDF", resolution=300.0)
    img.save(sub_png, "PNG")
    img.save(sub_pdf, "PDF", resolution=300.0)
    img.save(sub_tif, compression="tiff_lzw", dpi=(600, 600))
    print(f"Fig1B 3D gerada com sucesso: {out_png} ({img.size[0]}x{img.size[1]} px)")


if __name__ == "__main__":
    main()
