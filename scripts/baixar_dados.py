"""Baixa dados públicos de variantes e estruturas do CTSC (UniProt P53634).

Fontes:
  - ClinVar (NCBI E-utilities)
  - gnomAD v4 (API GraphQL, incluindo homozygote_count)
  - AlphaMissense e estrutura AlphaFold DB (P53634)
  - UniProt (P53634 JSON + FASTA)
  - PDB 1K3B (mmCIF + PDB)

Uso: py -I -X utf8 baixar_dados.py <pasta_saida>
"""
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

GENE = "CTSC"
UNIPROT = "P53634"
PDB_ID = "1K3B"
UA = {"User-Agent": "academic-variant-curation/0.1"}


def get(url, data=None, headers=None, tries=4):
    h = dict(UA)
    if headers:
        h.update(headers)
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def clinvar(gene, out):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    term = f"{gene}[gene] AND single_gene[prop]"
    q = urllib.parse.urlencode({"db": "clinvar", "term": term, "retmax": 5000, "retmode": "json"})
    ids = json.loads(get(base + "esearch.fcgi?" + q))["esearchresult"]["idlist"]
    rows = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        q = urllib.parse.urlencode({"db": "clinvar", "id": ",".join(chunk), "retmode": "json"})
        res = json.loads(get(base + "esummary.fcgi?" + q))["result"]
        for uid in res.get("uids", []):
            d = res[uid]
            germ = d.get("germline_classification") or d.get("clinical_significance") or {}
            traits = germ.get("trait_set") or d.get("trait_set") or []
            rows.append({
                "gene": gene,
                "clinvar_uid": uid,
                "accession": d.get("accession", ""),
                "title": d.get("title", ""),
                "protein_change": d.get("protein_change", ""),
                "molecular_consequence": ";".join(d.get("molecular_consequence_list", []) or []),
                "variant_type": d.get("obj_type", ""),
                "classification": germ.get("description", ""),
                "review_status": germ.get("review_status", ""),
                "conditions": ";".join(t.get("trait_name", "") for t in traits),
            })
        time.sleep(0.4)
    with open(out / f"clinvar_{gene}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["gene"])
        w.writeheader()
        w.writerows(rows)
    return len(rows)


GNOMAD_Q = """
query($g: String!) {
  meta {
    clinvar_release_date
  }
  gene(gene_symbol: $g, reference_genome: GRCh38) {
    canonical_transcript_id
    gnomad_constraint {
      pli
      lof_z
      mis_z
      oe_lof
      oe_mis
    }
    variants(dataset: gnomad_r4) {
      variant_id
      consequence
      hgvsp
      hgvsc
      transcript_id
      exome { ac an homozygote_count }
      genome { ac an homozygote_count }
    }
  }
}"""


def gnomad(gene, out):
    body = json.dumps({"query": GNOMAD_Q, "variables": {"g": gene}}).encode()
    res = json.loads(get("https://gnomad.broadinstitute.org/api", data=body,
                         headers={"Content-Type": "application/json"}))
    if "errors" in res:
        raise RuntimeError(f"gnomAD GraphQL error: {res['errors']}")
    g = res["data"]["gene"]
    meta_info = {
        "canonical_transcript_id": g.get("canonical_transcript_id"),
        "gnomad_constraint": g.get("gnomad_constraint"),
        "meta": res["data"].get("meta"),
        "dataset": "gnomad_r4",
    }
    (out / f"gnomad_{gene}_meta.json").write_text(json.dumps(meta_info, indent=2), encoding="utf-8")

    rows = []
    for v in g["variants"]:
        ex = v.get("exome") or {}
        ge = v.get("genome") or {}
        ac = (ex.get("ac") or 0) + (ge.get("ac") or 0)
        an = (ex.get("an") or 0) + (ge.get("an") or 0)
        hom = (ex.get("homozygote_count") or 0) + (ge.get("homozygote_count") or 0)
        rows.append({
            "gene": gene,
            "variant_id": v["variant_id"],
            "consequence": v["consequence"],
            "hgvsp": v.get("hgvsp") or "",
            "hgvsc": v.get("hgvsc") or "",
            "transcript": v.get("transcript_id") or "",
            "ac": ac,
            "an": an,
            "af": (ac / an) if an else "",
            "homozygote_count": hom,
        })
    with open(out / f"gnomad_{gene}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return len(rows), g["canonical_transcript_id"]


def alphafold_and_am(gene, uniprot, out):
    meta = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"))
    entry = meta[0]
    (out / f"alphafold_{gene}_meta.json").write_text(json.dumps(entry, indent=2), encoding="utf-8")
    am_url = entry.get("amAnnotationsUrl")
    am_count = 0
    if am_url:
        data = get(am_url)
        (out / f"alphamissense_{gene}.csv").write_bytes(data)
        am_count = data.count(b"\n") - 1
    pdb_url = entry.get("pdbUrl")
    if pdb_url:
        (out / f"AF-{uniprot}-F1.pdb").write_bytes(get(pdb_url))
    return am_count


def uniprot_and_pdb(uniprot, pdb_id, out):
    up_json = get(f"https://rest.uniprot.org/uniprotkb/{uniprot}.json")
    (out / f"uniprot_{uniprot}.json").write_bytes(up_json)
    up_fasta = get(f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta")
    (out / f"uniprot_{uniprot}.fasta").write_bytes(up_fasta)
    (out / f"{pdb_id}.cif").write_bytes(get(f"https://files.rcsb.org/download/{pdb_id}.cif"))
    (out / f"{pdb_id}.pdb").write_bytes(get(f"https://files.rcsb.org/download/{pdb_id}.pdb"))


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("dados_brutos")
    out.mkdir(parents=True, exist_ok=True)
    print(GENE, "ClinVar:", clinvar(GENE, out))
    print(GENE, "gnomAD:", gnomad(GENE, out))
    print(GENE, "AlphaMissense:", alphafold_and_am(GENE, UNIPROT, out))
    uniprot_and_pdb(UNIPROT, PDB_ID, out)
    print(GENE, "UniProt + PDB 1K3B: OK")


if __name__ == "__main__":
    main()
