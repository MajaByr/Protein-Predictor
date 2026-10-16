"""Step 3: merge GAF annotations and UniProt (Swiss-Prot) tables into one master file.

Usage: uv run python -m projekt.merge
"""

from pathlib import Path

import pandas as pd

DATA = Path("data")
OUT = DATA / "processed"

ORGANISMS = {
    # name: (GAF file, UniProt TSV file)
    "Human": ("HUMAN-uniprot.gaf.gz", "uniprotkb_organism_id_9606_AND_reviewed_2026_10_16.tsv"),
    "Zebrafish": ("DANRE-uniprot.gaf.gz", "uniprotkb_organism_id_7955_AND_reviewed_2026_10_16.tsv"),
    "Rat": ("RAT-uniprot.gaf.gz", "uniprotkb_organism_id_10116_AND_reviewe_2026_10_16.tsv"),
}

HQ_EVIDENCE = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP", "TAS", "IC"}

GAF_COLUMNS = [
    "db", "entry", "symbol", "qualifier", "go_id", "reference", "evidence", "with_from",
    "aspect", "name", "synonym", "type", "taxon", "date", "assigned_by", "extension", "gene_product_form",
]


def read_gaf(path: Path) -> pd.DataFrame:
    gaf = pd.read_csv(path, sep="\t", comment="!", header=None, names=GAF_COLUMNS,
                      usecols=["db", "entry", "qualifier", "go_id", "evidence"], dtype=str)
    # drop non-UniProt rows (e.g. ComplexPortal) and negative annotations ("NOT|enables")
    return gaf[(gaf["db"] == "UniProtKB") & ~gaf["qualifier"].str.startswith("NOT")]


def read_uniprot(path: Path) -> pd.DataFrame:
    up = pd.read_csv(path, sep="\t", dtype=str).rename(
        columns={"Entry": "entry", "Gene Ontology IDs": "go_ids", "Sequence": "sequence"})
    up["go_uniprot"] = up["go_ids"].fillna("").map(lambda s: {g.strip() for g in s.split(";") if g.strip()})
    return up[["entry", "go_uniprot", "sequence"]]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for organism, (gaf_file, tsv_file) in ORGANISMS.items():
        gaf = read_gaf(DATA / "GAF" / gaf_file)
        up = read_uniprot(DATA / "TSV" / tsv_file)

        go_all = gaf.groupby("entry")["go_id"].agg(set).rename("go_gaf_all")
        hq = gaf[gaf["evidence"].isin(HQ_EVIDENCE)]
        go_hq = hq.groupby("entry")["go_id"].agg(set).rename("go_gaf_hq")

        df = up.join(go_all, on="entry").join(go_hq, on="entry")
        df.insert(1, "organism", organism)
        print(f"{organism}: {len(up)} Swiss-Prot, {gaf['entry'].nunique()} in GAF, "
              f"{df['go_gaf_hq'].notna().sum()} Swiss-Prot with high-quality evidence")
        frames.append(df)

    df = pd.concat(frames, ignore_index=True)
    for col in ["go_gaf_all", "go_gaf_hq"]:
        df[col] = df[col].map(lambda s: s if isinstance(s, set) else set())

    # keep well-annotated proteins: at least one high-quality GAF annotation and UniProt GO terms
    df = df[(df["go_gaf_hq"].map(len) > 0) & (df["go_uniprot"].map(len) > 0)].copy()

    # compare GAF vs UniProt GO IDs per protein
    rows = []
    for gaf_col in ["go_gaf_all", "go_gaf_hq"]:
        only_gaf = [g - u for g, u in zip(df[gaf_col], df["go_uniprot"])]
        only_up = [u - g for g, u in zip(df[gaf_col], df["go_uniprot"])]
        for organism, sub in df.assign(only_gaf=only_gaf, only_up=only_up).groupby("organism"):
            rows.append({
                "gaf_set": gaf_col,
                "organism": organism,
                "proteins": len(sub),
                "identical_sets": int(((sub["only_gaf"].map(len) == 0) & (sub["only_up"].map(len) == 0)).sum()),
                "proteins_with_gaf_only_terms": int((sub["only_gaf"].map(len) > 0).sum()),
                "proteins_with_uniprot_only_terms": int((sub["only_up"].map(len) > 0).sum()),
                "mean_gaf_terms": round(sub[gaf_col].map(len).mean(), 1),
                "mean_uniprot_terms": round(sub["go_uniprot"].map(len).mean(), 1),
                "mean_gaf_only_terms": round(sub["only_gaf"].map(len).mean(), 2),
                "mean_uniprot_only_terms": round(sub["only_up"].map(len).mean(), 2),
            })
    comparison = pd.DataFrame(rows)
    comparison.to_csv(OUT / "go_comparison.tsv", sep="\t", index=False)
    print("\nGO comparison (GAF vs UniProt):")
    print(comparison.to_string(index=False))

    # UniProt GO IDs are the ground truth labels
    master = df[["entry", "organism", "go_uniprot", "go_gaf_hq", "sequence"]].copy()
    for col in ["go_uniprot", "go_gaf_hq"]:
        master[col] = master[col].map(lambda s: ";".join(sorted(s)))
    master = master.rename(columns={"go_uniprot": "go_ids"}).sort_values("entry")
    master.to_csv(OUT / "master.tsv", sep="\t", index=False)

    with open(OUT / "proteins.fasta", "w") as f:
        for entry, organism, seq in zip(master["entry"], master["organism"], master["sequence"]):
            f.write(f">{entry} {organism}\n{seq}\n")

    print(f"\nSaved {len(master)} proteins to {OUT}/master.tsv and {OUT}/proteins.fasta")
    print(master["organism"].value_counts().to_string())


if __name__ == "__main__":
    main()
