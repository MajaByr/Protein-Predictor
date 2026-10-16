from pathlib import Path

import pandas as pd

DATA = Path("data")
OUT = DATA / "processed"

ORGANISMS = {
    "Human": ("HUMAN-uniprot.gaf.gz", "uniprotkb_organism_id_9606_AND_reviewed_2026_10_16.tsv"),
    "Zebrafish": ("DANRE-uniprot.gaf.gz", "uniprotkb_organism_id_7955_AND_reviewed_2026_10_16.tsv"),
    "Rat": ("RAT-uniprot.gaf.gz", "uniprotkb_organism_id_10116_AND_reviewe_2026_10_16.tsv"),
}

#potwierdzenia wysokiej jakości 
HQ_EVIDENCE = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP", "TAS", "IC"}

#kolumny z pliku GAF, głównie potrzebne to entry, go_id, evidence, bo chcemy tylko wyfiltrować złe dowody 
GAF_COLUMNS = [
    "db", "entry", "symbol", "qualifier", "go_id", "reference", "evidence", "with_from",
    "aspect", "name", "synonym", "type", "taxon", "date", "assigned_by", "extension", "gene_product_form",
]

GAF_USEFUL_COLUMNS = ["db", "entry", "qualifier", "go_id", "evidence"]


def read_gaf(path: Path) -> pd.DataFrame:
    #tylko użyteczne kolumny     
    gaf = pd.read_csv(path, sep="\t", comment="!", header=None, names=GAF_COLUMNS,
                      usecols=GAF_USEFUL_COLUMNS, dtype=str)

    return gaf[(gaf["db"] == "UniProtKB") #tylko baza UniProtKB
               & ~gaf["qualifier"].str.startswith("NOT") #tylko te, które mają funkcję 
               & gaf["evidence"].isin(HQ_EVIDENCE)] #tylko mocne dowody 


def read_uniprot(path: Path) -> pd.DataFrame:
    
    up = pd.read_csv(path, sep="\t", dtype=str).rename(
        columns={"Entry": "entry", "Gene Ontology IDs": "go_ids", "Sequence": "sequence"})

    #zanieniamy wiele wartości GO na jeden zbiór, puste stringi to pusty zbiór
    up["go_uniprot"] = up["go_ids"].fillna("").map(lambda s: {g.strip() for g in s.split(";") if g.strip()})
    # usuwamy białka bez funkcji GO
    up = up[up["go_uniprot"].map(len) > 0]
    return up[["entry", "go_uniprot", "sequence"]]


def merge_organism(organism: str, gaf_file: str, tsv_file: str) -> pd.DataFrame:

    gaf = read_gaf(DATA / "GAF" / gaf_file)
    up = read_uniprot(DATA / "TSV" / tsv_file)

    # GAF ma jeden wiersz na (białko, funkcja, dowód) – grupujemy po białku i zbieramy
    # jego mocno potwierdzone funkcje w zbiór (duplikaty z różnych dowodów znikają same)
    go_hq = gaf.groupby("entry")["go_id"].agg(set).rename("go_gaf_hq")

    # Inner join po "entry": zostają tylko białka obecne w obu tabelach, czyli
    # białka, które mają mocny dowód i mają jakiekolwiek GO w UniProt
    df = up.join(go_hq, on="entry", how="inner")
    #dodanie kolumny z nazwą organizmu 
    df.insert(1, "organism", organism)
    
    return df


def merge_all() -> pd.DataFrame:
    return pd.concat([merge_organism(organism, gaf_file, tsv_file)
                      for organism, (gaf_file, tsv_file) in ORGANISMS.items()], ignore_index=True)


def split_functions(df: pd.DataFrame) -> pd.DataFrame:
    #w dataseciepo merge_organism zostały dwie kolumny z go_ids: go_gaf_hq i go_uniprot
    #chcemy zostawić jedynie funkcje obecne w obu źródłach.
    # go_confirmed – funkcje z UniProt, które mają mocny dowód w GAF (zachowujemy je jako etykiety)
    df["go_confirmed"] = [u & g for u, g in zip(df["go_uniprot"], df["go_gaf_hq"])]
    # go_dropped – funkcje z UniProt bez mocnego dowodu w GAF (odrzucamy je)
    df["go_dropped"] = [u - g for u, g in zip(df["go_uniprot"], df["go_gaf_hq"])]
    return df


def keep_confirmed_functions(df: pd.DataFrame) -> pd.DataFrame:
    # Białka, którym nie została żadna potwierdzona funkcja, wypadają ze zbioru
    return df[df["go_confirmed"].map(len) > 0].copy()


def compare_go(df: pd.DataFrame) -> pd.DataFrame:
    # Porównanie per organizm: ile funkcji z UniProt zachowaliśmy, a ile odrzuciliśmy
    # (liczone na wszystkich białkach przed odfiltrowaniem tych bez potwierdzonych funkcji)
    kept = df["go_confirmed"].map(len)
    dropped = df["go_dropped"].map(len)
    rows = []
    for organism, sub in df.assign(kept=kept, dropped=dropped).groupby("organism"):
        rows.append({
            "organism": organism,
            "proteins": len(sub),
            # średnia liczba funkcji na białko: wszystkie w UniProt / zachowane / odrzucone
            "mean_uniprot_terms": round(sub["go_uniprot"].map(len).mean(), 2),
            "mean_kept_terms": round(sub["kept"].mean(), 2),
            "mean_dropped_terms": round(sub["dropped"].mean(), 2),
            # łączna liczba par (białko, funkcja) zachowanych / odrzuconych
            "total_kept_terms": int(sub["kept"].sum()),
            "total_dropped_terms": int(sub["dropped"].sum()),
            # białka, którym zachowaliśmy wszystkie funkcje / żadnej (te wypadają z datasetu)
            "proteins_all_kept": int((sub["dropped"] == 0).sum()),
            "proteins_none_kept": int((sub["kept"] == 0).sum()),
        })
    return pd.DataFrame(rows)


def build_master(df: pd.DataFrame) -> pd.DataFrame:
    master = df[["entry", "organism", "go_confirmed", "sequence"]].copy()
    master["go_confirmed"] = master["go_confirmed"].map(lambda s: ";".join(sorted(s)))
    # Sortujemy po ID białka, żeby kolejność wierszy nie zależała od kolejności organizmów
    return master.rename(columns={"go_confirmed": "go_ids"}).sort_values("entry")


def write_fasta(master: pd.DataFrame, path: Path) -> None:
    with open(path, "w") as f:
        for entry, organism, seq in zip(master["entry"], master["organism"], master["sequence"]):
            f.write(f">{entry} {organism}\n{seq}\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    df = split_functions(merge_all())

    comparison = compare_go(df)
    comparison.to_csv(OUT / "go_comparison.tsv", sep="\t", index=False)

    df = keep_confirmed_functions(df)

    master = build_master(df)
    master.to_csv(OUT / "master.tsv", sep="\t", index=False)
    write_fasta(master, OUT / "proteins.fasta")
    
    print("Processed.")


if __name__ == "__main__":
    main()
