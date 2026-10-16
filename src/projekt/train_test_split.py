from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("data")
MASTER = DATA / "processed" / "master.tsv"
CLUSTERS = DATA / "clusters" / "clusters_cluster.tsv"
OUT = DATA / "splits"

SEED = 42
RATIOS = {"train": 9, "eval": 1, "test": 1}

# IMPORTANT!!
# first run mmseqs clustering (see docs/02_TRAIN_TEST.md) to get clusters_cluster.tsv, then run this script

def read_clusters(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, names=["cluster", "entry"], dtype=str)


def assign_splits(clusters: pd.DataFrame, seed: int) -> pd.DataFrame:

    # Liczba białek w każdym klastrze; sort_index daje stałą kolejność niezależną od kolejności w pliku
    sizes = clusters["cluster"].value_counts().sort_index()

    #dążymy do tego by białka były podzielone w proporcji 9:1:1, a zbiory zawierały całe klastry.
    #klastry bierzemy w losowej kolejności (bez sortowania po rozmiarze)
    #i przydzielamy je do zbioru, w którym najbardziej brakuje
    #ale robimy to względem wielkości celu, by train nie dostawał największych klastrów

    rng = np.random.default_rng(seed)
    sizes = sizes.iloc[rng.permutation(len(sizes))]

    total = sizes.sum()
    target = {s: total * r / sum(RATIOS.values()) for s, r in RATIOS.items()}
    current = dict.fromkeys(RATIOS, 0)
    assignment = {}

    for cluster, size in sizes.items():
        # względny brak: (cel - obecnie) / cel
        split = max(RATIOS, key=lambda s: (target[s] - current[s]) / target[s])
        assignment[cluster] = split
        current[split] += size

    return clusters.assign(split=clusters["cluster"].map(assignment))


def check_splits(df: pd.DataFrame) -> None: #sprawdzanie zasad
    assert (df.groupby("cluster")["split"].nunique() == 1).all(), "cluster split across sets"
    assert df["entry"].is_unique, "duplicated entries"


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    summary = df.groupby("split").agg(proteins=("entry", "size"), clusters=("cluster", "nunique"))
    summary["fraction"] = (summary["proteins"] / len(df)).round(3)
    # Rozkład organizmów w każdym zbiorze – sprawdzamy, czy podział jest zrównoważony
    by_organism = pd.crosstab(df["split"], df["organism"])
    return summary.join(by_organism).loc[list(RATIOS)]


def write_splits(df: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for split, sub in df.groupby("split"):
        sub.drop(columns="split").sort_values("entry").to_csv(OUT / f"{split}.tsv", sep="\t", index=False)


def main() -> None:
    master = pd.read_csv(MASTER, sep="\t", dtype=str)
    clusters = read_clusters(CLUSTERS)
    assert set(clusters["entry"]) == set(master["entry"]), "clusters do not match master.tsv, rerun mmseqs"

    df = master.merge(assign_splits(clusters, SEED), on="entry")
    check_splits(df)
    write_splits(df)

    print(summarize(df).to_string())
    print(f"\nSaved {', '.join(f'{s}.tsv' for s in RATIOS)} to {OUT}/")


if __name__ == "__main__":
    main()
