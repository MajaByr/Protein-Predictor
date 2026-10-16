# :dna: Protein Function Predictor (PFP) :dna:

*End-to-end framework that predicts protein function (in the form of Gene Ontology term probabilities) from an amino acid protein sequence.*

---

## How to run

### Requirements

- [uv](https://docs.astral.sh/uv/) (Python package manager)
- [Docker](https://www.docker.com/) (for MMseqs2 clustering)

All input data (GAF and UniProt TSV files) are already in the repository – no manual download needed.

### Steps

```bash
# 1. Clone the repository
git clone https://github.com/MajaByr/Protein-Predictor.git
cd Protein-Predictor

# 2. Install the Python environment (creates .venv from uv.lock)
uv sync

# 3. Merge GAF and UniProt data, keep only GO terms with high-quality evidence
uv run python -m projekt.merge

# 4. Cluster sequences with MMseqs2 (min-seq-id 0.5, coverage 0.8, cov-mode 0)
#    run step by step inside a Docker container – see docs/02_TRAIN_TEST.md

# 5. Assign whole clusters to train/eval/test at 9:1:1 ratio
uv run python -m projekt.train_test_split
```

### Output

| File | Created by | Content |
|------|-----------|---------|
| `data/processed/master.tsv` | step 3 | one protein per row: `entry`, `organism`, `go_ids` (labels), `sequence` |
| `data/processed/proteins.fasta` | step 3 | sequences of the proteins in `master.tsv` |
| `data/processed/go_comparison.tsv` | step 3 | per organism: GO terms kept vs dropped per protein |
| `data/clusters/clusters_cluster.tsv` | step 4 | MMseqs2 clusters: representative, member |
| `data/splits/{train,eval,test}.tsv` | step 5 | `master.tsv` rows split by cluster, plus `cluster` column |

Detailed description (in Polish): [docs/01_DATA.md](docs/01_DATA.md), [docs/02_TRAIN_TEST.md](docs/02_TRAIN_TEST.md).

## Data sources and versions

### GO annotations (GAF, UniProt-centric)

Source: https://geneontology.org/docs/download-go-annotations/downloads/

Stored gzip-compressed (`gzip -9`) to stay under GitHub's 100 MB file limit; pandas reads them directly.

| Organism  | File                     | Taxon ID | Generated        | GO release |
|-----------|--------------------------|----------|------------------|------------|
| Human     | `HUMAN-uniprot.gaf.gz`   | 9606     | 2026-07-28 19:42 | 2026-07-26 |
| Zebrafish | `DANRE-uniprot.gaf.gz`   | 7955     | 2026-07-28 19:38 | 2026-07-26 |
| Rat       | `RAT-uniprot.gaf.gz`     | 10116    | 2026-07-28 19:51 | 2026-07-26 |

### UniProt ID-to-sequence mappings

Source: https://www.uniprot.org/ — UniProtKB release **2026_03** (released 2026-09-02), accessed 2026-10-16.

Query: `(organism_id:<taxon>) AND (reviewed:true)` (Swiss-Prot only), columns: Entry, Gene Ontology IDs, Sequence, format TSV.

| Organism  | Taxon ID | Swiss-Prot entries |
|-----------|----------|--------------------|
| Human     | 9606     | 20431              |
| Zebrafish | 7955     | 3385               |
| Rat       | 10116    | 8233               |
