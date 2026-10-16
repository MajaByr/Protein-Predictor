# :dna: Protein Function Predictor (PFP) :dna:

*End-to-end framework that predicts protein function (in the form of Gene Ontology term probabilities) from an amino acid protein sequence.*

---


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
