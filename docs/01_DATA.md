# Dane

Cel: na podstawie **sekwencji** białka przewidzieć jego **funkcje** (zbiór terminów GO). To klasyfikacja wieloetykietowa, bo jedno białko ma wiele funkcji naraz.

## W skrócie: co zawiera GAF, a co TSV

- **GAF**: `entry` (unikalne ID białka), funkcja (`go_id`) i dowód na nią (`evidence`). Jeden wiersz to jedna para białko–funkcja, więc to samo `entry` powtarza się w wielu wierszach.
- **TSV (UniProt)**: również `entry` i funkcje, ale tym razem w jednym wierszu jest cały zbiór funkcji białka, a do tego jego sekwencja. Jedno `entry` to jeden wiersz.

GAF ma dowody (evidence), a TSV ich nie ma. Właśnie dlatego GAF jest nam potrzebny do filtrowania funkcji.
TSV ma sekwencje, a GAF nie. Dlatego TSV jest potrzebny jako wejście modelu.

## 1. TSV z UniProt: sekwencje i etykiety (`data/TSV/`)

Tabela pobrana z [uniprot.org](https://www.uniprot.org/) (tylko Swiss-Prot, czyli białka zrecenzowane ręcznie). Kolumny są rozdzielone tabulatorem, jeden wiersz to jedno białko.

```
Entry    Gene Ontology IDs                     Sequence
P04637   GO:0000122; GO:0000785; GO:0003677…   MEEPQSDPSVEPPLSQETFSDLWKLL…
```

- **Entry**: identyfikator białka w UniProt. Jest kluczem łączenia z plikami GAF.
- **Gene Ontology IDs**: wszystkie funkcje białka według UniProt. Z nich wybieramy nasze **etykiety** (tylko te z mocnym dowodem w GAF).
- **Sequence**: sekwencja aminokwasów (jedna litera to jeden aminokwas). To jest **wejście** modelu.

Jedno ID ma jedną sekwencję i wiele terminów GO.

## 2. GAF: adnotacje GO z dowodami (`data/GAF/`)

Pliki GAF (Gene Association File) pochodzą z [geneontology.org](https://geneontology.org/docs/download-go-annotations/downloads/). Jeden wiersz to jedna adnotacja („białko X pełni funkcję Y i dowodem na to jest Z”), więc jedno białko ma wiele wierszy. Ta sama funkcja danego białka może mieć wiele dowodów, np. z różnych eksperymentów albo publikacji, i każdy z nich jest osobnym wierszem. Linie zaczynające się od `!` to nagłówek z metadanymi.

> Pliki są spakowane gzipem (`.gaf.gz`), bo surowe GAF-y przekraczają limit 100 MB na GitHubie. Pandas czyta je bez rozpakowywania.

Najważniejsze kolumny:

| Kolumna     | Znaczenie                                                                                   |
|-------------|---------------------------------------------------------------------------------------------|
| `db`        | baza źródłowa (zostawiamy tylko `UniProtKB`)                                                |
| `entry`     | ID białka w UniProt                                                                         |
| `qualifier` | relacja białka do terminu, np. `enables`; `NOT\|...` oznacza, że białko tej funkcji **nie** pełni |
| `go_id`     | termin GO, np. `GO:0003677` (DNA binding)                                                   |
| `evidence`  | kod dowodowy, czyli skąd wiemy o tej funkcji (np. `IDA` to eksperyment, `IEA` to przewidywanie automatyczne) |
| `aspect`    | `F` (Molecular Function), `P` (Biological Process), `C` (Cellular Component)                |

GAF służy jako **filtr funkcji**: z funkcji białka wypisanych w UniProt zostawiamy tylko te, które w GAF mają co najmniej jeden dowód wysokiej jakości (`EXP`, `IDA`, `IPI`, `IMP`, `IGI`, `IEP`, `HTP`, `HDA`, `HMP`, `HGI`, `HEP`, `TAS`, `IC`). Funkcje potwierdzone tylko słabymi dowodami (np. `IEA` – automatyczne przewidywanie, `ISS` – homologia) odrzucamy. Filtr działa **per funkcja, a nie per białko** – tak samo definiuje się ground truth w benchmarku CAFA.

## 3. Co robimy i co chcemy dostać

```
GAF ──(które funkcje mają mocny dowód?)──┐
                                         ├──> master.tsv + proteins.fasta ──> podział train/eval/test (docs/02_TRAIN_TEST.md)
TSV ──(sekwencja + funkcje z UniProt)────┘
```

**Krok 3: łączenie danych** (`uv run python -m projekt.merge`, plik [`src/projekt/merge.py`](../src/projekt/merge.py))

1. Wczytujemy GAF i zostawiamy tylko adnotacje białek z UniProtKB, bez `NOT`, z mocnym dowodem.
2. Wczytujemy TSV z UniProt i zostawiamy tylko białka z co najmniej jednym terminem GO.
3. Łączymy obie tabele po `entry` (inner join) i dodajemy kolumnę z organizmem (Human / Zebrafish / Rat).
4. Dla każdego białka dzielimy funkcje z UniProt na:
   - **zachowane** – mają mocny dowód w GAF (to są etykiety),
   - **odrzucone** – nie mają mocnego dowodu.
5. Białka, którym nie została żadna zachowana funkcja, wypadają.
6. Zapisujemy wyniki do `data/processed/`:
   - `master.tsv`: jedno białko na wiersz, kolumny `entry`, `organism`, `go_ids` (etykiety – tylko zachowane funkcje), `sequence`;
   - `proteins.fasta`: same sekwencje w formacie FASTA (`>entry organism`, a pod spodem sekwencja);
   - `go_comparison.tsv`: per organizm – średnia liczba funkcji zachowanych vs odrzuconych na białko oraz ile białek wypadło.
