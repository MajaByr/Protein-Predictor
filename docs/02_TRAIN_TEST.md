# Podział train/eval/test

Cel: podzielić białka z `master.tsv` na zbiory treningowy, walidacyjny i testowy w proporcji **9:1:1** tak, żeby podobne białka nie trafiły do różnych zbiorów.

## Dlaczego nie losowo?

Białka o podobnych sekwencjach zwykle pełnią podobne funkcje (np. to samo białko u człowieka i szczura). Przy losowym podziale w zbiorze testowym byłyby prawie identyczne kopie białek z treningu, wynik testu byłby zawyżony.

Dlatego najpierw **grupujemy podobne sekwencje w klastry**, a potem przypisujemy do zbiorów **całe klastry**.

## Krok 4.1: klastrowanie (mmseqs2 w Dockerze)

MMseqs2 uruchamiamy w Dockerze, w oficjalnym obrazie przypiętym do wersji `18-8cc5c` (powtarzalność wyników). Wymaga uruchomionego Docker Desktop.

**1. W PowerShellu, w katalogu repo** – uruchamiamy kontener z terminalem w środku:

```
docker run -it --rm -v "${PWD}/data:/data" --entrypoint sh ghcr.io/soedinglab/mmseqs2:18-8cc5c
```

`-v` podłącza folder `data/` z repo pod `/data` w kontenerze, `--rm` usuwa kontener po wyjściu. Pojawi się znak zachęty `#`.

**2. W kontenerze** – po kolei, czekając aż każda komenda się skończy:

```
mkdir -p /data/clusters /tmp/work
E=/usr/local/bin/entrypoint
$E createdb /data/processed/proteins.fasta /tmp/work/db
$E cluster /tmp/work/db /tmp/work/clu /tmp/work/tmp --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads 4
$E createtsv /tmp/work/db /tmp/work/db /tmp/work/clu /data/clusters/clusters_cluster.tsv
ls -l /data/clusters
exit
```

| Komenda | Co robi |
|---|---|
| `E=/usr/local/bin/entrypoint` | w tym obrazie program nie nazywa się `mmseqs` – `entrypoint` sam wybiera wersję pod procesor (avx2/sse41/sse2) |
| `createdb` | zamienia FASTA na wewnętrzną bazę mmseqs |
| `cluster` | właściwe klastrowanie |
| `createtsv` | zapisuje wynik jako tabelę do `data/clusters/` |

Parametry klastrowania (z zadania):

| Parametr | Znaczenie |
|---|---|
| `--min-seq-id 0.5` | białka w jednym klastrze mają co najmniej 50% identycznych aminokwasów |
| `-c 0.8` | dopasowanie musi pokrywać co najmniej 80% długości sekwencji |
| `--cov-mode 0` | pokrycie 80% liczone dla obu sekwencji (zapytania i celu) |
| `--threads 4` | liczba wątków – nie wpływa na wynik, tylko na obciążenie laptopa |

> **Dlaczego nie `mmseqs easy-cluster`?** `easy-cluster` to skrót na te same kroki plus zapis sekwencji reprezentantów (`result2flat`). Ten ostatni krok zawieszał się w Dockerze na Windowsie, a jego wynik nie jest nam potrzebny – dlatego uruchamiamy tylko potrzebne kroki. Wynik klastrowania jest identyczny.

Wynik: `data/clusters/clusters_cluster.tsv` – dwie kolumny bez nagłówka:

```
P04637   P04637      ← reprezentant klastra, członek
P04637   Q9Z0J8
...
```

Pierwsza kolumna (reprezentant) służy jako ID klastra.

## Krok 4.2: przypisanie klastrów do zbiorów

```
uv run python -m projekt.train_test_split
```

Plik [`src/projekt/train_test_split.py`](../src/projekt/train_test_split.py). Klastry mają bardzo różne rozmiary (od 1 do kilkuset białek), więc nie da się ich po prostu podzielić 9:1:1 według liczby klastrów. Stosujemy procedurę iteracyjną:

1. Tasujemy klastry z ustalonym seedem (`SEED = 42`).
2. Po kolei przypisujemy każdy klaster do zbioru o największym **względnym** braku: `(docelowa liczba białek − obecna) / docelowa liczba białek`.
3. Dzięki temu wszystkie zbiory zapełniają się w tym samym tempie, a proporcje na końcu wychodzą prawie dokładnie 9:1:1.

Losowość jest więc zachowana, tylko na poziomie **klastrów**, a nie pojedynczych białek – każdy zbiór dostaje losową mieszankę rodzin białek, ale jedna rodzina nigdy nie jest rozdzielona (tzw. *grouped split*).

> **Dlaczego nie sortujemy klastrów od największych?** Przy sortowaniu i bezwzględnym braku train (który ma 9× większy cel) zgarniał wszystkie duże klastry, a do eval/test trafiały prawie same pojedyncze białka. Zbiory miałyby wtedy inny charakter biologicznie (duże rodziny białek tylko w treningu) i wynik testu byłby zafałszowany.

Proporcje liczymy w **białkach**, a nie w klastrach.

Na koniec skrypt sprawdza, że:
- żaden klaster nie jest rozdzielony między zbiory,
- każde białko jest w dokładnie jednym zbiorze,
- klastry pasują do aktualnego `master.tsv` (jeśli nie – trzeba ponownie uruchomić mmseqs).

## Wynik

`data/splits/train.tsv`, `eval.tsv`, `test.tsv` – kolumny jak w `master.tsv` (`entry`, `organism`, `go_ids`, `sequence`) plus `cluster` (ID klastra).

Skrypt wypisuje podsumowanie: liczbę białek i klastrów w każdym zbiorze, ich udział (oczekiwane ok. 0.818 / 0.091 / 0.091) oraz rozkład organizmów, żeby sprawdzić, czy podział jest zrównoważony.
