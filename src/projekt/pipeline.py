import skbio
import torch
import torch.nn as nn
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from pathlib import Path

# ==========================================
# 1. KONFIGURACJA I SŁOWNIK
# ==========================================

torch.manual_seed(42)


# 20 standardowych aminokwasów + X (nieznany) + PAD (wypełnienie)
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY") + ["X", "PAD"]
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
VOCAB_SIZE = len(AMINO_ACIDS)  # N = 22
PAD_IDX = AA_TO_IDX["PAD"]

# Przykładowe (zahardcodowane) terminy GO z wymaganego Spotkania 2
# W docelowym rozwiązaniu wczytasz je z prawdziwych list GO.
GO_VOCABULARY = {
    "BP": ["GO:0008150", "GO:0009987", "GO:0008152"],  # Biological Process
    "MF": ["GO:0003674", "GO:0005488"],                # Molecular Function
    "CC": ["GO:0005575", "GO:0005623", "GO:0005886", "GO:0043226"]  # Cellular Component
}

# ==========================================
# 2. DATASET (Ładowanie z użyciem scikit-bio)
# ==========================================

class ProteinDataset(Dataset):
    def __init__(self, fasta_file: str, max_len: int = 1024):
        """
        Wczytuje plik FASTA, koduje sekwencje i przygotowuje maski.
        """
        self.max_len = max_len
        self.proteins = []
        
        # Użycie scikit-bio do parsowania FASTA zgodnie z poleceniem
        for seq in skbio.io.read(fasta_file, format='fasta'):
            self.proteins.append({
                'id': seq.metadata['id'],
                'sequence': str(seq)
            })

    def __len__(self):
        return len(self.proteins)

    def __getitem__(self, idx):
        prot = self.proteins[idx]
        seq_str = prot['sequence'].upper()  # Upewniamy się, że wszystkie aminokwasy są wielkimi literami
        
        # Przycinanie długich sekwencji (Truncation)
        seq_str = seq_str[:self.max_len]
        
        # Inicjalizacja tensora dla sekwencji o kształcie (N, D)
        # Gdzie N to VOCAB_SIZE, a D to max_len
        seq_tensor = torch.zeros((VOCAB_SIZE, self.max_len))
        mask = torch.zeros(self.max_len, dtype=torch.bool)
        
        # Kodowanie One-Hot i ustawianie maski na True dla ważnych aminokwasów
        for i, aa in enumerate(seq_str):
            aa_idx = AA_TO_IDX.get(aa, AA_TO_IDX["X"]) # X jeśli nieznany
            seq_tensor[aa_idx, i] = 1.0
            mask[i] = True 
            
        # Padding (wypełnianie) dla krótszych sekwencji
        for i in range(len(seq_str), self.max_len):
            seq_tensor[PAD_IDX, i] = 1.0
            mask[i] = False # False oznacza, że to tylko wypełniacz (padding)
            
        return prot['id'], seq_tensor, mask

# ==========================================
# 3. DUMMY MODEL (Atrapa architektury)
# ==========================================

class DummyPredictor(nn.Module):
    def __init__(self, vocab_size: int, max_len: int, num_go_terms: int):
        super().__init__()
        self.num_go_terms = num_go_terms
        
        # Puste warstwy - model atrapa. 
        # Zapisujemy parametry, by interfejs zgadzał się z docelowym modelem.
        self.vocab_size = vocab_size
        self.max_len = max_len

    def forward(self, x, mask=None):
        """
        Oczekiwany input: 
        x kształt: (B, N, D) -> Batch Size, Vocab Size, Max Length
        mask kształt: (B, D)
        
        Zgodnie z poleceniem, interfejs przyjmuje maskę, ale dummy model 
        zwraca po prostu losowe wyniki w kształcie (B, P).
        """
        batch_size = x.size(0)
        
        # Generowanie losowych przewidywań z zakresu 0.0 - 1.0
        # Kształt wyjściowy: (B, P), gdzie P to liczba terminów GO dla danej ontologii
        random_scores = torch.rand((batch_size, self.num_go_terms), device=x.device)
        
        return random_scores

# ==========================================
# 4. GŁÓWNY RUROCIĄG (End-to-End Pipeline)
# ==========================================

def run_prediction_pipeline(fasta_path: str, output_csv: str, batch_size: int = 32, max_len: int = 1024):
    print("Ładowanie danych...")
    dataset = ProteinDataset(fasta_path, max_len=max_len)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    
    # Inicjalizacja trzech modeli - po jednym dla każdej ontologii (BP, MF, CC)
    models = {
        ontology: DummyPredictor(
            vocab_size=VOCAB_SIZE, 
            max_len=max_len, 
            num_go_terms=len(terms)
        )
        for ontology, terms in GO_VOCABULARY.items()
    }
    
    # Przełączamy modele w tryb ewaluacji (choć dla atrapy nie ma to znaczenia, to dobra praktyka)
    for model in models.values():
        model.eval()

    results = []
    
    print("Rozpoczęcie predykcji...")
    with torch.no_grad():
        for batch_ids, batch_x, batch_mask in dataloader:
            
            # batch_x ma kształt: (B, 22, 1024)
            
            for ontology, model in models.items():
                go_terms = GO_VOCABULARY[ontology]
                
                # Zwraca losowe wyniki (B, P)
                scores = model(batch_x, mask=batch_mask)
                
                # Zapisujemy każdy wynik (Białko -> GO Term -> Score)
                for i, prot_id in enumerate(batch_ids):
                    for j, term in enumerate(go_terms):
                        results.append({
                            "ProteinID": prot_id,
                            "GOterm": term,
                            "Score": round(scores[i, j].item(), 4) # Zaokrąglenie do 4 miejsc
                        })
                        
    # Eksport do pliku CSV
    print(f"Zapisywanie wyników do {output_csv}...")
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print("Gotowe!")

# ==========================================
# 5. URUCHOMIENIE
# ==========================================

if __name__ == "__main__":
    # Testowe uruchomienie rurociągu.
    # Upewnij się, że masz plik fasta o nazwie 'test_proteins.fasta' (np. ten z poprzedniego zadania)
    INPUT_FASTA = "data/processed/proteins.fasta" 
    OUTPUT_CSV = "predictions.csv"
    
    if Path(INPUT_FASTA).exists():
        run_prediction_pipeline(
            fasta_path=INPUT_FASTA,
            output_csv=OUTPUT_CSV,
            batch_size=32,
            max_len=1024
        )
    else:
        print(f"Błąd: Plik {INPUT_FASTA} nie istnieje. Podaj poprawną ścieżkę do pliku FASTA.")