import skbio
import torch
import torch.nn as nn
import pandas as pd
from torch.utils.data import Dataset, DataLoader
from pathlib import Path


torch.manual_seed(42)

#CONFIG

# 20 standardowych aminokwasów + X (nieznany) + PAD (padding)
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY") + ["X", "PAD"]
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
VOCAB_SIZE = len(AMINO_ACIDS)  # N = 22
PAD_IDX = AA_TO_IDX["PAD"]

# Sample, hardocoded GO terms for each ontology (for demonstration purposes). SHOULD BE CHANGED LATER
GO_VOCABULARY = {
    "BP": ["GO:0008150", "GO:0009987", "GO:0008152"],  # Biological Process
    "MF": ["GO:0003674", "GO:0005488"],                # Molecular Function
    "CC": ["GO:0005575", "GO:0005623", "GO:0005886", "GO:0043226"]  # Cellular Component
}

#LOAD DATASET

class ProteinDataset(Dataset):
    def __init__(self, fasta_file: str, max_len: int = 1024):
        """
        Reads a FASTA file, encodes sequences, and prepares masks.
        """
        self.max_len = max_len
        self.proteins = []

        for seq in skbio.io.read(fasta_file, format='fasta'):
            self.proteins.append({
                'id': seq.metadata['id'],
                'sequence': str(seq)
            })

    def __len__(self):
        return len(self.proteins)

    def __getitem__(self, idx):
        prot = self.proteins[idx]
        seq_str = prot['sequence'].upper()  #make sure all amino acids are uppercase
        
        # Truncating too long sequences
        seq_str = seq_str[:self.max_len]

        #Initializing tensor for sequence of shape (N, D)
        #Where N is VOCAB_SIZE, and D is max_len
        seq_tensor = torch.zeros((VOCAB_SIZE, self.max_len))
        mask = torch.zeros(self.max_len, dtype=torch.bool)

        # One-Hot Encoding and setting mask to True for valid amino acids
        for i, aa in enumerate(seq_str):
            aa_idx = AA_TO_IDX.get(aa, AA_TO_IDX["X"]) # X if amino acid is unknown
            seq_tensor[aa_idx, i] = 1.0
            mask[i] = True 
            
        # Padding for shorter sequences 
        for i in range(len(seq_str), self.max_len):
            seq_tensor[PAD_IDX, i] = 1.0
            mask[i] = False #False means that this position is padding and should be ignored in the model
            
        return prot['id'], seq_tensor, mask


# DUMMY MODEL 

class DummyPredictor(nn.Module):
    def __init__(self, vocab_size: int, max_len: int, num_go_terms: int):
        super().__init__()
        self.num_go_terms = num_go_terms
        
        # Empty layers - dummy model.
        # We save parameters to ensure the interface matches the target model.

        self.vocab_size = vocab_size
        self.max_len = max_len

    def forward(self, x, mask=None):
        """
        
        Expected input: 
        x shape: (B, N, D) -> Batch Size, Vocab Size, Max Length
        mask shape: (B, D)
        
        As per the instruction, the interface accepts a mask, but the dummy model
        simply returns random results in shape (B, P).
        """
        batch_size = x.size(0)
        
        # Generating random predictions in the range 0.0 - 1.0
        # Output shape: (B, P), where P is the number of GO terms for the given ontology
        random_scores = torch.rand((batch_size, self.num_go_terms), device=x.device)
        
        return random_scores


# END TO END PIPELINE


def run_prediction_pipeline(fasta_path: str, output_csv: str, batch_size: int = 32, max_len: int = 1024):
    print("Loading data...")
    dataset = ProteinDataset(fasta_path, max_len=max_len)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    
    # Initialization of three models - one for each ontology (BP, MF, CC)
    models = {
        ontology: DummyPredictor(
            vocab_size=VOCAB_SIZE, 
            max_len=max_len, 
            num_go_terms=len(terms)
        )
        for ontology, terms in GO_VOCABULARY.items()
    }
    
    # We set the models to evaluation mode (though for the dummy model it doesn't matter, it's good practice)
    for model in models.values():
        model.eval()

    results = []
    
    print("Prediction starts...")
    with torch.no_grad():
        for batch_ids, batch_x, batch_mask in dataloader:
            
            # batch_x has shape: (B, 22, 1024)
            
            for ontology, model in models.items():
                go_terms = GO_VOCABULARY[ontology]
                
                # Returns random scores (B,P) where B is batch size and P is number of GO terms for the ontology
                scores = model(batch_x, mask=batch_mask)
                
                # We save each result (Protein -> GO Term -> Score)
                for i, prot_id in enumerate(batch_ids):
                    for j, term in enumerate(go_terms):
                        results.append({
                            "ProteinID": prot_id,
                            "GOterm": term,
                            "Score": round(scores[i, j].item(), 4) # Rounding to 4 decimal places for better readability
                        })
                        
    # Export to CSV file
    print(f"Writing results to {output_csv}...")
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print("Results written to CSV")

# RUN PIPELINE

if __name__ == "__main__":
    # Test run of the pipeline

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
        print(f"Error: File {INPUT_FASTA} does not exist. Please provide a valid path to the FASTA file.")