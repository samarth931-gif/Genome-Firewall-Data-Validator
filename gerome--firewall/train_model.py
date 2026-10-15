# train_model.py
import os
import glob
import numpy as np
import pandas as pd
from collections import Counter
from Bio import SeqIO
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import joblib

K_SIZE = 6
TOP_KMERS_COUNT = 500

def get_kmers(sequence, k=K_SIZE):
    """Extract k-mer counts from a DNA sequence."""
    seq = str(sequence).upper()
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1) if "N" not in seq[i:i+k]]
    return Counter(kmers)

def build_dataset(fasta_dir, labels_csv):
    """
    fasta_dir: folder containing .fasta files
    labels_csv: CSV mapping filename to drug status (1 = Resistant, 0 = Susceptible)
    """
    labels_df = pd.read_csv(labels_csv) # Columns: ["filename", "is_resistant"]
    
    dataset = []
    for _, row in labels_df.iterrows():
        filepath = os.path.join(fasta_dir, row["filename"])
        if not os.path.exists(filepath):
            continue
            
        full_seq = ""
        for record in SeqIO.parse(filepath, "fasta"):
            full_seq += str(record.seq)
            
        counts = get_kmers(full_seq, k=K_SIZE)
        counts["_target_"] = row["is_resistant"]
        dataset.append(counts)
        
    df = pd.DataFrame(dataset).fillna(0)
    return df

# --- Dummy/Synthetic Training Workflow for Fast Iteration ---
if __name__ == "__main__":
    print("Initializing trainer...")
    
    # Install dependencies: pip install scikit-learn joblib numpy
    # Build synthetic k-mer feature matrix (100 sample genomes, 500 top k-mer features)
    X_synthetic = np.random.randint(0, 50, size=(100, TOP_KMERS_COUNT))
    y_synthetic = np.random.choice([0, 1], size=(100,))
    
    # Generate canonical feature column names
    feature_names = [f"kmer_{i}" for i in range(TOP_KMERS_COUNT)]
    
    X_train, X_test, y_train, y_test = train_test_split(X_synthetic, y_synthetic, test_size=0.2, random_state=42)
    
    clf = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=42)
    clf.fit(X_train, y_train)
    
    preds = clf.predict(X_test)
    print("Model Validation Report:\n", classification_report(y_test, preds))
    
    # Save model and feature schema
    model_payload = {
        "model": clf,
        "feature_names": feature_names,
        "k_size": K_SIZE
    }
    joblib.dump(model_payload, "amr_model.pkl")
    print("Model saved to amr_model.pkl")