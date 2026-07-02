import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
import numpy as np
import random

from deep_models import LegNetPredictor, PromoterDataset, seq_to_one_hot
from deep_generator import PromoterGeneticOptimizer, has_tata_box

def generate_synthetic_data(num_samples=100, seq_len=1000):
    """Generate mock promoter dataset representing S. cerevisiae biological constraints."""
    random.seed(42)
    np.random.seed(42)
    
    sequences = []
    expressions = []
    
    bases = ['A', 'C', 'G', 'T']
    
    for _ in range(num_samples):
        # 1. Start with a random background sequence
        seq = "".join(random.choice(bases) for _ in range(seq_len))
        
        # 2. Add structural parameters that affect expression level
        expr_val = 20.0  # weak default
        
        # Implant TATA Box at random position in the core promoter region
        tata_pos = random.randint(850, 920)
        seq_list = list(seq)
        
        # 40% chance of high expression due to active TATA Box
        if random.random() < 0.5:
            seq_list[tata_pos:tata_pos + 8] = list("TATATAAA")
            expr_val += 50.0
            
        # Implant Activator motif (GCN4: TGACTC) in upstream region
        if random.random() < 0.6:
            act_pos = random.randint(400, 700)
            seq_list[act_pos:act_pos + 6] = list("TGACTC")
            expr_val += 25.0
            
        seq = "".join(seq_list)
        
        # Normalize range to [5.0, 95.0]
        expr_val = max(5.0, min(95.0, expr_val))
        
        sequences.append(seq)
        expressions.append(expr_val)
        
    return sequences, expressions

def train_predictor(model, dataset, epochs=5, batch_size=8, device="cpu"):
    """Train the soft-classification predictor engine using KL-Divergence loss."""
    model.to(device)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    # We train soft-classification distribution using KL-Divergence loss
    criterion = nn.KLDivLoss(reduction="batchmean")
    optimizer = optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    
    print(f"[*] Training LegNetPredictor on {device} for {epochs} epochs...")
    model.train()
    for epoch in range(epochs):
        epoch_loss = 0.0
        for x, y_soft, y_scalar in dataloader:
            x, y_soft = x.to(device), y_soft.to(device)
            
            optimizer.zero_grad()
            log_probs, _ = model(x)
            
            loss = criterion(log_probs, y_soft)
            loss.backward()
            optimizer.step()
            
            epoch_loss += loss.item() * x.size(0)
            
        avg_loss = epoch_loss / len(dataset)
        print(f"  - Epoch {epoch+1:02d}/{epochs:02d} | Avg KL-Loss: {avg_loss:.4f}")
        
    print("[*] Training completed successfully.")

import sys

def main():
    print("=========================================================")
    print("    Yeast Promoter In Silico Deep Learning Toolkit       ")
    print("=========================================================\n")
    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # Check if fast diagnostic test mode is requested
    fast_test = "--fast-test" in sys.argv or "-t" in sys.argv
    
    num_samples = 4 if fast_test else 120
    epochs = 1 if fast_test else 5
    pop_size = 3 if fast_test else 30
    generations = 2 if fast_test else 30
    
    if fast_test:
        print("[!] Running in Fast Diagnostic Test Mode...")
    
    # 1. Generate synthetic dataset
    print("[+] Simulating yeast promoter training dataset...")
    seqs, exprs = generate_synthetic_data(num_samples=num_samples, seq_len=1000)
    dataset = PromoterDataset(seqs, exprs, max_len=1000, num_bins=18)
    print(f"  - Successfully created {len(dataset)} promoter-expression sample pairs.")
    
    # 2. Construct model and train
    model = LegNetPredictor(seq_len=1000, num_bins=18)
    train_predictor(model, dataset, epochs=epochs, batch_size=16 if not fast_test else 2, device=device)
    
    # 3. Test scalar prediction (inference dot product)
    print("\n[+] Testing continuous expression prediction engine...")
    test_seq = seqs[0]
    x_test = seq_to_one_hot(test_seq, max_len=1000).unsqueeze(0).to(device)
    
    model.eval()
    predicted_val = model.predict(x_test).item()
    print(f"  - True simulated expression: {exprs[0]:.2f}%")
    print(f"  - Model predicted expression (Softmax Dot Product): {predicted_val:.2f}%")
    
    # 4. Run Sequence Generator via Genetic Algorithm
    target_expression = 75.0  # Goal: design a promoter yielding exactly 75.0% strength
    print(f"\n[+] Initializing sequence generator towards target strength: {target_expression}%")
    
    # Use a random background sequence as starter template
    initial_seq = generate_synthetic_data(num_samples=1)[0][0]
    
    optimizer = PromoterGeneticOptimizer(
        predictor=model,
        target_expression=target_expression,
        base_sequence=initial_seq,
        pop_size=pop_size,
        generations=generations,
        device=device
    )
    
    best_seq, best_pred, best_mse = optimizer.run_optimization()
    
    # Print results
    print("\n=================== Design Output ===================")
    print(f"Target Expression   : {target_expression:.2f}%")
    print(f"Designed Expression : {best_pred:.2f}% (MSE: {best_mse:.4f})")
    print(f"Has TATA Box        : {'YES' if has_tata_box(best_seq) else 'NO'}")
    print(f"First 100bp of Seq  : {best_seq[:100]}...")
    print(f"Core TATA region    : ...{best_seq[800:970]}...")
    print("=====================================================")

if __name__ == "__main__":
    main()
