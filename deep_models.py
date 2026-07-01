import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset
import numpy as np

# 1. Data Preprocessing & Encoding Module
DNA_MAP = {
    'A': [1.0, 0.0, 0.0, 0.0],
    'C': [0.0, 1.0, 0.0, 0.0],
    'G': [0.0, 0.0, 1.0, 0.0],
    'T': [0.0, 0.0, 0.0, 1.0],
    'N': [0.0, 0.0, 0.0, 0.0]
}

def seq_to_one_hot(seq, max_len=1000):
    """Convert DNA sequence to a 4 x L One-hot PyTorch Tensor with Zero-padding.
    A=[1,0,0,0], C=[0,1,0,0], G=[0,0,1,0], T=[0,0,0,1], others/padding=[0,0,0,0].
    """
    seq = seq.upper()
    seq_len = len(seq)
    
    # Initialize zero-matrix
    encoding = np.zeros((4, max_len), dtype=np.float32)
    
    # Fill in base values up to max_len
    for i in range(min(seq_len, max_len)):
        base = seq[i]
        vector = DNA_MAP.get(base, DNA_MAP['N'])
        for channel in range(4):
            encoding[channel, i] = vector[channel]
            
    return torch.tensor(encoding, dtype=torch.float32)

class PromoterDataset(Dataset):
    """Custom PyTorch Dataset for yeast promoter sequences and expression strengths."""
    def __init__(self, sequences, expressions, max_len=1000, num_bins=18):
        self.sequences = sequences
        self.expressions = expressions  # Continuous expressions (e.g. 0.0 to 100.0)
        self.max_len = max_len
        self.num_bins = num_bins
        
        # Pre-calculate bin edges (for 0 to 100 range)
        self.bin_edges = np.linspace(0.0, 100.0, num_bins + 1)
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2.0
        
    def __len__(self):
        return len(self.sequences)
        
    def _expression_to_soft_target(self, val):
        """Convert continuous value to soft-target distribution over bins (for soft classification)."""
        # Simple Gaussian smoothing around the target value to create soft targets
        sigma = 5.0
        diffs = self.bin_centers - val
        probs = np.exp(-(diffs ** 2) / (2 * sigma ** 2))
        probs = probs / np.sum(probs)
        return torch.tensor(probs, dtype=torch.float32)
        
    def __getitem__(self, idx):
        seq = self.sequences[idx]
        expr = self.expressions[idx]
        
        x = seq_to_one_hot(seq, max_len=self.max_len)
        y_soft = self._expression_to_soft_target(expr)
        y_scalar = torch.tensor(expr, dtype=torch.float32)
        
        return x, y_soft, y_scalar

# 2. Squeeze-and-Excitation (SE) Block for channel recalibration
class SqueezeExcitation1d(nn.Module):
    def __init__(self, channels, reduction=4):
        super(SqueezeExcitation1d, self).__init__()
        self.fc1 = nn.Conv1d(channels, channels // reduction, kernel_size=1)
        self.fc2 = nn.Conv1d(channels // reduction, channels, kernel_size=1)
        
    def forward(self, x):
        # Squeeze: Global Average Pooling
        w = F.adaptive_avg_pool1d(x, 1)
        # Excitation: MLP with ReLU + Sigmoid
        w = F.relu(self.fc1(w))
        w = torch.sigmoid(self.fc2(w))
        # Recalibrate original input
        return x * w

# 3. LegNet-inspired 1D-CNN Predictor Engine
class LegNetPredictor(nn.Module):
    def __init__(self, seq_len=1000, num_bins=18):
        super(LegNetPredictor, self).__init__()
        self.seq_len = seq_len
        self.num_bins = num_bins
        
        # Initial stem block
        self.stem = nn.Sequential(
            nn.Conv1d(in_channels=4, out_channels=64, kernel_size=7, padding=3),
            nn.BatchNorm1d(64),
            nn.SiLU()
        )
        
        # Dilated residual blocks with Squeeze-and-Excitation
        self.block1 = nn.Sequential(
            nn.Conv1d(64, 64, kernel_size=5, padding=4, dilation=2),
            nn.BatchNorm1d(64),
            nn.SiLU(),
            SqueezeExcitation1d(64)
        )
        
        self.block2 = nn.Sequential(
            nn.Conv1d(64, 128, kernel_size=5, padding=8, dilation=4),
            nn.BatchNorm1d(128),
            nn.SiLU(),
            SqueezeExcitation1d(128)
        )
        
        self.block3 = nn.Sequential(
            nn.Conv1d(128, 256, kernel_size=5, padding=16, dilation=8),
            nn.BatchNorm1d(256),
            nn.SiLU(),
            SqueezeExcitation1d(256)
        )
        
        # Global pooling and classification head
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.fc = nn.Sequential(
            nn.Linear(256, 128),
            nn.SiLU(),
            nn.Dropout(0.3),
            nn.Linear(128, num_bins)
        )
        
        # Reference bin weights (center values for 0.0 to 100.0)
        bin_edges = np.linspace(0.0, 100.0, num_bins + 1)
        self.bin_centers = torch.tensor((bin_edges[:-1] + bin_edges[1:]) / 2.0, dtype=torch.float32)
        
    def forward(self, x):
        """Returns log-probabilities (for KL-div loss) and raw softmax probabilities."""
        x = self.stem(x)
        x = x + self.block1(x)  # Residual connection
        x = self.block2(x)
        x = self.block3(x)
        
        x = self.pool(x).squeeze(-1)
        logits = self.fc(x)
        
        # Return log softmax for stable KL-Divergence loss training
        log_probs = F.log_softmax(logits, dim=-1)
        # Return standard softmax for inference
        probs = F.softmax(logits, dim=-1)
        
        return log_probs, probs
        
    def predict(self, x):
        """Predict continuous expression level via dot product of softmax probabilities and bin weights."""
        self.eval()
        with torch.no_grad():
            _, probs = self.forward(x)
            # Dot product of probabilities and bin centers
            # Ensure bin centers tensor is on the same device as input
            bin_centers = self.bin_centers.to(x.device)
            predicted_vals = torch.sum(probs * bin_centers, dim=-1)
            return predicted_vals
