import random
import re
import torch
from deep_models import seq_to_one_hot

# Core TATA Box pattern for S. cerevisiae (typically TATAAA / TATATAAA)
TATA_PATTERN = r"TATA[AT]A[AT][AG]"

def has_tata_box(seq, start_offset=800, end_offset=970):
    """Check if the promoter sequence contains a functional TATA box in the core region.
    For S. cerevisiae, this is typically between 30bp and 120bp upstream of the TSS (positions 800 to 970 in a 1000bp sequence).
    """
    core_region = seq[start_offset:end_offset].upper()
    return bool(re.search(TATA_PATTERN, core_region))

def get_random_sequence(length=1000):
    """Generate a random DNA sequence of given length."""
    return "".join(random.choice("ACGT") for _ in range(length))

def mutate_sequence(seq, mutation_count=random.randint(1, 3)):
    """Introduce random point mutations to the DNA sequence, with a guided chance of TATA-implantation."""
    seq_list = list(seq)
    seq_len = len(seq)
    
    # Guided biological mutation (15% chance to implant TATA box in core region 850-920)
    # This helps GA find active promoters without brute-forcing a 8bp sequence
    if random.random() < 0.15 and seq_len >= 950:
        tata_pos = random.randint(850, 920)
        seq_list[tata_pos:tata_pos+8] = list("TATATAAA")
        # Reduce subsequent random mutation count
        mutation_count = max(0, mutation_count - 1)
        
    if mutation_count > 0:
        # Select unique random positions to mutate
        positions = random.sample(range(seq_len), min(mutation_count, seq_len))
        bases = ['A', 'C', 'G', 'T']
        
        for pos in positions:
            orig = seq_list[pos].upper()
            choices = [b for b in bases if b != orig]
            seq_list[pos] = random.choice(choices)
        
    return "".join(seq_list)

def crossover_sequences(parent1, parent2):
    """Perform single-point crossover between two parent sequences."""
    if len(parent1) != len(parent2):
        raise ValueError("Parents must have identical length for crossover.")
    
    seq_len = len(parent1)
    pt = random.randint(100, seq_len - 100)  # Avoid crossover at the extreme edges
    
    child1 = parent1[:pt] + parent2[pt:]
    child2 = parent2[:pt] + parent1[pt:]
    return child1, child2

def evaluate_fitness(seq, predictor, target_expr, device="cpu"):
    """Evaluate fitness of a single promoter sequence.
    Fitness is based on MSE difference from target expression, with severe TATA-loss penalty.
    """
    # 1. One-hot encode and run predictor inference
    x = seq_to_one_hot(seq, max_len=len(seq)).unsqueeze(0).to(device)
    pred_val = predictor.predict(x).item()
    
    # 2. Check biological constraints: TATA-loss penalty
    # If the sequence lacks a functional TATA box in the core promoter region,
    # it receives a severe penalty (we drop its predicted strength by 50%)
    if not has_tata_box(seq):
        pred_val *= 0.5
        
    # Fitness is the negative Mean Squared Error (closer to 0 is better)
    mse = (pred_val - target_expr) ** 2
    return pred_val, mse

class PromoterGeneticOptimizer:
    def __init__(self, predictor, target_expression, base_sequence=None, pop_size=50, generations=100, device="cpu"):
        self.predictor = predictor
        self.target_expr = target_expression
        self.base_seq = base_sequence if base_sequence else get_random_sequence(1000)
        self.seq_len = len(self.base_seq)
        self.pop_size = pop_size
        self.generations = generations
        self.device = device
        
    def run_optimization(self):
        """Execute Genetic Algorithm to evolve promoter sequence towards the target expression level."""
        print(f"[*] Starting GA optimization loop for Target Expression: {self.target_expr}%")
        
        # 1. Initialize population (clone mutated variants of the base sequence)
        population = [self.base_seq]
        for _ in range(self.pop_size - 1):
            # Mutate base sequence to create diversity
            population.append(mutate_sequence(self.base_seq, mutation_count=random.randint(5, 15)))
            
        best_overall_seq = self.base_seq
        best_overall_mse = float("inf")
        best_overall_pred = 0.0
        
        for gen in range(self.generations):
            # Evaluate fitness of all individuals
            evals = []
            for seq in population:
                pred, mse = evaluate_fitness(seq, self.predictor, self.target_expr, self.device)
                evals.append((seq, pred, mse))
                
            # Sort population by MSE (ascending - lower is better)
            evals.sort(key=lambda x: x[2])
            
            # Record best of generation
            best_seq, best_pred, best_mse = evals[0]
            if best_mse < best_overall_mse:
                best_overall_seq = best_seq
                best_overall_mse = best_mse
                best_overall_pred = best_pred
                
            if (gen + 1) % 10 == 0 or gen == 0:
                has_tata = "YES" if has_tata_box(best_seq) else "NO"
                print(f"  - Gen {gen+1:02d}/{self.generations:02d} | Best MSE: {best_mse:.2f} | Predicted: {best_pred:.2f}% | Has TATA: {has_tata}")
                
            # Elitism: carry top 20% directly to next generation
            elite_count = max(2, int(self.pop_size * 0.20))
            elites = [x[0] for x in evals[:elite_count]]
            
            next_generation = list(elites)
            
            # Fill the rest of the population
            while len(next_generation) < self.pop_size:
                # Select parents from elite pool
                p1 = random.choice(elites)
                p2 = random.choice(elites)
                
                # Crossover
                if random.random() < 0.7:
                    c1, c2 = crossover_sequences(p1, p2)
                else:
                    c1, c2 = p1, p2
                    
                # Mutation
                c1 = mutate_sequence(c1, mutation_count=random.randint(1, 3))
                c2 = mutate_sequence(c2, mutation_count=random.randint(1, 3))
                
                next_generation.append(c1)
                if len(next_generation) < self.pop_size:
                    next_generation.append(c2)
                    
            population = next_generation[:self.pop_size]
            
        print("[*] Genetic optimization finished.")
        return best_overall_seq, best_overall_pred, best_overall_mse
