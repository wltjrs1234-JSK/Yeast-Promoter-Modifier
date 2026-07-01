import sys
sys.path.append(".")
import data_loader
import analyzer
import numpy as np

def sim_index(seq, sites, max_act, ka, max_rep, kr, baseline, tata_weight):
    seq_len = len(seq)
    
    tata_sites = [s for s in sites if s["tf_id"] == "TATA"]
    if tata_sites:
        best_tata = max(tata_sites, key=lambda x: x["score"])
        tata_pos = best_tata["start"]
        tata_score = best_tata["score"]
        tata_dist_mult = analyzer.get_tata_distance_multiplier(tata_pos, seq_len)
        tata_eff = tata_score * tata_dist_mult
    else:
        tata_pos = seq_len - 120
        tata_eff = 0.15
        
    kozak_sites = [s for s in sites if s["tf_id"] == "KOZAK"]
    kozak_score = kozak_sites[0]["score"] if kozak_sites else 0.5
    
    act_sum = 0.0
    rep_sum = 0.0
    
    rap1_sites = [s for s in sites if s["tf_id"] == "RAP1"]
    gcr1_sites = [s for s in sites if s["tf_id"] == "GCR1"]
    has_synergy = any(abs(r["start"] - g["start"]) <= 80 for r in rap1_sites for g in gcr1_sites)
    synergy_multiplier = 1.5 if has_synergy else 1.0
    
    for s in sites:
        if s["tf_id"] in ["TATA", "KOZAK"]:
            continue
        acc = analyzer.get_chromatin_accessibility(s["start"], seq, seq_len)
        decay = analyzer.get_activator_decay(s["start"], tata_pos) if s["type"] == "activator" else 1.0
        eff_score = s["score"] * acc * decay
        
        if s["type"] == "activator":
            act_sum += s["weight"] * eff_score * synergy_multiplier
        elif s["type"] == "repressor":
            rep_sum += s["weight"] * eff_score
            
    # Saturation functions
    sat_act = max_act * (act_sum / (act_sum + ka)) if act_sum > 0 else 0.0
    sat_rep = max_rep * (rep_sum / (rep_sum + kr)) if rep_sum > 0 else 0.0
    
    transcription = baseline + (tata_weight * tata_eff) + sat_act - sat_rep
    transcription = max(1.0, transcription)
    
    kozak_mult = 0.5 + (0.8 * kozak_score)
    return transcription * kozak_mult

def run_simulation():
    genes = {
        "YGR192C": "TDH3 (Strong: 100)",
        "YAL003W": "TEF1 (Strong: 80)",
        "YCR012W": "PGK1 (Strong: 70)",
        "YOL086C": "ADH1 (Medium: 40)",
        "YFL039C": "ACT1 (Weak: 15)",
        "YJR048W": "CYC1 (Weak: 5)"
    }
    
    # Load all promoter data
    data_dict = {}
    for sys_id in genes:
        d = data_loader.get_promoter_data(sys_id)
        if d:
            data_dict[sys_id] = {
                "seq": d["seq"],
                "sites": analyzer.scan_promoter_motifs(d["seq"])
            }
            
    # Try different parameter sets
    # Params: max_act, ka, max_rep, kr, baseline, tata_weight
    param_sets = [
        (80.0, 300.0, 50.0, 150.0, 5.0, 15.0),
        (50.0, 400.0, 30.0, 200.0, 5.0, 30.0),
        (40.0, 500.0, 20.0, 300.0, 1.0, 50.0), # Heavy TATA dependence
        (100.0, 600.0, 40.0, 200.0, 2.0, 20.0),
    ]
    
    for idx, p in enumerate(param_sets):
        max_act, ka, max_rep, kr, baseline, tata_weight = p
        print(f"\n--- Parameter Set {idx+1}: max_act={max_act}, ka={ka}, baseline={baseline}, tata={tata_weight} ---")
        scores = []
        for sys_id, name in genes.items():
            g_data = data_dict[sys_id]
            score = sim_index(g_data["seq"], g_data["sites"], max_act, ka, max_rep, kr, baseline, tata_weight)
            scores.append((name, score))
            
        # Sort by score descending to see the predicted hierarchy
        scores.sort(key=lambda x: x[1], reverse=True)
        for name, score in scores:
            print(f"  {name}: {score:.3f}")

if __name__ == "__main__":
    run_simulation()
