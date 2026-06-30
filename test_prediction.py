import analyzer
import random

def run_tests():
    print("--- Running Promoter Analyzer Tests ---")
    
    # 1. Generate random background sequence (A/C/G/T mixed) using a fixed seed for repeatability
    random.seed(42)
    bg = "".join(random.choice("ACGT") for _ in range(200))
    
    # Implant TATA box and GCN4 site
    bg_list = list(bg)
    bg_list[50:58] = list("TATATAAA")
    bg_list[120:126] = list("TGACTC")
    seq = "".join(bg_list)
    
    print("Base Sequence Length:", len(seq))
    
    # Scan motifs
    sites = analyzer.scan_promoter_motifs(seq)
    print(f"Scanned {len(sites)} sites:")
    for s in sites:
        print(f" - [{s['tf_id']}] Start: {s['start']}, End: {s['end']}, Strand: {s['strand']}, Seq: {s['sequence']}, Score: {s['score']}")
        
    assert any(s["tf_id"] == "TATA" for s in sites), "TATA box not detected!"
    assert any(s["tf_id"] == "GCN4" for s in sites), "GCN4 site not detected!"
    
    # 2. Test mutation prediction
    # Mutate TATA box (TATATAAA -> TATAGAAA) at pos 54 (third A to G)
    mutations = [{"pos": 54, "to": "G"}]
    mutant_seq = analyzer.apply_mutations(seq, mutations)
    
    print("\nPredicting mutation effect (TATA box mutation):")
    result = analyzer.predict_expression_change(seq, mutant_seq, wild_sites=sites)
    print("Predicted Expression:", result["predicted_value"], "%")
    print("Change:", result["change_percentage"], "%")
    print("TATA Destroyed:", result["tata_destroyed"])
    
    # Since TATA is destroyed, expected expression should drop heavily
    assert result["predicted_value"] < 30.0, "TATA destruction should cause expression drop!"
    
    # 3. Test recommendation engine
    print("\nGenerating Recommendations:")
    recs = analyzer.get_mutation_recommendations(seq)
    print("UP recommendations:")
    for r in recs["up"]:
        print(f" - Effect: {r['effect']}, Pos: {r['pos']}, From: {r['from']} -> To: {r['to']}")
        
    print("DOWN recommendations:")
    for r in recs["down"]:
        print(f" - Effect: {r['effect']}, Pos: {r['pos']}, From: {r['from']} -> To: {r['to']}")
        
    print("\nAll tests passed successfully!")

if __name__ == "__main__":
    run_tests()
