import sys
sys.path.append(".")
import data_loader
import analyzer

def debug_gene(sys_id):
    data = data_loader.get_promoter_data(sys_id)
    if not data:
        print(f"Failed to load {sys_id}")
        return
    seq = data["seq"]
    sites = analyzer.scan_promoter_motifs(seq)
    print(f"\n=== DEBUG {data['symbol']} ({sys_id}) ===")
    print(f"Total sites scanned: {len(sites)}")
    
    tata_sites = [s for s in sites if s["tf_id"] == "TATA"]
    kozak_sites = [s for s in sites if s["tf_id"] == "KOZAK"]
    activators = [s for s in sites if s["type"] == "activator" and s["tf_id"] not in ["TATA", "KOZAK"]]
    repressors = [s for s in sites if s["type"] == "repressor"]
    
    print(f"TATA Sites: {len(tata_sites)}")
    for t in tata_sites:
        print(f"  TATA at {t['start']}..{t['end']}, score={t['score']:.2f}")
        
    print(f"Kozak Sites: {len(kozak_sites)}")
    for k in kozak_sites:
        print(f"  Kozak score={k['score']:.2f}")
        
    print(f"Activators: {len(activators)}")
    # Print top 5 activators by contribution (weight * score)
    activators.sort(key=lambda x: x["weight"] * x["score"], reverse=True)
    for a in activators[:8]:
        print(f"  {a['tf_name']} at {a['start']}, weight={a['weight']}, score={a['score']:.2f}")
        
    print(f"Repressors: {len(repressors)}")
    for r in repressors:
        print(f"  {r['tf_name']} at {r['start']}, weight={r['weight']}, score={r['score']:.2f}")
        
    raw_score = analyzer.calculate_absolute_expression_index(seq, sites)
    print(f"=> Calculated Raw Score: {raw_score:.3f}")

if __name__ == "__main__":
    debug_gene("YFL039C") # ACT1
    debug_gene("YOL086C") # ADH1
