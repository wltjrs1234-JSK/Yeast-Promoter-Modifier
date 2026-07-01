import sys
sys.path.append(".")
import data_loader
import analyzer
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def main():
    standard_promoters = {
        "YGR192C": {"symbol": "TDH3", "ref": 100.0},
        "YAL003W": {"symbol": "TEF1 (EFB1)", "ref": 80.0},
        "YCR012W": {"symbol": "PGK1", "ref": 70.0},
        "YOL086C": {"symbol": "ADH1", "ref": 40.0},
        "YFL039C": {"symbol": "ACT1", "ref": 15.0},
        "YJR048W": {"symbol": "CYC1", "ref": 5.0}
    }
    
    print("Calculating raw scores for standard promoters...")
    for sys_id, info in standard_promoters.items():
        data = data_loader.get_promoter_data(sys_id)
        if data:
            seq = data["seq"]
            sites = analyzer.scan_promoter_motifs(seq)
            raw_score = analyzer.calculate_absolute_expression_index(seq, sites)
            print(f"{info['symbol']} ({sys_id}): Raw Score = {raw_score:.3f}, Expected Ref = {info['ref']}")
        else:
            print(f"Failed to load promoter data for {sys_id}")

if __name__ == "__main__":
    main()
