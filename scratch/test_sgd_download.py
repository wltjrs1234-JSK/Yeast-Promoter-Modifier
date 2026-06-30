import os
import json
import gzip
import requests

CACHE_DIR = "data_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

def test_sgd():
    url = "https://downloads.yeastgenome.org/sequence/S288C_reference/orf_dna/orf_genomic_1000_upstream.fasta.gz"
    print("Downloading SGD Fasta...")
    res = requests.get(url, timeout=30, verify=False)
    if res.status_code == 200:
        print("Success! Parsing...")
        decompressed = gzip.decompress(res.content).decode("utf-8")
        lines = decompressed.splitlines()
        print(f"Total lines: {len(lines)}")
        
        sgd_promoters = {}
        current_id = None
        current_seq = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id and current_seq:
                    sgd_promoters[current_id] = "".join(current_seq).upper()
                parts = line[1:].split()
                if parts:
                    current_id = parts[0].strip().upper()
                current_seq = []
            else:
                current_seq.append(line)
        if current_id and current_seq:
            sgd_promoters[current_id] = "".join(current_seq).upper()
            
        print(f"Total parsed: {len(sgd_promoters)}")
        # Check YGR192C (TDH3)
        if "YGR192C" in sgd_promoters:
            print("YGR192C found!")
            print(sgd_promoters["YGR192C"][:100] + "...")
        else:
            print("YGR192C NOT found!")
    else:
        print(f"Failed to download: {res.status_code}")

if __name__ == "__main__":
    test_sgd()
