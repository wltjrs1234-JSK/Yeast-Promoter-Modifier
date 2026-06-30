import requests

def check_ensembl_params():
    systematic_id = "YGR192C" # TDH3
    url = f"https://rest.ensembl.org/sequence/id/{systematic_id}"
    params = {
        "upstream": 1000,
        "type": "genomic"
    }
    headers = {"Content-Type": "application/json"}
    res = requests.get(url, headers=headers, params=params, timeout=10, verify=False)
    if res.status_code == 200:
        data = res.json()
        seq = data.get("seq", "").upper()
        print(f"Ensembl YGR192C genomic upstream 1000bp (length {len(seq)}):")
        print(seq[:100] + "..." + seq[-100:])
        
        # Check if it starts with the expected promoter sequence (TTCATCCTTT...)
        expected_start = "TTCATCCTTT"
        if seq.startswith(expected_start):
            print("MATCH! The sequence correctly starts with the TDH3 promoter.")
        else:
            print("MISMATCH! The sequence does not start with the TDH3 promoter.")
    else:
        print(f"Error: {res.status_code}, {res.text}")

if __name__ == "__main__":
    check_ensembl_params()
