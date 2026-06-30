import requests

def check_ensembl_genomic():
    systematic_id = "YGR192C" # TDH3
    url = f"https://rest.ensembl.org/sequence/id/{systematic_id}?upstream=1000;type=genomic"
    headers = {"Content-Type": "application/json"}
    res = requests.get(url, headers=headers, timeout=10, verify=False)
    if res.status_code == 200:
        data = res.json()
        seq = data.get("seq", "").upper()
        print(f"Ensembl YGR192C genomic upstream 1000bp (length {len(seq)}):")
        print(seq[:100] + "..." + seq[-100:])
    else:
        print(f"Error: {res.status_code}, {res.text}")

if __name__ == "__main__":
    check_ensembl_genomic()
