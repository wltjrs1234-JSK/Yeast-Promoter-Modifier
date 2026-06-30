import requests

def check_ensembl():
    systematic_id = "YGR192C"
    url = f"https://rest.ensembl.org/sequence/id/{systematic_id}?upstream=1000"
    headers = {"Content-Type": "application/json"}
    res = requests.get(url, headers=headers, timeout=10, verify=False)
    if res.status_code == 200:
        data = res.json()
        seq = data.get("seq", "").upper()
        print(f"Ensembl YGR192C upstream 1000bp (length {len(seq)}):")
        print(seq[:100] + "..." + seq[-100:])
    else:
        print(f"Error: {res.status_code}, {res.text}")

if __name__ == "__main__":
    check_ensembl()
