import requests

def print_tdh3():
    # YGR192C (TDH3) promoter range: VII:883811..884810:-1
    headers = {"Content-Type": "application/json"}
    region_url = "https://rest.ensembl.org/sequence/region/saccharomyces_cerevisiae/VII:883811..884810:-1"
    res = requests.get(region_url, headers=headers, timeout=10, verify=False)
    if res.status_code == 200:
        seq = res.json().get("seq", "").upper()
        print("TDH3_SEQ_START")
        print(seq)
        print("TDH3_SEQ_END")

if __name__ == "__main__":
    print_tdh3()
