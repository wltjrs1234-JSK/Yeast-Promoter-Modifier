import requests
import json
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def test():
    url = "https://rest.ensembl.org/sequence/id/YGR192C?upstream=1000"
    headers = {"Content-Type": "application/json"}
    try:
        res = requests.get(url, headers=headers, timeout=10, verify=False)
        print("Status Code:", res.status_code)
        if res.status_code == 200:
            data = res.json()
            print("Sequence Length:", len(data.get("seq", "")))
            print("Sequence Preview:", data.get("seq", "")[:50] + "...")
            print("ID:", data.get("id"))
        else:
            print("Error response:", res.text)
    except Exception as e:
        print("Exception:", e)

if __name__ == "__main__":
    test()
