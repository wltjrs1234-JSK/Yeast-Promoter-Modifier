import requests
import re

def list_root():
    url = "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/"
    res = requests.get(url, timeout=15)
    if res.status_code == 200:
        html = res.text
        links = re.findall(r'href="([^"]+)"', html)
        print("Links found on root:")
        for link in links:
            print(link)
    else:
        print(f"Error: {res.status_code}")

if __name__ == "__main__":
    list_root()
