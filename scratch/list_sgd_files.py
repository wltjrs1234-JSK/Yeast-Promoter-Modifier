import requests
import re

def list_files():
    url = "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/"
    res = requests.get(url, timeout=15)
    if res.status_code == 200:
        html = res.text
        # Find all href links
        links = re.findall(r'href="([^"]+)"', html)
        print("Files in http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_dna/:")
        for link in links:
            if not link.startswith("..") and not link.startswith("?"):
                print(link)
    else:
        print(f"Error: {res.status_code}")

if __name__ == "__main__":
    list_files()
