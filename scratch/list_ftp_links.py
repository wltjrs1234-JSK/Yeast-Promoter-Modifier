import requests
import re

def list_ftp():
    url = "https://ftp.yeastgenome.org/pub/yeast/sequence/S288C_reference/"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
    }
    res = requests.get(url, headers=headers, timeout=15, verify=False)
    if res.status_code == 200:
        html = res.text
        links = re.findall(r'href="([^"]+)"', html)
        print("Links found on FTP:")
        for link in links:
            print(link)
    else:
        print(f"Error: {res.status_code}")

if __name__ == "__main__":
    list_ftp()
