import requests
import re

def list_root_chrome():
    url = "http://sgd-archive.yeastgenome.org/sequence/S288C_reference/"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8",
        "Accept-Language": "ko-KR,ko;q=0.9,en-US;q=0.8,en;q=0.7",
        "Connection": "keep-alive"
    }
    res = requests.get(url, headers=headers, timeout=15)
    if res.status_code == 200:
        html = res.text
        links = re.findall(r'href="([^"]+)"', html)
        print("Links found on root:")
        for link in links:
            print(link)
    else:
        print(f"Error: {res.status_code}")
        # Print first 200 chars of response to see what's wrong
        print(res.text[:200])

if __name__ == "__main__":
    list_root_chrome()
