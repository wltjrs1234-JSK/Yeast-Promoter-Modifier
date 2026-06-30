import os
import json
import requests
import urllib3

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

CACHE_DIR = "data_cache"
os.makedirs(CACHE_DIR, exist_ok=True)

PROMOTER_CACHE_FILE = os.path.join(CACHE_DIR, "promoter_cache.json")
GENE_MAP_FILE = os.path.join(CACHE_DIR, "gene_mapping.json")

# Default pre-defined strong promoters for quick offline usage
DEFAULT_PROMOTERS = {
    "YGR192C": {
        "symbol": "TDH3",
        "description": "Glyceraldehyde-3-phosphate dehydrogenase, subunit 3; strong constitutive promoter",
        "seq": "TTCATCCTTTTTTTCTTCCTTACTAACCTACCGTGATATGAGCGCATCGCGCTCACATCGAAAGAGACGGCTTTTCCAACAATTTTGTCAAGCGACGGCCATTGAGATTGAATGGAAGAAAAGTGGGTGTGTGGAGTGGGGGTTTGGGTGACGAAAGTGAAAGCGATTGGGGTTGGTGATGAGCGTGACCAGCAGCAGCAGGACAGGTCATGGCCGGGACGACAGCAACACACGCTCGGGCACGGTTGTTCTGGGGGCGGGGGCGGGGGTGATGGCGTTGGGGGTCGTGGTGCGGCGGCTACTTTCGTCCTCAATCGGTATTTAGAAATCCCGCAATTCGTGTTTTCCTCCTTTTCTTGTTCTTTTCTTTTGGTTTTTATTATGACGTTACAAACAACAAGAACAACAACGACAACAACAACAACAAACAACAACGACAACAACAACAACAAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAACAAACAACAACAACAAC"
    }
}

def load_promoter_cache():
    if os.path.exists(PROMOTER_CACHE_FILE):
        try:
            with open(PROMOTER_CACHE_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return DEFAULT_PROMOTERS.copy()
    else:
        # Save default to create the file
        save_promoter_cache(DEFAULT_PROMOTERS)
        return DEFAULT_PROMOTERS.copy()

def save_promoter_cache(cache):
    with open(PROMOTER_CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache, f, indent=2, ensure_ascii=False)

def download_gene_mapping():
    """Download SGD_features.tab and construct standard gene symbol to systematic name map."""
    print("Downloading SGD_features.tab for gene name mapping...")
    url = "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
    backup_url = "https://ftp.yeastgenome.org/pub/yeast/curation/chromosomal_feature/SGD_features.tab"
    
    features_path = os.path.join(CACHE_DIR, "SGD_features.tab")
    
    # Try downloading
    success = False
    for download_url in [url, backup_url]:
        try:
            res = requests.get(download_url, timeout=30, verify=False)
            res.raise_for_status()
            with open(features_path, "wb") as f:
                f.write(res.content)
            success = True
            break
        except Exception as e:
            print(f"Failed to download from {download_url}: {e}")
            
    if not success:
        print("Could not download SGD_features.tab. Standard mappings might be limited.")
        return {}
        
    gene_map = {}
    try:
        with open(features_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 16:
                    feature_type = parts[1]
                    # We only care about ORFs (genes)
                    if feature_type == "ORF":
                        systematic_name = parts[3].strip()
                        standard_name = parts[4].strip()
                        aliases = [a.strip() for a in parts[5].split(",") if a.strip()]
                        desc = parts[15].strip()
                        
                        gene_data = {
                            "systematic_name": systematic_name,
                            "symbol": standard_name if standard_name else systematic_name,
                            "description": desc,
                            "aliases": aliases
                        }
                        
                        # Map by systematic name
                        gene_map[systematic_name.upper()] = gene_data
                        # Map by standard name
                        if standard_name:
                            gene_map[standard_name.upper()] = gene_data
                        # Map by aliases
                        for alias in aliases:
                            gene_map[alias.upper()] = gene_data
                            
        with open(GENE_MAP_FILE, "w", encoding="utf-8") as f:
            json.dump(gene_map, f, indent=2, ensure_ascii=False)
            
        print(f"Constructed gene mapping with {len(gene_map)} keys.")
        return gene_map
    except Exception as e:
        print(f"Error parsing SGD_features.tab: {e}")
        return {}

def get_gene_mapping():
    if os.path.exists(GENE_MAP_FILE):
        try:
            with open(GENE_MAP_FILE, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return download_gene_mapping()

def fetch_promoter_sequence(systematic_id):
    """Fetch 1000bp upstream sequence for a given yeast gene from Ensembl REST API."""
    url = f"https://rest.ensembl.org/sequence/id/{systematic_id}?upstream=1000"
    headers = {"Content-Type": "application/json"}
    try:
        res = requests.get(url, headers=headers, timeout=10, verify=False)
        if res.status_code == 200:
            data = res.json()
            return data.get("seq", "").upper()
        else:
            print(f"Ensembl API error {res.status_code}: {res.text}")
    except Exception as e:
        print(f"Exception fetching sequence: {e}")
    return None

def get_promoter_data(gene_query):
    """Retrieve promoter data (seq, symbol, description) for a gene query (Symbol or Systematic ID)."""
    gene_query = gene_query.strip().upper()
    gene_map = get_gene_mapping()
    
    if gene_query not in gene_map:
        return None
        
    gene_info = gene_map[gene_query]
    systematic_name = gene_info["systematic_name"]
    symbol = gene_info["symbol"]
    description = gene_info["description"]
    
    # Check cache first
    cache = load_promoter_cache()
    if systematic_name in cache:
        cached_data = cache[systematic_name]
        return {
            "systematic_name": systematic_name,
            "symbol": symbol,
            "description": description,
            "seq": cached_data["seq"]
        }
        
    # Fetch from API
    seq = fetch_promoter_sequence(systematic_name)
    if seq:
        # Cache the result
        cache[systematic_name] = {
            "symbol": symbol,
            "description": description,
            "seq": seq
        }
        save_promoter_cache(cache)
        
        return {
            "systematic_name": systematic_name,
            "symbol": symbol,
            "description": description,
            "seq": seq
        }
        
    return None
