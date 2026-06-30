import requests

def check_lookup_region():
    systematic_id = "YGR192C" # TDH3
    
    # 1. Lookup gene coordinates
    lookup_url = f"https://rest.ensembl.org/lookup/id/{systematic_id}"
    headers = {"Content-Type": "application/json"}
    
    print(f"Looking up coordinates for {systematic_id}...")
    res = requests.get(lookup_url, headers=headers, timeout=10, verify=False)
    if res.status_code != 200:
        print(f"Lookup failed: {res.status_code}, {res.text}")
        return
        
    gene_info = res.json()
    chrom = gene_info.get("seq_region_name")
    start = gene_info.get("start")
    end = gene_info.get("end")
    strand = gene_info.get("strand")
    
    print(f"Coordinates found: Chromosome {chrom}, Start {start}, End {end}, Strand {strand}")
    
    # 2. Calculate promoter range (1000bp upstream)
    # Strand 1: upstream is before start
    # Strand -1: upstream is after end
    if strand == 1:
        promoter_start = start - 1000
        promoter_end = start - 1
    else:
        promoter_start = end + 1
        promoter_end = end + 1000
        
    print(f"Calculated 1000bp upstream range: {promoter_start}..{promoter_end}")
    
    # 3. Retrieve genomic region sequence
    # Format: /sequence/region/saccharomyces_cerevisiae/Chromosome:start..end:strand
    region_url = f"https://rest.ensembl.org/sequence/region/saccharomyces_cerevisiae/{chrom}:{promoter_start}..{promoter_end}:{strand}"
    
    print(f"Fetching sequence from: {region_url}")
    res = requests.get(region_url, headers=headers, timeout=10, verify=False)
    if res.status_code != 200:
        print(f"Region query failed: {res.status_code}, {res.text}")
        return
        
    seq_data = res.json()
    seq = seq_data.get("seq", "").upper()
    print(f"Fetched sequence (length {len(seq)}):")
    print(seq[:100] + "..." + seq[-100:])
    
    # Check match with expected TDH3 promoter (TTCATCCTTT...)
    expected_start = "TTCATCCTTT"
    if seq.startswith(expected_start):
        print("SUCCESS! The sequence matches the TDH3 promoter perfectly!")
    else:
        print("FAILURE! The sequence does not match.")
        
if __name__ == "__main__":
    check_lookup_region()
