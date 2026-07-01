import re
import numpy as np

# Define transcription factor binding site details: IUPAC consensus and weights
# IUPAC Code mappings: 
# R: A/G, Y: C/T, S: G/C, W: A/T, K: G/T, M: A/C, B: C/G/T, D: A/G/T, H: A/C/T, V: A/C/G, N: A/C/G/T
IUPAC_MAP = {
    'A': '[A]', 'C': '[C]', 'G': '[G]', 'T': '[T]',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

TF_MOTIFS = {
    "TATA": {
        "name": "TATA Box",
        "consensus": "TATAWAWR",  # e.g., TATA[AT]A[AT][AG]
        "weight": 100,
        "type": "activator",
        "threshold": 0.99,  # Requires perfect match (8/8)
        "desc": "Essential core promoter element for transcription initiation."
    },
    "GAL4": {
        "name": "GAL4",
        "consensus": "CGGNNNNNNNNNNNCCG",  # CGG[ACGT]{11}CCG
        "weight": 80,
        "type": "activator",
        "threshold": 0.99,  # Requires perfect match (17/17)
        "desc": "Galactose-induction master activator. Binds galactose-responsive promoters."
    },
    "GCN4": {
        "name": "GCN4",
        "consensus": "TGACTC",  # Also TGANTCA
        "weight": 25,
        "type": "activator",
        "threshold": 0.99,  # Requires perfect match (6/6)
        "desc": "General amino acid control activator, induced under nutrient starvation."
    },
    "MSN2_4": {
        "name": "MSN2/4 (STRE)",
        "consensus": "AGGGG",  # STRE element (AGGGG / CCCCT)
        "weight": 20,
        "type": "activator",
        "threshold": 0.99,  # Requires perfect match (5/5)
        "desc": "Stress Response Element (STRE) activator, induced by environmental stresses."
    },
    "MIG1": {
        "name": "MIG1",
        "consensus": "SYGGRG",  # [GC][CT]GG[AG]G
        "weight": 60,
        "type": "repressor",
        "threshold": 0.99,  # Requires perfect match (6/6)
        "desc": "Glucose-induced repressor. Shuts down alternative carbon source utilization in glucose."
    },
    "RAP1": {
        "name": "RAP1",
        "consensus": "RMACCCANCATTG",
        "weight": 50,
        "type": "activator",
        "threshold": 0.95,  # Requires near-perfect match (12/13)
        "desc": "Constitutive strong transcriptional activator, drives ribosomal and glycolytic genes."
    },
    "GCR1": {
        "name": "GCR1",
        "consensus": "CTTCC",
        "weight": 30,
        "type": "activator",
        "threshold": 0.99,  # Requires perfect match (5/5)
        "desc": "Glycolytic genes regulator, acts cooperatively with Rap1p."
    },
    "HSF1": {
        "name": "HSF1",
        "consensus": "NGAAN",
        "weight": 25,
        "type": "activator",
        "threshold": 0.99,  # Requires perfect match (5/5)
        "desc": "Heat shock transcription factor, binds Heat Shock Elements (HSE)."
    }
}

# TF PWMs (Log-odds or relative affinity weights per position)
# Mapping: position -> nucleotide -> score
TF_PWMS = {
    "TATA": [  # consensus: TATAWAWR
        {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}, # T
        {"A": 0.8, "C": 0.1, "G": 0.0, "T": 0.1}, # A
        {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}, # T
        {"A": 0.8, "C": 0.1, "G": 0.0, "T": 0.1}, # A
        {"A": 0.5, "C": 0.0, "G": 0.0, "T": 0.5}, # W (A/T)
        {"A": 0.8, "C": 0.1, "G": 0.0, "T": 0.1}, # A
        {"A": 0.5, "C": 0.0, "G": 0.0, "T": 0.5}, # W (A/T)
        {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0}, # R (A/G)
    ],
    "GAL4": [  # consensus: CGGNNNNNNNNNNNCCG (17bp)
        {"A": 0.0, "C": 0.1, "G": 0.8, "T": 0.1},
        {"A": 0.0, "C": 0.8, "G": 0.1, "T": 0.1},
        {"A": 0.0, "C": 0.8, "G": 0.1, "T": 0.1},
        *([{"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}] * 11),
        {"A": 0.1, "C": 0.1, "G": 0.8, "T": 0.0},
        {"A": 0.1, "C": 0.1, "G": 0.8, "T": 0.0},
        {"A": 0.0, "C": 0.8, "G": 0.1, "T": 0.1},
    ],
    "GCN4": [  # consensus: TGACTC
        {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}, # T
        {"A": 0.1, "C": 0.0, "G": 0.8, "T": 0.1}, # G
        {"A": 0.8, "C": 0.1, "G": 0.1, "T": 0.0}, # A
        {"A": 0.1, "C": 0.8, "G": 0.1, "T": 0.0}, # C
        {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}, # T
        {"A": 0.1, "C": 0.8, "G": 0.1, "T": 0.0}, # C
    ],
    "MSN2_4": [ # consensus: AGGGG
        {"A": 0.8, "C": 0.0, "G": 0.1, "T": 0.1},
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1},
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1},
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1},
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1},
    ],
    "MIG1": [  # consensus: SYGGRG
        {"A": 0.0, "C": 0.5, "G": 0.5, "T": 0.0}, # S (G/C)
        {"A": 0.0, "C": 0.5, "G": 0.0, "T": 0.5}, # Y (C/T)
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1}, # G
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1}, # G
        {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0}, # R (A/G)
        {"A": 0.0, "C": 0.0, "G": 0.9, "T": 0.1}, # G
    ],
    "RAP1": [  # consensus: RMACCCANCATTG
        {"A": 0.5, "C": 0.0, "G": 0.5, "T": 0.0}, # R (A/G)
        {"A": 0.5, "C": 0.5, "G": 0.0, "T": 0.0}, # M (A/C)
        {"A": 0.8, "C": 0.1, "G": 0.1, "T": 0.0}, # A
        {"A": 0.0, "C": 0.9, "G": 0.0, "T": 0.1}, # C
        {"A": 0.0, "C": 0.9, "G": 0.0, "T": 0.1}, # C
        {"A": 0.0, "C": 0.9, "G": 0.0, "T": 0.1}, # C
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, # N
        {"A": 0.1, "C": 0.8, "G": 0.1, "T": 0.0}, # C
        {"A": 0.8, "C": 0.1, "G": 0.1, "T": 0.0}, # A
        {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}, # T
        {"A": 0.1, "C": 0.0, "G": 0.1, "T": 0.8}, # T
        {"A": 0.1, "C": 0.0, "G": 0.8, "T": 0.1}, # G
        {"A": 0.1, "C": 0.8, "G": 0.1, "T": 0.0}, # C
    ],
    "GCR1": [  # consensus: CTTCC
        {"A": 0.0, "C": 0.9, "G": 0.1, "T": 0.0}, # C
        {"A": 0.0, "C": 0.0, "G": 0.1, "T": 0.9}, # T
        {"A": 0.0, "C": 0.0, "G": 0.1, "T": 0.9}, # T
        {"A": 0.0, "C": 0.9, "G": 0.1, "T": 0.0}, # C
        {"A": 0.0, "C": 0.9, "G": 0.1, "T": 0.0}, # C
    ],
    "HSF1": [  # consensus: NGAAN
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, # N
        {"A": 0.1, "C": 0.0, "G": 0.8, "T": 0.1}, # G
        {"A": 0.8, "C": 0.1, "G": 0.1, "T": 0.0}, # A
        {"A": 0.8, "C": 0.1, "G": 0.1, "T": 0.0}, # A
        {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}, # N
    ],
    "KOZAK": [  # consensus: AAAAAA
        {"A": 0.9, "C": 0.02, "G": 0.05, "T": 0.03},
        {"A": 0.9, "C": 0.02, "G": 0.05, "T": 0.03},
        {"A": 0.9, "C": 0.02, "G": 0.05, "T": 0.03},
        {"A": 0.9, "C": 0.02, "G": 0.05, "T": 0.03},
        {"A": 0.9, "C": 0.02, "G": 0.05, "T": 0.03},
        {"A": 0.9, "C": 0.02, "G": 0.05, "T": 0.03},
    ]
}

# S. cerevisiae 대표 유전자들의 실험적/문헌적 기본 발현 강도 (0 ~ 100 스케일)
GENE_BASELINES = {
    "YGR192C": 100.0, "TDH3": 100.0,
    "YAL003W": 80.0,  "TEF1": 80.0, "EFB1": 80.0,
    "YCR012W": 70.0,  "PGK1": 70.0,
    "YOL086C": 40.0,  "ADH1": 40.0,
    "YFL039C": 15.0,  "ACT1": 15.0,
    "YJR048W": 5.0,   "CYC1": 5.0,
    
    # 주요 개량 및 대사 관련 유전자 baseline
    "YOL049W": 35.0,  "GSH1": 35.0,
    "YBR029C": 30.0,  "GSH2": 30.0,
    "YEL046C": 25.0,  "GLR1": 25.0,
    "YLR303W": 50.0,  "MET17": 50.0,
    "YDL054C": 8.0,   "MIG1": 8.0,
    "YPL075W": 12.0,  "GCN4": 12.0
}

def get_baseline_strength(gene_symbol, systematic_name, raw_wt_score):
    """Retrieve or estimate the baseline expression strength (0 to 100 scale) of the wild-type promoter."""
    for key in [gene_symbol, systematic_name]:
        if key and key.upper() in GENE_BASELINES:
            return GENE_BASELINES[key.upper()]
            
    # Estimate baseline based on raw score distribution if not found in baseline map
    if raw_wt_score < 650:
        return 12.0   # Weak promoter baseline
    elif raw_wt_score < 850:
        return 40.0   # Medium promoter baseline
    else:
        return 75.0   # Strong promoter baseline


def get_pwm_score(subseq, tf_id):
    """Calculate PWM score for a subsequence. Returns value between 0.0 and 1.0."""
    subseq = subseq.upper()
    if tf_id not in TF_PWMS:
        return 0.0
    
    pwm = TF_PWMS[tf_id]
    cmp_len = min(len(subseq), len(pwm))
    if cmp_len == 0:
        return 0.0
        
    score_sum = 0.0
    max_possible = 0.0
    min_possible = 0.0
    
    for idx in range(cmp_len):
        char = subseq[idx]
        freq = pwm[idx].get(char, 0.25)
        if freq < 0.01:
            freq = 0.01
        
        score_sum += np.log2(freq / 0.25)
        
        max_freq = max(pwm[idx].values())
        max_possible += np.log2(max_freq / 0.25)
        
        min_freq = min(pwm[idx].values())
        if min_freq < 0.01:
            min_freq = 0.01
        min_possible += np.log2(min_freq / 0.25)
        
    if max_possible == min_possible:
        return 1.0
        
    normalized = (score_sum - min_possible) / (max_possible - min_possible)
    return max(0.0, min(1.0, normalized))

def get_tata_distance_multiplier(pos, seq_len):
    """Evaluate TATA Box positioning efficiency relative to TSS (at the downstream end)."""
    dist = seq_len - pos
    # Optimum distance is typically between 30bp and 120bp upstream of TSS in S. cerevisiae
    if 30 <= dist <= 120:
        return 1.0
    elif dist < 30:
        return max(0.1, dist / 30.0)
    else:
        return max(0.1, float(np.exp(-(dist - 120) / 100.0)))

def get_activator_decay(act_pos, tata_pos):
    """Activator effect decays exponentially as its distance from TATA Box/TSS increases."""
    dist = abs(act_pos - tata_pos)
    if dist <= 100:
        return 1.0
    return max(0.05, float(np.exp(-(dist - 100) / 300.0)))

def get_translation_efficiency_mfe(seq):
    """Estimate translation efficiency based on 5' UTR mRNA secondary structure MFE approximation (Strategy 3).
    S. cerevisiae UTR is at the 3' end of the 1000bp promoter sequence (typically final 60bp).
    Stable secondary structures (low MFE due to high GC / self-pairing) drop translation rate.
    """
    seq_len = len(seq)
    # Extract the typical 5' UTR region (final 60 bp of the promoter)
    utr_seq = seq[max(0, seq_len - 60):].upper()
    if not utr_seq:
        return 1.0
        
    # 1. GC content penalty on mRNA secondary structure
    gc_count = utr_seq.count('G') + utr_seq.count('C')
    gc_ratio = gc_count / len(utr_seq)
    
    # Base penalty: high GC in 5' UTR forms rigid hairpin loops easily
    gc_penalty = 0.0
    if gc_ratio > 0.40:
        # Penalize up to -20% translation rate if GC ratio is very high
        gc_penalty = min(0.20, (gc_ratio - 0.40) * 0.8)
        
    # 2. Self-complementarity check (approximate stem-loop formation propensity)
    # We scan for small complementary stems (length 4) within the 60bp UTR
    stem_count = 0
    utr_len = len(utr_seq)
    for i in range(utr_len - 8):
        stem_candidate = utr_seq[i:i+4]
        if 'N' in stem_candidate:
            continue
        # Find reverse complement
        rev_comp = get_reverse_complement(stem_candidate)
        # Search downstream for this reverse complement to see if they form a loop
        if rev_comp in utr_seq[i+8:]:
            # Found a pairing stem-loop target
            stem_count += 1
            
    # Self-complementary loops penalty: max 15% reduction
    loop_penalty = min(0.15, stem_count * 0.03)
    
    # Net Translation Efficiency multiplier (range [0.65, 1.0])
    translation_mult = 1.0 - gc_penalty - loop_penalty
    return max(0.65, translation_mult)

def get_chromatin_accessibility(pos, seq, seq_len):
    """Calculate local accessibility based on Poly-A NDR tracts and GC-dependent bendability (Strategy 1 & 4)."""
    accessibility = 1.0
    
    # 1. Poly-dA:dT NDR (Nucleosome Depleted Region) effect (Strategy 1)
    # Exclude nucleosomes with continuous A/T tracts. Effect decays with distance.
    ndr_bonus = 0.0
    for m in re.finditer(r"[A]{5,}|[T]{5,}", seq.upper()):
        track_len = len(m.group(0))
        # Center position of the poly-A/T tract
        track_center = (m.start() + m.end()) / 2.0
        dist = abs(pos - track_center)
        
        if dist <= 150:
            # Longer tracts (>=7) provide stronger nucleosome displacement than short ones (5-6)
            intensity = 0.4 if track_len >= 7 else 0.15
            # Decay exponentially with distance from tract center
            ndr_bonus += intensity * float(np.exp(-dist / 80.0))
            
    # Maximum accessibility bonus is capped at +40% (1.4 multiplier)
    accessibility += min(0.4, ndr_bonus)
            
    # 2. GC content & bendability penalty (Strategy 4)
    # High or extremely low GC ratio decreases DNA bendability and increases nucleosome occupancy
    w_start = max(0, pos - 30)
    w_end = min(seq_len, pos + 30)
    window = seq[w_start:w_end]
    if len(window) > 0:
        gc_count = window.upper().count('G') + window.upper().count('C')
        gc_ratio = gc_count / len(window)
        
        # S. cerevisiae optimal GC content for bendability is 30% - 42%
        if gc_ratio > 0.42:
            # Linear decay down to 0.7 for GC = 0.65
            gc_mult = max(0.7, 1.0 - (gc_ratio - 0.42) * 1.3)
            accessibility *= gc_mult
        elif gc_ratio < 0.30:
            # Linear decay down to 0.8 for GC = 0.15
            gc_mult = max(0.8, 1.0 - (0.30 - gc_ratio) * 1.3)
            accessibility *= gc_mult
            
    return accessibility

def calculate_absolute_expression_index(seq, sites, reference_tata_pos=None):
    """Compute the absolute transcriptional + translational expression strength index."""
    seq_len = len(seq)
    
    # 1. TATA Box score & distance effect
    tata_sites = [s for s in sites if s["tf_id"] == "TATA" and s["strand"] == "+"]
    if tata_sites:
        best_tata = max(tata_sites, key=lambda x: x["score"])
        tata_pos = best_tata["start"]
        tata_score = best_tata["score"]
        tata_dist_mult = get_tata_distance_multiplier(tata_pos, seq_len)
        tata_eff = tata_score * tata_dist_mult
    else:
        # If TATA-box is mutated/absent, but we have a reference WT TATA position,
        # we preserve that coordinate to avoid spatial warp artifacts in activator decay calculations.
        tata_pos = reference_tata_pos if reference_tata_pos is not None else (seq_len - 120)
        tata_eff = 0.15 # Weak baseline for TATA-less promoters
        
    # 2. Kozak Sequence Translation multiplier
    kozak_sites = [s for s in sites if s["tf_id"] == "KOZAK"]
    kozak_score = kozak_sites[0]["score"] if kozak_sites else 0.5
    
    # 3. Activator and Repressor contributions
    act_sum = 0.0
    rep_sum = 0.0
    
    # Cooperative binding (e.g. RAP1 and GCR1 synergize within 80bp)
    rap1_sites = [s for s in sites if s["tf_id"] == "RAP1"]
    gcr1_sites = [s for s in sites if s["tf_id"] == "GCR1"]
    
    has_synergy = False
    for r in rap1_sites:
        for g in gcr1_sites:
            if abs(r["start"] - g["start"]) <= 80:
                has_synergy = True
                break
        if has_synergy:
            break
            
    synergy_multiplier = 1.5 if has_synergy else 1.0
    
    # To maintain spatial coordinate reference for decay, use the reference TATA pos if available
    decay_anchor_pos = reference_tata_pos if reference_tata_pos is not None else tata_pos
    
    for s in sites:
        if s["tf_id"] in ["TATA", "KOZAK"]:
            continue
            
        acc = get_chromatin_accessibility(s["start"], seq, seq_len)
        decay = get_activator_decay(s["start"], decay_anchor_pos) if s["type"] == "activator" else 1.0
        
        eff_score = s["score"] * acc * decay
        
        if s["type"] == "activator":
            act_sum += s["weight"] * eff_score * synergy_multiplier
        elif s["type"] == "repressor":
            rep_sum += s["weight"] * eff_score
            
    # Integrated Polynomial formula for Transcription Index
    transcription_index = 20.0 + (80.0 * tata_eff) + act_sum - rep_sum
    transcription_index = max(5.0, transcription_index) # Capped minimum
    
    # Translation scaling factor based on Kozak score & 5' UTR MFE (Strategy 3)
    # Range of translation multiplier is [0.5, 1.3]
    kozak_mult = 0.5 + (0.8 * kozak_score)
    utr_mfe_mult = get_translation_efficiency_mfe(seq)
    
    return transcription_index * kozak_mult * utr_mfe_mult

def iupac_to_regex(consensus):
    regex = ""
    # Process IUPAC codes, check for potential N{count} style loops
    i = 0
    while i < len(consensus):
        char = consensus[i]
        regex += IUPAC_MAP.get(char, char)
        i += 1
    return regex

def get_match_score(subseq, consensus):
    """Calculate matching score of a subsequence against a consensus string (0.0 to 1.0)."""
    if len(subseq) != len(consensus):
        return 0.0
    
    matches = 0
    for s_char, c_char in zip(subseq, consensus):
        pattern = IUPAC_MAP.get(c_char, c_char)
        if re.match(f"^{pattern}$", s_char):
            matches += 1
            
    return matches / len(consensus)

def scan_promoter_motifs(seq):
    """Scan a promoter sequence for all TF motifs. Returns a list of found sites."""
    found_sites = []
    
    for tf_key, info in TF_MOTIFS.items():
        motif_len = len(info["consensus"])
        regex_pattern = iupac_to_regex(info["consensus"])
        threshold = info.get("threshold", 0.75)
        
        # Scan forward strand
        for i in range(len(seq) - motif_len + 1):
            subseq = seq[i:i+motif_len]
            score = get_match_score(subseq, info["consensus"])
            
            # Threshold for binding
            if score >= threshold:
                # Add to results
                found_sites.append({
                    "tf_id": tf_key,
                    "tf_name": info["name"],
                    "start": i,
                    "end": i + motif_len,
                    "strand": "+",
                    "sequence": subseq,
                    "consensus": info["consensus"],
                    "score": round(get_pwm_score(subseq, tf_key), 3),
                    "weight": info["weight"],
                    "type": info["type"],
                    "desc": info["desc"]
                })
                
        # Scan reverse strand (except for palindromic or symmetric cases if redundant, but standard is to check both)
        # To avoid overly complicating coordinates, we map reverse complement back to forward coordinates.
        rev_seq = get_reverse_complement(seq)
        for i in range(len(rev_seq) - motif_len + 1):
            subseq = rev_seq[i:i+motif_len]
            score = get_match_score(subseq, info["consensus"])
            
            if score >= threshold:
                # Forward coordinates mapping
                # i on reverse strand corresponds to:
                fwd_start = len(seq) - (i + motif_len)
                fwd_end = len(seq) - i
                
                # Check for redundancy (don't add duplicates at same location if palindromic/symmetric score matches)
                is_duplicate = False
                for s in found_sites:
                    if s["tf_id"] == tf_key and s["start"] == fwd_start:
                        is_duplicate = True
                        break
                        
                if not is_duplicate:
                    found_sites.append({
                        "tf_id": tf_key,
                        "tf_name": info["name"],
                        "start": fwd_start,
                        "end": fwd_end,
                        "strand": "-",
                        "sequence": get_reverse_complement(subseq),
                        "consensus": info["consensus"],
                        "score": round(get_pwm_score(subseq, tf_key), 3),
                        "weight": info["weight"],
                        "type": info["type"],
                        "desc": info["desc"]
                    })
                    
    # Append Kozak sequence at the very end (-6bp to -1bp)
    kozak_subseq = seq[-6:]
    kozak_score = get_pwm_score(kozak_subseq.upper(), "KOZAK")
    found_sites.append({
        "tf_id": "KOZAK",
        "tf_name": "Kozak Sequence",
        "start": len(seq) - 6,
        "end": len(seq),
        "strand": "+",
        "sequence": kozak_subseq,
        "consensus": "AAAAAA",
        "score": round(kozak_score, 3),
        "weight": 40,
        "type": "activator",
        "desc": "번역 개시 효율 조절 인자. 최적의 Kozak 서열(AAAAAA)은 단백질 번역 수준을 극대화합니다."
    })

    # Sort sites by start position
    found_sites.sort(key=lambda x: x["start"])
    return found_sites

def get_reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def apply_mutations(base_seq, mutations):
    """Apply a list of mutations to base_seq. mutations: [{'pos': int, 'to': str}]"""
    seq_list = list(base_seq)
    for mut in mutations:
        pos = mut['pos']
        to_base = mut['to'].upper()
        if 0 <= pos < len(seq_list):
            seq_list[pos] = to_base
    return "".join(seq_list)

def predict_expression_change(wild_seq, mutant_seq, wild_sites=None, gene_symbol=None, systematic_name=None):
    """Predict the expression change (%) of a mutant promoter relative to WT, and calibrate absolute strengths."""
    if wild_sites is None:
        wild_sites = scan_promoter_motifs(wild_seq)
        
    mut_sites = scan_promoter_motifs(mutant_seq)
    
    # 1. Identify WT TATA box position to keep coordinate reference constant
    wt_tatas = [s for s in wild_sites if s["tf_id"] == "TATA" and s["strand"] == "+"]
    wt_tata_pos = max(wt_tatas, key=lambda x: x["score"])["start"] if wt_tatas else None
    
    wt_index = calculate_absolute_expression_index(wild_seq, wild_sites, reference_tata_pos=wt_tata_pos)
    mut_index = calculate_absolute_expression_index(mutant_seq, mut_sites, reference_tata_pos=wt_tata_pos)
    
    # 2. Check if TATA box was present in WT but got completely destroyed/mutated in mutant
    mut_tatas = [s for s in mut_sites if s["tf_id"] == "TATA" and s["strand"] == "+"]
    tata_destroyed = len(wt_tatas) > 0 and len(mut_tatas) == 0
    
    # 3. Apply a biological penalty if a functional TATA box is completely knocked out.
    # In S. cerevisiae TATA-containing promoters, destroying the TATA box typically results
    # in a massive (70-90%) drop in transcription initiation efficiency.
    if tata_destroyed:
        mut_index = min(mut_index, wt_index * 0.3)
        
    # 3.1 Epistasis Nonlinear Penalty for multi-mutations (Strategy 2)
    # Count the total number of single nucleotide variants (SNVs) between WT and Mut
    snv_count = sum(1 for w, m in zip(wild_seq, mutant_seq) if w.upper() != m.upper())
    if snv_count >= 2:
        # Apply exponential epistasis decay on the mutant score
        # e.g., 2 mutations -> -3% extra penalty, 5 mutations -> -11%, 10 mutations -> -24%
        epistasis_mult = float(np.exp(-0.03 * (snv_count - 1)))
        epistasis_mult = max(0.70, epistasis_mult)  # Cap epistasis penalty at 30% reduction
        mut_index *= epistasis_mult
        
    relative_ratio = mut_index / wt_index
    raw_final_expr = relative_ratio * 100.0
    
    # 3.2 Ensure exactly 100% when there is no sequence difference
    if snv_count == 0 or abs(raw_final_expr - 100.0) < 1e-5:
        final_expr = 100.0
    else:
        # Soft-classification/Calibration Curve to model biological saturation and noise (Strategy 5)
        # Instead of hard clipping, use sigmoid-like soft saturation for extreme changes
        if raw_final_expr > 100.0:
            # Soft-cap hyper-activation at 350.0% using a logarithmic/tanh compression
            final_expr = 100.0 + 250.0 * float(np.tanh((raw_final_expr - 100.0) / 250.0))
        else:
            # Soft-cap severe inactivation down to 5.0% with zero-point alignment at 100%
            final_expr = 100.0 + 95.0 * float(np.tanh((raw_final_expr - 100.0) / 50.0))
    
    # 캘리브레이션 연산 적용 (Yeast 표준 및 타겟 유전자 baseline 대조)
    baseline = get_baseline_strength(gene_symbol, systematic_name, wt_index)
    calibrated_wt = baseline
    calibrated_mut = baseline * relative_ratio
    
    # Calibrated score capping: 0 to 150
    calibrated_wt = max(0.0, min(150.0, calibrated_wt))
    calibrated_mut = max(0.0, min(150.0, calibrated_mut))
    
    return {
        "predicted_value": round(final_expr, 1),
        "change_percentage": round(final_expr - 100.0, 1),
        "tata_destroyed": tata_destroyed,
        "calibrated_wt": round(calibrated_wt, 2),
        "calibrated_mut": round(calibrated_mut, 2)
    }

def get_mutation_recommendations(seq):
    """Generate potential point mutations that increase or decrease promoter expression."""
    sites = scan_promoter_motifs(seq)
    
    up_recommendations = []
    down_recommendations = []
    
    # 1. UP recommendations
    # Idea A: Look for repressor sites (MIG1) and recommend breaking them
    mig1_sites = [s for s in sites if s["tf_id"] == "MIG1"]
    for site in mig1_sites:
        break_pos = site["start"] + 2 # pos of G
        orig_base = seq[break_pos]
        up_recommendations.append({
            "pos": break_pos,
            "from": orig_base,
            "to": "T" if orig_base != "T" else "A",
            "effect": "MIG1 Repressor Site 파괴",
            "desc": "Glucose Repression(포도당 발현 억제)을 해제하여 배지 내 포도당이 있는 조건 하에서 발현량을 향상시킵니다.",
            "impact_type": "high_increase"
        })
        
    # Idea B: Look for weak activator sites (score < 1.0) and propose making them consensus
    weak_activators = [s for s in sites if s["type"] == "activator" and s["score"] < 1.0]
    for site in weak_activators:
        consensus = site["consensus"]
        for idx, (s_char, c_char) in enumerate(zip(site["sequence"], consensus)):
            pattern = IUPAC_MAP.get(c_char, c_char)
            if not re.match(f"^{pattern}$", s_char):
                suggested_base = c_char
                if c_char in IUPAC_MAP and c_char not in ['A', 'C', 'G', 'T']:
                    if c_char == 'W': suggested_base = 'A'
                    elif c_char == 'S': suggested_base = 'G'
                    elif c_char == 'R': suggested_base = 'G'
                    elif c_char == 'Y': suggested_base = 'C'
                    elif c_char == 'M': suggested_base = 'A'
                    elif c_char == 'K': suggested_base = 'G'
                    else: suggested_base = 'A'
                    
                up_recommendations.append({
                    "pos": site["start"] + idx,
                    "from": s_char,
                    "to": suggested_base,
                    "effect": f"{site['tf_name']} Activator Site 최적화",
                    "desc": f"전사인자 {site['tf_name']}와의 결합력을 극대화하여 발현 활성을 높입니다.",
                    "impact_type": "increase"
                })
                break
                
    # Idea C: Propose creating TATA box if none exists or if it is weak
    tata_sites = [s for s in sites if s["tf_id"] == "TATA"]
    if not tata_sites:
        # TSS region insertion (around -120bp)
        tata_pos = len(seq) - 120
        if 0 <= tata_pos < len(seq) - 8:
            up_recommendations.append({
                "pos": tata_pos,
                "from": seq[tata_pos:tata_pos+8],
                "to": "TATATAAA",
                "effect": "TATA Box 신규 도입",
                "desc": "핵심 프로모터 영역에 TATA Box를 생성하여 전사 개시 효율을 비약적으로 상승시킵니다.",
                "is_multi_base": True,
                "impact_type": "high_increase"
            })
    else:
        best_tata = max(tata_sites, key=lambda x: x["score"])
        if best_tata["score"] < 1.0:
            up_recommendations.append({
                "pos": best_tata["start"],
                "from": best_tata["sequence"],
                "to": "TATATAAA",
                "effect": "TATA Box Consensus화",
                "desc": "TATA Box 서열을 완벽한 Consensus(TATATAAA)로 전환하여 RNA 중합효소 복합체 유인력을 높입니다.",
                "is_multi_base": True,
                "impact_type": "increase"
            })

    # Idea D (NEW): De novo introduction of strong constitutive activators (RAP1, GCR1, GCN4)
    # Check if we already have RAP1 site in the promoter
    rap1_present = any(s["tf_id"] == "RAP1" for s in sites)
    if not rap1_present:
        # Propose inserting RAP1 at upstream region (around -450bp)
        rap1_pos = len(seq) - 450
        if 0 <= rap1_pos < len(seq) - 13:
            up_recommendations.append({
                "pos": rap1_pos,
                "from": seq[rap1_pos:rap1_pos+13],
                "to": "AAACCCAGCATTG", # RAP1 consensus match
                "effect": "RAP1 Activator Site 신규 도입",
                "desc": "구성적 강력 활성인자인 RAP1 결합 사이트를 상류에 신규로 도입하여 전체 전사 속도를 대폭 증대시킵니다.",
                "is_multi_base": True,
                "impact_type": "high_increase"
            })

    gcr1_present = any(s["tf_id"] == "GCR1" for s in sites)
    if not gcr1_present:
        # Propose inserting GCR1 at mid-stream region (around -350bp)
        gcr1_pos = len(seq) - 350
        if 0 <= gcr1_pos < len(seq) - 5:
            up_recommendations.append({
                "pos": gcr1_pos,
                "from": seq[gcr1_pos:gcr1_pos+5],
                "to": "CTTCC", # GCR1 consensus match
                "effect": "GCR1 Activator Site 신규 도입",
                "desc": "글리콜리시스 유전자들의 활성인자인 GCR1 사이트를 신규 도입하여 프로모터 기본 활성을 끌어올립니다.",
                "is_multi_base": True,
                "impact_type": "increase"
            })

    gcn4_present = any(s["tf_id"] == "GCN4" for s in sites)
    if not gcn4_present:
        # Propose inserting GCN4 (around -250bp)
        gcn4_pos = len(seq) - 250
        if 0 <= gcn4_pos < len(seq) - 6:
            up_recommendations.append({
                "pos": gcn4_pos,
                "from": seq[gcn4_pos:gcn4_pos+6],
                "to": "TGACTC", # GCN4 consensus match
                "effect": "GCN4 Activator Site 신규 도입",
                "desc": "아미노산 starvation 조건 및 스트레스 조건 하에서 발현을 촉진하는 GCN4 결합 모티프를 새로 삽입합니다.",
                "is_multi_base": True,
                "impact_type": "increase"
            })

    # Idea E (NEW): Kozak Sequence Optimization
    kozak_sites = [s for s in sites if s["tf_id"] == "KOZAK"]
    if kozak_sites:
        kozak_site = kozak_sites[0]
        if kozak_site["score"] < 1.0:
            up_recommendations.append({
                "pos": kozak_site["start"],
                "from": kozak_site["sequence"],
                "to": "AAAAAA",
                "effect": "Kozak Sequence 최적화",
                "desc": "ATG 직전의 Kozak 서열을 최적의 컨센서스(AAAAAA)로 보완하여 단백질 번역(Translation) 효율을 대폭 향상시킵니다.",
                "is_multi_base": True,
                "impact_type": "increase"
            })
            
    # 2. DOWN recommendations
    # Idea A: Destroy TATA box
    for site in tata_sites:
        # TSS로부터 발현 기여도가 있는 유효한 TATA Box만 무력화 추천 대상으로 삼습니다.
        if get_tata_distance_multiplier(site["start"], len(seq)) > 0.15:
            down_recommendations.append({
                "pos": site["start"] + 2, # Mutate third base of TATA (T->G)
                "from": seq[site["start"] + 2],
                "to": "G",
                "effect": "TATA Box 활성 무력화",
                "desc": "TATA Box 서열을 변형시켜 RNA 중합효소 결합을 원천 차단하고 프로모터 발현 강도를 급격히 억제합니다.",
                "impact_type": "high_decrease"
            })
        
    # Idea B: Destroy Activator sites (RAP1, GCR1, GAL4, GCN4 etc.)
    strong_activators = [s for s in sites if s["type"] == "activator" and s["tf_id"] != "TATA"]
    for site in strong_activators:
        break_pos = site["start"] + (len(site["sequence"]) // 2)
        orig_base = seq[break_pos]
        down_recommendations.append({
            "pos": break_pos,
            "from": orig_base,
            "to": "G" if orig_base != "G" else "A",
            "effect": f"{site['tf_name']} Activator Site 파괴",
            "desc": f"전사인자 {site['tf_name']}의 활성화 결합을 방해하여 프로모터 활성을 하향시킵니다.",
            "impact_type": "decrease"
        })

    # Idea C (NEW): De novo introduction of MIG1 Repressor Site (Glucose Repression)
    mig1_present = any(s["tf_id"] == "MIG1" for s in sites)
    if not mig1_present:
        # Propose inserting MIG1 repressor site near core promoter (around -180bp)
        mig1_insert_pos = len(seq) - 180
        if 0 <= mig1_insert_pos < len(seq) - 6:
            down_recommendations.append({
                "pos": mig1_insert_pos,
                "from": seq[mig1_insert_pos:mig1_insert_pos+6],
                "to": "GCGGGG", # MIG1 consensus match (SYGGRG)
                "effect": "MIG1 Repressor Site 신규 도입",
                "desc": "TATA box 상류에 MIG1 결합 서열을 생성하여 포도당 배지 상에서 유전자의 전사 수준을 능동적으로 억제시킵니다.",
                "is_multi_base": True,
                "impact_type": "high_decrease"
            })

    # Idea D (NEW): Kozak Sequence Weakening
    kozak_sites = [s for s in sites if s["tf_id"] == "KOZAK"]
    if kozak_sites:
        kozak_site = kozak_sites[0]
        if kozak_site["score"] >= 0.5:
            down_recommendations.append({
                "pos": kozak_site["start"],
                "from": kozak_site["sequence"],
                "to": "CGTGGG",
                "effect": "Kozak Sequence 약화",
                "desc": "번역 개시 신호 서열을 무력화하여 리보솜 결합 능력을 방해하고 단백질 발현량을 급격히 감쇄시킵니다.",
                "is_multi_base": True,
                "impact_type": "decrease"
            })

    # Sort and slice recommendations to top 5
    return {
        "up": up_recommendations[:5],
        "down": down_recommendations[:5]
    }
