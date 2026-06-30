import os
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel
from typing import List, Dict, Any

import data_loader
import analyzer

app = FastAPI(title="Yeast Promoter Modifier API")

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Request schemas
class Mutation(BaseModel):
    pos: int
    to: str

class PredictionRequest(BaseModel):
    wild_seq: str
    mutations: List[Mutation]

class RecommendationRequest(BaseModel):
    seq: str

@app.get("/api/search")
def search_genes(q: str = Query(..., min_length=1)):
    """Search for Yeast genes matching standard symbol, systematic name, or alias."""
    q = q.strip().upper()
    gene_map = data_loader.get_gene_mapping()
    
    results = []
    seen_systematic = set()
    
    # 1. Direct match first
    if q in gene_map:
        info = gene_map[q]
        sys_name = info["systematic_name"]
        seen_systematic.add(sys_name)
        results.append({
            "systematic_name": sys_name,
            "symbol": info["symbol"],
            "description": info["description"]
        })
        
    # 2. Substring match to populate autocomplete/search results
    for key, info in gene_map.items():
        sys_name = info["systematic_name"]
        if sys_name in seen_systematic:
            continue
            
        # Match systematic name, standard symbol, or aliases
        if (q in sys_name or 
            q in info["symbol"].upper() or 
            any(q in alias.upper() for alias in info["aliases"])):
            
            seen_systematic.add(sys_name)
            results.append({
                "systematic_name": sys_name,
                "symbol": info["symbol"],
                "description": info["description"]
            })
            
            if len(results) >= 20: # Limit autocomplete results
                break
                
    return {"results": results}

@app.get("/api/promoter")
def get_promoter(gene: str):
    """Retrieve promoter sequence and scan its core structural features."""
    promoter_data = data_loader.get_promoter_data(gene)
    if not promoter_data:
        raise HTTPException(status_code=404, detail="유전자를 찾을 수 없거나 프로모터 서열을 가져오는데 실패했습니다.")
        
    # Scan features
    seq = promoter_data["seq"]
    sites = analyzer.scan_promoter_motifs(seq)
    
    return {
        "systematic_name": promoter_data["systematic_name"],
        "symbol": promoter_data["symbol"],
        "description": promoter_data["description"],
        "seq": seq,
        "sites": sites
    }

@app.post("/api/predict")
def predict_mutation(payload: PredictionRequest):
    """Predict expression levels after applying point mutations."""
    wild_seq = payload.wild_seq.upper()
    muts_dict = [{"pos": m.pos, "to": m.to} for m in payload.mutations]
    
    # Apply mutations to get mutant sequence
    mutant_seq = analyzer.apply_mutations(wild_seq, muts_dict)
    
    # Analyze WT and Mutant
    wt_sites = analyzer.scan_promoter_motifs(wild_seq)
    mutant_sites = analyzer.scan_promoter_motifs(mutant_seq)
    
    # Run prediction
    prediction = analyzer.predict_expression_change(wild_seq, mutant_seq, wild_sites=wt_sites)
    
    return {
        "mutant_seq": mutant_seq,
        "predicted_value": prediction["predicted_value"],
        "change_percentage": prediction["change_percentage"],
        "tata_destroyed": prediction["tata_destroyed"],
        "mutant_sites": mutant_sites
    }

@app.post("/api/recommend")
def recommend_mutations(payload: RecommendationRequest):
    """Recommend point mutations to tune expression levels up or down."""
    seq = payload.seq.upper()
    recs = analyzer.get_mutation_recommendations(seq)
    return recs

# Serve Web UI
static_path = os.path.join(os.path.dirname(__file__), "static")
os.makedirs(static_path, exist_ok=True)

app.mount("/static", StaticFiles(directory=static_path), name="static")

@app.get("/")
def read_index():
    index_file = os.path.join(static_path, "index.html")
    if os.path.exists(index_file):
        return FileResponse(index_file)
    return JSONResponse(content={"message": "Yeast Promoter Modifier API is running. UI is not yet built."}, status_code=200)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
