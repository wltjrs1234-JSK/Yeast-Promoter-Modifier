// State Management
let currentPromoter = null; // Stores { systematic_name, symbol, description, seq, sites }
let activeMutations = [];   // Stores [{ pos, from, to }]
let recommendedMutations = { up: [], down: [] };
let activeRecTab = 'up';    // 'up' or 'down'
let selectedBaseIndex = null; // For nucleotide mutator selector positioning
let zoomLevel = 1.0;
let panOffset = 0.0;

// DOM Elements
const searchInput = document.getElementById('gene-search-input');
const searchBtn = document.getElementById('search-btn');
const autocompleteList = document.getElementById('autocomplete-list');
const geneInfoBox = document.getElementById('selected-gene-info');
const infoSymbol = document.getElementById('info-symbol');
const infoSysname = document.getElementById('info-sysname');
const infoDesc = document.getElementById('info-desc');

const trackSvg = document.getElementById('promoter-track-svg');
const trackTooltip = document.getElementById('track-tooltip');
const btnZoomIn = document.getElementById('btn-zoom-in');
const btnZoomOut = document.getElementById('btn-zoom-out');
const btnZoomReset = document.getElementById('btn-zoom-reset');

const sequenceGrid = document.getElementById('sequence-grid');
const mutatorSelector = document.getElementById('mutator-selector');
const closeMutatorBtn = document.getElementById('close-mutator-btn');
const selectedBaseCoord = document.getElementById('selected-base-coord');
const resetMutationsBtn = document.getElementById('reset-mutations-btn');

const predictionVal = document.getElementById('prediction-val');
const gaugeFillArc = document.getElementById('gauge-fill-arc');
const statusTag = document.getElementById('status-tag');
const statusDesc = document.getElementById('status-desc');
const predictionStatusBox = document.getElementById('prediction-status');

const activeMutationsList = document.getElementById('active-mutations-list');
const mutationsCount = document.getElementById('mutations-count');
const clearAllMutationsBtn = document.getElementById('clear-all-mutations-btn');

const tabUpBtn = document.getElementById('tab-up-btn');
const tabDownBtn = document.getElementById('tab-down-btn');
const recommendationsList = document.getElementById('recommendations-list');

// Helper to update CSS pointer based on zoomLevel
function updateTrackCursor() {
    const trackWrapper = document.querySelector('.map-track-wrapper');
    if (trackWrapper) {
        trackWrapper.style.cursor = zoomLevel > 1.0 ? 'grab' : 'default';
    }
}

// Initialize Events
document.addEventListener('DOMContentLoaded', () => {
    setupAutocomplete();
    setupMutator();
    setupRecommendationTabs();
    setupTrackZoomPan();
    
    if (btnZoomIn) {
        btnZoomIn.addEventListener('click', () => {
            if (!currentPromoter) return;
            zoomLevel = Math.min(15.0, zoomLevel * 1.3);
            updateTrackCursor();
            updateTrackViewBox();
        });
    }
    if (btnZoomOut) {
        btnZoomOut.addEventListener('click', () => {
            if (!currentPromoter) return;
            zoomLevel = Math.max(1.0, zoomLevel / 1.3);
            updateTrackCursor();
            updateTrackViewBox();
        });
    }
    if (btnZoomReset) {
        btnZoomReset.addEventListener('click', () => {
            zoomLevel = 1.0;
            panOffset = 0.0;
            updateTrackCursor();
            updateTrackViewBox();
        });
    }
    
    resetMutationsBtn.addEventListener('click', resetAllMutations);
    if (clearAllMutationsBtn) {
        clearAllMutationsBtn.addEventListener('click', resetAllMutations);
    }
    searchBtn.addEventListener('click', () => triggerSearch(searchInput.value));
    
    // Load default TDH3 promoter on start
    loadPromoter('TDH3');
});

// Setup Gene Search & Autocomplete dropdown
function setupAutocomplete() {
    let debounceTimer;
    
    searchInput.addEventListener('input', () => {
        clearTimeout(debounceTimer);
        const query = searchInput.value.trim();
        
        if (query.length < 1) {
            autocompleteList.classList.add('hidden');
            return;
        }
        
        debounceTimer = setTimeout(async () => {
            try {
                const res = await fetch(`/api/search?q=${encodeURIComponent(query)}`);
                const data = await res.json();
                renderAutocomplete(data.results);
            } catch (err) {
                console.error("Error searching genes:", err);
            }
        }, 250);
    });
    
    // Hide dropdown when clicking outside
    document.addEventListener('click', (e) => {
        if (!searchInput.contains(e.target) && !autocompleteList.contains(e.target)) {
            autocompleteList.classList.add('hidden');
        }
    });
}

function renderAutocomplete(results) {
    if (!results || results.length === 0) {
        autocompleteList.classList.add('hidden');
        return;
    }
    
    autocompleteList.innerHTML = '';
    results.forEach(gene => {
        const div = document.createElement('div');
        div.className = 'autocomplete-item';
        div.innerHTML = `
            <span class="symbol">${gene.symbol}</span>
            <span class="sysname">${gene.systematic_name}</span>
        `;
        div.addEventListener('click', () => {
            searchInput.value = gene.symbol;
            autocompleteList.classList.add('hidden');
            loadPromoter(gene.systematic_name);
        });
        autocompleteList.appendChild(div);
    });
    
    autocompleteList.classList.remove('hidden');
}

async function triggerSearch(query) {
    if (!query.trim()) return;
    try {
        const res = await fetch(`/api/search?q=${encodeURIComponent(query)}`);
        const data = await res.json();
        if (data.results && data.results.length > 0) {
            loadPromoter(data.results[0].systematic_name);
            autocompleteList.classList.add('hidden');
        } else {
            alert("유전자를 찾을 수 없습니다.");
        }
    } catch (err) {
        console.error("Search failed:", err);
    }
}

// Load promoter sequence and structures
async function loadPromoter(geneId) {
    try {
        // Show loading state in sequence editor
        sequenceGrid.innerHTML = `
            <div class="empty-state">
                <i class="fa-solid fa-spinner fa-spin empty-icon"></i>
                <p>프로모터 서열을 분석하는 중...</p>
            </div>
        `;
        
        const res = await fetch(`/api/promoter?gene=${encodeURIComponent(geneId)}`);
        if (!res.ok) throw new Error("API request failed");
        
        currentPromoter = await res.json();
        activeMutations = []; // Clear mutations on new gene load
        zoomLevel = 1.0;      // Reset zoom
        panOffset = 0.0;      // Reset pan
        
        // Clear highlights
        document.querySelectorAll('.highlighted-site').forEach(el => el.classList.remove('highlighted-site'));
        
        // Update Info Box
        infoSymbol.textContent = currentPromoter.symbol;
        infoSysname.textContent = currentPromoter.systematic_name;
        infoDesc.textContent = currentPromoter.description;
        geneInfoBox.classList.remove('hidden');
        
        // Render UI
        renderPromoterTrack(currentPromoter.sites);
        renderSequenceGrid();
        
        // Render initial calibration spectrum before mutation calculation
        renderSpectrum(currentPromoter.calibrated_wt, currentPromoter.calibrated_wt, currentPromoter.references);
        
        updatePrediction();
        fetchRecommendations();
        
    } catch (err) {
        console.error("Error loading promoter:", err);
        sequenceGrid.innerHTML = `
            <div class="empty-state">
                <i class="fa-solid fa-triangle-exclamation empty-icon text-red"></i>
                <p>프로모터 데이터를 가져오는데 실패했습니다. 유전자 이름을 다시 확인하세요.</p>
            </div>
        `;
    }
}

// Render Promoter 1D track map (SVG) with Zoom & Pan support
function renderPromoterTrack(sites) {
    trackSvg.innerHTML = '';
    const seqLen = currentPromoter ? currentPromoter.seq.length : 1000;
    
    // Draw baseline
    const baseline = document.createElementNS("http://www.w3.org/2000/svg", "line");
    baseline.setAttribute("x1", "0");
    baseline.setAttribute("y1", "35");
    baseline.setAttribute("x2", seqLen.toString());
    baseline.setAttribute("y2", "35");
    baseline.setAttribute("class", "track-base-line");
    trackSvg.appendChild(baseline);
    
    // Draw scale ticks and numbers (dynamic loop based on actual length)
    // Decreased spacing to 50bp to have detailed grid when zoomed in, but labels every 100bp
    for (let i = 0; i <= seqLen; i += 50) {
        const tick = document.createElementNS("http://www.w3.org/2000/svg", "line");
        tick.setAttribute("x1", i.toString());
        tick.setAttribute("y1", "30");
        tick.setAttribute("x2", i.toString());
        tick.setAttribute("y2", "40");
        tick.setAttribute("class", "track-tick");
        trackSvg.appendChild(tick);
        
        if (i % 100 === 0) {
            const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
            text.setAttribute("x", i.toString());
            text.setAttribute("y", "58");
            text.setAttribute("text-anchor", "middle");
            text.setAttribute("class", "track-tick-text");
            
            // Map 0 -> -seqLen bp, seqLen -> -1bp
            const coord = i - seqLen;
            text.textContent = coord === 0 ? "-1 bp" : `${coord} bp`;
            trackSvg.appendChild(text);
        }
    }
    
    // Draw sites (TATA, Activators, Repressors)
    sites.forEach(site => {
        const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
        rect.setAttribute("x", site.start.toString());
        rect.setAttribute("y", "20");
        rect.setAttribute("width", (site.end - site.start).toString());
        rect.setAttribute("height", "30");
        rect.setAttribute("rx", "4");
        rect.setAttribute("ry", "4");
        
        let blockClass = "block-activator";
        if (site.tf_id === "TATA") blockClass = "block-tata";
        else if (site.tf_id === "KOZAK") blockClass = "block-kozak";
        else if (site.type === "repressor") blockClass = "block-repressor";
        
        rect.setAttribute("class", `track-block ${blockClass}`);
        
        // Add events for floating tooltip
        rect.addEventListener('mouseenter', (e) => {
            const rectBounds = trackSvg.getBoundingClientRect();
            const xPos = e.clientX - rectBounds.left + 10;
            const yPos = e.clientY - rectBounds.top - 70;
            
            showTooltip(site, xPos, yPos);
        });
        
        rect.addEventListener('mouseleave', hideTooltip);
        
        rect.addEventListener('click', () => {
            // Scroll to the respective sequence position
            const targetEl = document.getElementById(`base-${site.start}`);
            if (targetEl) {
                targetEl.scrollIntoView({ behavior: 'smooth', block: 'center' });
                
                // Clear any previous highlights
                document.querySelectorAll('.highlighted-site').forEach(el => {
                    el.classList.remove('highlighted-site');
                });
                
                // Highlight all bases in the site range (persistent until cleared)
                for (let pos = site.start; pos < site.end; pos++) {
                    const el = document.getElementById(`base-${pos}`);
                    if (el) {
                        el.classList.add('highlighted-site');
                    }
                }
            }
        });
        
        trackSvg.appendChild(rect);
    });
    
    // Apply ViewBox scale/zoom initially
    updateTrackViewBox();
}

function updateTrackViewBox() {
    const seqLen = currentPromoter ? currentPromoter.seq.length : 1000;
    const viewWidth = seqLen / zoomLevel;
    const maxOffset = seqLen - viewWidth;
    panOffset = Math.max(0, Math.min(panOffset, maxOffset));
    
    trackSvg.setAttribute("viewBox", `${panOffset} 0 ${viewWidth} 80`);
}

function showTooltip(site, x, y) {
    let typeText = "Activator Binding Site";
    let typeClass = "indicator-boost";
    
    if (site.tf_id === "TATA") {
        typeText = "Core TATA Box";
        typeClass = "indicator-normal";
    } else if (site.type === "repressor") {
        typeText = "Repressor Binding Site";
        typeClass = "indicator-drop";
    }
    
    trackTooltip.innerHTML = `
        <div class="tooltip-title">${site.tf_name}</div>
        <span class="tooltip-type ${typeClass}">${typeText}</span>
        <div class="tooltip-desc">
            <strong>위치:</strong> -${1000 - site.start} bp ~ -${1000 - site.end} bp (${site.strand} 가닥)<br>
            <strong>서열:</strong> <span style="font-family: monospace; font-weight:bold;">${site.sequence}</span><br>
            <strong>점수:</strong> ${site.score} (Consensus: ${site.consensus})<br>
            <span style="font-size:0.7rem; color:var(--text-muted); display:block; margin-top:4px;">${site.desc}</span>
        </div>
    `;
    
    trackTooltip.style.left = `${x}px`;
    trackTooltip.style.top = `${y}px`;
    trackTooltip.classList.remove('hidden');
}

function hideTooltip() {
    trackTooltip.classList.add('hidden');
}

// Render the 1D grid layout of nucleotides
function renderSequenceGrid() {
    sequenceGrid.innerHTML = '';
    const seq = currentPromoter.seq;
    const sites = currentPromoter.sites;
    
    // Create base annotation mapping
    // We tag each position index (0-999) with its TF affiliation
    const posMeta = new Array(seq.length).fill(null);
    sites.forEach(site => {
        for (let i = site.start; i < site.end; i++) {
            // Priority: TATA > Repressor > Activator
            if (site.tf_id === "TATA") {
                posMeta[i] = { type: 'tata', name: site.tf_name };
            } else if (site.type === 'repressor' && (!posMeta[i] || posMeta[i].type !== 'tata')) {
                posMeta[i] = { type: 'repressor', name: site.tf_name };
            } else if (!posMeta[i]) {
                posMeta[i] = { type: 'activator', name: site.tf_name };
            }
        }
    });
    
    // Build rows of 50 bases
    const rowSize = 50;
    for (let rowIndex = 0; rowIndex < seq.length; rowIndex += rowSize) {
        const rowDiv = document.createElement('div');
        rowDiv.className = 'seq-row';
        
        // Add coordinate label
        const labelDiv = document.createElement('div');
        labelDiv.className = 'row-coord';
        labelDiv.textContent = `-${1000 - rowIndex}`;
        rowDiv.appendChild(labelDiv);
        
        // Add bases
        for (let baseIndex = rowIndex; baseIndex < rowIndex + rowSize && baseIndex < seq.length; baseIndex++) {
            const baseChar = seq[baseIndex];
            const span = document.createElement('div');
            span.className = 'base-block';
            span.id = `base-${baseIndex}`;
            span.dataset.index = baseIndex;
            span.dataset.base = baseChar;
            span.textContent = baseChar;
            
            // Apply TFBS class if applicable
            const meta = posMeta[baseIndex];
            if (meta) {
                if (meta.type === 'tata') span.classList.add('base-tata');
                else if (meta.type === 'repressor') span.classList.add('base-repressor');
                else if (meta.type === 'activator') span.classList.add('base-activator');
                
                span.title = meta.name;
            }
            
            // Check if mutated
            const mut = activeMutations.find(m => m.pos === baseIndex);
            if (mut) {
                span.classList.add('mutated');
                span.textContent = mut.to;
                span.dataset.base = mut.to;
            }
            
            // Event
            span.addEventListener('click', (e) => openMutator(e, baseIndex, baseChar));
            rowDiv.appendChild(span);
        }
        
        sequenceGrid.appendChild(rowDiv);
    }
}

// Nucleotide Mutator Selector popup handling
function setupMutator() {
    closeMutatorBtn.addEventListener('click', () => {
        mutatorSelector.classList.add('hidden');
        selectedBaseCoord.textContent = '-';
    });
    
    const mutButtons = mutatorSelector.querySelectorAll('.mut-btn:not(.btn-cancel)');
    mutButtons.forEach(btn => {
        btn.addEventListener('click', () => {
            const targetBase = btn.dataset.base;
            applyPointMutation(selectedBaseIndex, targetBase);
        });
    });
}

function openMutator(e, index, currentBase) {
    selectedBaseIndex = index;
    const coordVal = 1000 - index;
    selectedBaseCoord.textContent = `-${coordVal} bp (야생형: ${currentBase})`;
    
    // Position mutator selector popup
    const gridBounds = sequenceGrid.getBoundingClientRect();
    const elemBounds = e.target.getBoundingClientRect();
    
    // Relative coordinates
    const topOffset = elemBounds.bottom - gridBounds.top + sequenceGrid.scrollTop + 8;
    const leftOffset = elemBounds.left - gridBounds.left - 40;
    
    mutatorSelector.style.top = `${topOffset}px`;
    mutatorSelector.style.left = `${leftOffset}px`;
    mutatorSelector.classList.remove('hidden');
}

// Apply Mutation
function applyPointMutation(index, toBase) {
    const wildBase = currentPromoter.seq[index];
    
    // Remove mutator popup
    mutatorSelector.classList.add('hidden');
    selectedBaseCoord.textContent = '-';
    
    // If returning back to wild type, delete mutation entry
    if (wildBase === toBase) {
        activeMutations = activeMutations.filter(m => m.pos !== index);
    } else {
        // Remove existing mutation at this position if any
        activeMutations = activeMutations.filter(m => m.pos !== index);
        activeMutations.push({
            pos: index,
            from: wildBase,
            to: toBase
        });
    }
    
    // Re-render
    renderSequenceGrid();
    updatePrediction();
    fetchRecommendations();
}

function resetAllMutations() {
    activeMutations = [];
    document.querySelectorAll('.highlighted-site').forEach(el => el.classList.remove('highlighted-site'));
    renderSequenceGrid();
    updatePrediction();
    fetchRecommendations();
}

// Update prediction results (API call)
async function updatePrediction() {
    if (!currentPromoter) return;
    
    try {
        const payload = {
            wild_seq: currentPromoter.seq,
            mutations: activeMutations.map(m => ({ pos: m.pos, to: m.to })),
            gene_symbol: currentPromoter.symbol,
            systematic_name: currentPromoter.systematic_name
        };
        
        const res = await fetch('/api/predict', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        
        const data = await res.json();
        
        // Update Gauge
        renderGauge(data.predicted_value);
        
        // Update status box
        updateStatusBox(data.predicted_value, data.change_percentage, data.tata_destroyed);
        
        // Update Active mutations list
        renderActiveMutations();
        
        // Update Calibration Spectrum Bar
        renderSpectrum(data.calibrated_wt, data.calibrated_mut, data.references);
        
        // Re-draw track based on mutant sites if needed
        renderPromoterTrack(data.mutant_sites);
        
    } catch (err) {
        console.error("Error predicting expression:", err);
    }
}

// Circular progress gauge
function renderGauge(value) {
    predictionVal.textContent = `${value}%`;
    
    // Circumference = 2 * PI * r = 2 * 3.14159 * 80 = 502.65
    const circ = 502.65;
    
    // Map value (0 to 200%) to dashoffset
    // 0% -> offset = 502.65
    // 100% -> offset = 251.3
    // 200% -> offset = 0
    const maxVal = 200;
    const clampedVal = Math.min(value, maxVal);
    const offset = circ - (circ * clampedVal / maxVal);
    
    gaugeFillArc.style.strokeDashoffset = offset;
    
    // Color transitions based on level
    if (value > 110) {
        gaugeFillArc.style.stroke = 'var(--color-green)';
        gaugeFillArc.style.filter = 'drop-shadow(0 0 8px var(--color-green-glow))';
    } else if (value < 90) {
        gaugeFillArc.style.stroke = 'var(--color-red)';
        gaugeFillArc.style.filter = 'drop-shadow(0 0 8px var(--color-red-glow))';
    } else {
        gaugeFillArc.style.stroke = 'var(--color-cyan)';
        gaugeFillArc.style.filter = 'none';
    }
}

function updateStatusBox(value, change, isTataDestroyed) {
    statusTag.className = 'status-indicator';
    
    if (isTataDestroyed) {
        statusTag.classList.add('indicator-drop');
        statusTag.textContent = '핵심 프로모터 파괴';
        statusDesc.textContent = '핵심 요소인 TATA Box가 손상되었습니다. 전사 개시가 극도로 억제되어 유전자 발현이 거의 일어나지 않을 것입니다.';
    } else if (value > 110) {
        statusTag.classList.add('indicator-boost');
        statusTag.textContent = `발현 강도 향상 (${change > 0 ? '+' : ''}${change}%)`;
        statusDesc.textContent = 'Point mutation에 의해 Activator 모티프가 최적화되었거나 Repressor 억제 부위가 파괴되어 발현 강도가 강화되었습니다.';
    } else if (value < 90) {
        statusTag.classList.add('indicator-drop');
        statusTag.textContent = `발현 강도 하락 (${change}%)`;
        statusDesc.textContent = '핵심 Activator 결합 부위가 손상되어 전사인자 결합률이 떨어지고 프로모터 발현 강도가 하향되었습니다.';
    } else {
        statusTag.classList.add('indicator-normal');
        statusTag.textContent = '활성 변화 없음';
        statusDesc.textContent = '야생형과 생물학적으로 동등한 수준의 프로모터 발현 강도를 나타냅니다.';
    }
}

// Active Mutations list
function renderActiveMutations() {
    mutationsCount.textContent = activeMutations.length;
    
    if (activeMutations.length === 0) {
        activeMutationsList.innerHTML = '<p class="no-mutations-text">적용된 point mutation이 없습니다. 서열 에디터에서 염기를 클릭해 변경해 보세요.</p>';
        return;
    }
    
    activeMutationsList.innerHTML = '';
    activeMutations.forEach(mut => {
        const item = document.createElement('div');
        item.className = 'mut-item';
        
        const coord = 1000 - mut.pos;
        item.innerHTML = `
            <div class="mut-info">
                -${coord} bp: <span>${mut.from}</span> <i class="fa-solid fa-arrow-right" style="font-size:0.7rem; color:var(--text-muted);"></i> <span>${mut.to}</span>
            </div>
            <button class="btn-delete-mut" onclick="deleteMutation(${mut.pos})">
                <i class="fa-solid fa-trash-can"></i>
            </button>
        `;
        activeMutationsList.appendChild(item);
    });
}

function deleteMutation(pos) {
    activeMutations = activeMutations.filter(m => m.pos !== pos);
    renderSequenceGrid();
    updatePrediction();
    fetchRecommendations();
}

// Setup Recommendations Tabs and Load Recommendations
function setupRecommendationTabs() {
    tabUpBtn.addEventListener('click', () => {
        tabUpBtn.classList.add('active');
        tabDownBtn.classList.remove('active');
        activeRecTab = 'up';
        renderRecommendations();
    });
    
    tabDownBtn.addEventListener('click', () => {
        tabDownBtn.classList.add('active');
        tabUpBtn.classList.remove('active');
        activeRecTab = 'down';
        renderRecommendations();
    });
}

async function fetchRecommendations() {
    if (!currentPromoter) return;
    
    try {
        // Calculate the current active sequence
        let currentSeq = currentPromoter.seq;
        currentSeq = analyzer_apply_mutations_js(currentSeq, activeMutations);
        
        const res = await fetch('/api/recommend', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ seq: currentSeq })
        });
        
        recommendedMutations = await res.json();
        renderRecommendations();
        
    } catch (err) {
        console.error("Error fetching recommendations:", err);
    }
}

function renderRecommendations() {
    const list = recommendedMutations[activeRecTab] || [];
    recommendationsList.innerHTML = '';
    
    if (list.length === 0) {
        recommendationsList.innerHTML = '<p class="no-recs-text">더 이상 추천할 추천 변이가 없습니다.</p>';
        return;
    }
    
    list.forEach(rec => {
        const card = document.createElement('div');
        card.className = 'rec-card';
        
        let impactClass = "increase";
        let sign = "+";
        if (rec.impact_type === "high_increase") { impactClass = "high_increase"; sign = "++"; }
        else if (rec.impact_type === "high_decrease") { impactClass = "high_decrease"; sign = "--"; }
        else if (rec.impact_type === "decrease") { impactClass = "decrease"; sign = "-"; }
        
        // Multi-base mutation coordinates differ
        let coordText = `-${1000 - rec.pos} bp`;
        let mutText = `${rec.from} → ${rec.to}`;
        if (rec.is_multi_base) {
            coordText = `-${1000 - rec.pos} ~ -${1000 - (rec.pos + rec.from.length)} bp`;
        }
        
        card.innerHTML = `
            <div class="rec-main">
                <div class="rec-title-row">
                    <span class="rec-effect ${impactClass}">[${sign}] ${rec.effect}</span>
                    <span class="rec-mutation">${coordText}: ${mutText}</span>
                </div>
                <div class="rec-desc">${rec.desc}</div>
            </div>
            <button class="btn-apply-rec" onclick="applyRecommendation(${rec.pos}, '${rec.from}', '${rec.to}', ${rec.is_multi_base || false})">적용</button>
        `;
        recommendationsList.appendChild(card);
    });
}

function applyRecommendation(pos, fromBase, toBase, isMulti) {
    if (isMulti) {
        // Multi base substitution
        for (let i = 0; i < toBase.length; i++) {
            const singlePos = pos + i;
            const singleFrom = fromBase[i];
            const singleTo = toBase[i];
            activeMutations = activeMutations.filter(m => m.pos !== singlePos);
            if (singleFrom !== singleTo) {
                activeMutations.push({
                    pos: singlePos,
                    from: singleFrom,
                    to: singleTo
                });
            }
        }
    } else {
        // Single base mutation
        activeMutations = activeMutations.filter(m => m.pos !== pos);
        activeMutations.push({
            pos: pos,
            from: fromBase,
            to: toBase
        });
    }
    
    renderSequenceGrid();
    updatePrediction();
    fetchRecommendations();
}

// Client-side simulation helpers to save extra requests during recommendation computing
function analyzer_apply_mutations_js(baseSeq, mutations) {
    let arr = baseSeq.split('');
    mutations.forEach(m => {
        if (m.pos >= 0 && m.pos < arr.length) {
            arr[m.pos] = m.to;
        }
    });
    return arr.join('');
}

// Zoom & Drag Pan setup for Promoter Track SVG
function setupTrackZoomPan() {
    const trackWrapper = document.querySelector('.map-track-wrapper');
    if (!trackWrapper) return;

    // 1. Mouse wheel zoom (focus on mouse cursor position in sequence coordinate)
    trackWrapper.addEventListener('wheel', (e) => {
        if (!currentPromoter) return;
        e.preventDefault();
        
        const seqLen = currentPromoter.seq.length;
        const rect = trackSvg.getBoundingClientRect();
        const mouseXRel = (e.clientX - rect.left) / rect.width;
        
        const oldViewWidth = seqLen / zoomLevel;
        const mousePosInSeq = panOffset + mouseXRel * oldViewWidth;
        
        const zoomFactor = 1.2;
        if (e.deltaY < 0) {
            // Zoom in
            zoomLevel = Math.min(15.0, zoomLevel * zoomFactor);
        } else {
            // Zoom out
            zoomLevel = Math.max(1.0, zoomLevel / zoomFactor);
        }
        
        const newViewWidth = seqLen / zoomLevel;
        panOffset = mousePosInSeq - mouseXRel * newViewWidth;
        
        trackWrapper.style.cursor = zoomLevel > 1.0 ? 'grab' : 'default';
        updateTrackViewBox();
    }, { passive: false });
    
    // 2. Mouse drag panning
    let isPanning = false;
    let startX = 0;
    let startPanOffset = 0;
    
    trackWrapper.addEventListener('mousedown', (e) => {
        if (zoomLevel <= 1.0) return; // Only pan when zoomed in
        isPanning = true;
        trackWrapper.style.cursor = 'grabbing';
        startX = e.clientX;
        startPanOffset = panOffset;
        e.preventDefault();
    });
    
    document.addEventListener('mousemove', (e) => {
        if (!isPanning || !currentPromoter) return;
        const seqLen = currentPromoter.seq.length;
        const rect = trackSvg.getBoundingClientRect();
        
        const dxPixels = e.clientX - startX;
        const viewWidth = seqLen / zoomLevel;
        const dxSeq = -(dxPixels / rect.width) * viewWidth;
        
        panOffset = startPanOffset + dxSeq;
        updateTrackViewBox();
    });
    
    document.addEventListener('mouseup', () => {
        if (isPanning) {
            isPanning = false;
            trackWrapper.style.cursor = zoomLevel > 1.0 ? 'grab' : 'default';
        }
    });
}

// Render calibrated standard promoter spectrum bar
function renderSpectrum(wtVal, mutVal, references) {
    const pointer = document.getElementById('spectrum-pointer');
    const tooltip = document.getElementById('pointer-tooltip');
    const wtLabel = document.getElementById('wt-cal-val');
    const mutLabel = document.getElementById('mut-cal-val');
    
    if (!pointer || !tooltip || !wtLabel || !mutLabel) return;
    
    // Calculate Left % coordinates using piecewise linear interpolation
    const leftPercent = getLeftPercentForScore(mutVal);
    
    // Position Pointer
    pointer.style.left = `${leftPercent}%`;
    
    // Update tooltip text
    tooltip.textContent = `현재: ${mutVal}%`;
    
    // Update footer labels
    wtLabel.textContent = `${wtVal}%`;
    mutLabel.textContent = `${mutVal}%`;
    
    // Dynamic pointer color based on mutation vs wt comparison
    const pin = pointer.querySelector('.pointer-pin');
    if (pin) {
        if (mutVal > wtVal + 1) {
            pin.style.background = '#10b981'; // Green for boost
            pin.style.boxShadow = '0 0 10px rgba(16, 185, 129, 0.8)';
        } else if (mutVal < wtVal - 1) {
            pin.style.background = '#ef4444'; // Red for drop
            pin.style.boxShadow = '0 0 10px rgba(239, 68, 68, 0.8)';
        } else {
            pin.style.background = '#7c3aed'; // Purple for WT equivalent
            pin.style.boxShadow = '0 0 8px rgba(124, 58, 237, 0.7)';
        }
    }
}

function getLeftPercentForScore(score) {
    // Piecewise linear interpolation mapping score (0 to 150) to left percent (0% to 100%)
    const points = [
        { score: 0.0, pct: 0.0 },
        { score: 5.0, pct: 5.0 },     // CYC1
        { score: 15.0, pct: 15.0 },   // ACT1
        { score: 40.0, pct: 40.0 },   // ADH1
        { score: 70.0, pct: 70.0 },   // PGK1
        { score: 80.0, pct: 80.0 },   // TEF1
        { score: 100.0, pct: 95.0 },  // TDH3
        { score: 150.0, pct: 100.0 }  // Max
    ];
    
    if (score <= points[0].score) return points[0].pct;
    if (score >= points[points.length - 1].score) return points[points.length - 1].pct;
    
    for (let i = 0; i < points.length - 1; i++) {
        const p1 = points[i];
        const p2 = points[i+1];
        if (score >= p1.score && score <= p2.score) {
            const ratio = (score - p1.score) / (p2.score - p1.score);
            return p1.pct + ratio * (p2.pct - p1.pct);
        }
    }
    return 50.0;
}

