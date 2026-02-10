#!/usr/bin/env bash
#SBATCH --job-name=MzTab_Lollipop_V7_Grid
#SBATCH --partition=noninterruptive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=%x_%j.log

# --- ENVIRONMENT SETUP ---
source /data/nasif12/modules_if12/SL7/i12g/miniforge/24.9.0-0/etc/profile.d/conda.sh
conda activate denovo_align_blastpDiamond

# Start Python Heredoc
python3 - << 'EOF_PY'
import pandas as pd
import numpy as np
import os
import glob
import re
import json
import sys
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================

TARGET_MAP = {
    'P00533': 'EGFR',
    'P06737': 'PYGL',
    'P08670': 'VIM'
}

REF_FASTA_PATH = "/s/raw/kusterlab/internal_projects/active/CLINSPECT_WP4/topas_pipeline/MQ_searches/GLIOMA/uniprot_proteome_up000005640_03112020.fasta"

ALIGNMENT_DIRS = [
    "/s/project/denovo-prosit/MoritzRothhaar/cancerAppliKuster/alignment/fullAlign/alignFastas12.12/alignflow06.01/step4enrich/run2fixedIdenti/4runsFilteredenriched_v2/",
    "/data/nasif12/home_if12/rotm/Moe/cancerAppliKuster/runCkpt1Beam/run11.01/align/enrich/4runsFilteredenriched_v2/"
]

PSM_DIR = "/s/project/denovo-prosit/data/TMT-cancer-kusterlab20251119/dnps_glioma_predictions/ckpt1PP_fullrun/ckpt1mzTabMQmerged/extended_tsvs/beam5_rerun/"

OUTPUT_HTML = "Lollipop_Coverage_Report_V7_Grid.html"

ALIGN_FILES = [
    "mq_unique_qcov80max1top100_filtered_enriched.tsv",
    "dnps_unique_COMPLETE_enriched.tsv",
    "shared_peptides_qcov80max1top100_filtered_enriched.tsv",
    "dnps_beam5_new_unique_qcov80max1top100_filtered_enriched.tsv"
]

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def parse_fasta(path, target_ids):
    sequences = {}
    current_id = None
    keep = False
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                match = re.search(r"\|(.*?)\|", line)
                if match:
                    uid = match.group(1)
                    if uid in target_ids:
                        current_id = uid
                        sequences[current_id] = []
                        keep = True
                    else:
                        keep = False
            elif keep:
                sequences[current_id].append(line)
    return {k: "".join(v) for k, v in sequences.items()}

def get_mod_status_detailed(row):
    dnps_pred = str(row.get('DNPS_pred', ''))
    mq_pred = str(row.get('MQ_pred', ''))
    
    dnps_p = '[+79.966]' in dnps_pred
    mq_p = '(Phospho (STY))' in mq_pred
    
    return {
        'dnps_p': dnps_p,
        'mq_p': mq_p,
        'any_p': dnps_p or mq_p,
        'dnps_pred_str': dnps_pred,
        'mq_pred_str': mq_pred
    }

def extract_psite_info(peptide_seq, protein_seq, start_pos, psm_list):
    """
    Robust P-site extraction. 
    1. Finds the mod tag location in the string.
    2. Counts actual amino acids (A-Z) before that tag to get the index.
    3. Maps to protein sequence using start_pos.
    """
    p_sites_found = set()
    
    for psm in psm_list:
        mod_str = ""
        if psm.get('dnps_pred') and '[+79.966]' in psm['dnps_pred']:
            mod_str = psm['dnps_pred']
        elif psm.get('mq_pred') and '(Phospho' in psm['mq_pred']:
            mod_str = psm['mq_pred']
            
        if not mod_str: continue

        # Identify all occurrences of Phospho mods
        # Regex to find mod blocks. 
        # DNPS: [+79.966]
        # MQ: (Phospho (STY))
        
        # We iterate over matches to handle multiple P-sites in one peptide
        # Find all mod patterns
        iter_matches = re.finditer(r'(\[\+.*?\]|\(.*?\))', mod_str)
        
        for match in iter_matches:
            tag = match.group(1)
            # Check if this specific tag is a phospho tag
            if '79.9' in tag or 'Phospho' in tag:
                # The index of the modification start in the string
                mod_start_idx = match.start()
                
                # Get the substring BEFORE this modification
                substring_before = mod_str[:mod_start_idx]
                
                # Count how many Amino Acids (A-Z) are in this substring
                # This ignores characters like '-', '[', ']', numbers, etc.
                aa_count = len(re.findall(r'[A-Z]', substring_before))
                
                # The modified residue is the one immediately preceding the tag.
                # If aa_count is 5, it means there are 5 AAs before the mod.
                # The 5th AA (index 4) is the one modified.
                res_index = aa_count - 1
                
                if res_index >= 0:
                    p_sites_found.add(res_index)

    labels = []
    for pep_idx in sorted(list(p_sites_found)):
        if pep_idx < 0: continue # Should be handled by logic above, but safety first
        
        # Calculate Absolute Position
        # start_pos is 1-based start of the peptide in protein
        abs_pos = start_pos + pep_idx 
        
        # Verify against protein sequence (0-based index is abs_pos - 1)
        if abs_pos <= len(protein_seq):
            aa = protein_seq[abs_pos - 1] 
            labels.append(f"p.{aa}{abs_pos}")
            
    return list(set(labels))

# ==========================================
# MAIN EXECUTION
# ==========================================

print(">>> 1. Setup Targets...")
target_uniprots = set(TARGET_MAP.keys())
print(f"Targets: {TARGET_MAP}")

print(">>> 2. Loading Reference Proteome...")
ref_seqs = parse_fasta(REF_FASTA_PATH, target_uniprots)

print(">>> 3. Processing Alignment Files...")
relevant_sequences = set()
protein_alignments = defaultdict(list)

found_files = [os.path.join(d, f) for d in ALIGNMENT_DIRS for f in ALIGN_FILES if os.path.exists(os.path.join(d, f))]

for fpath in found_files:
    print(f"   Reading {os.path.basename(fpath)}...")
    df = pd.read_csv(fpath, sep='\t', low_memory=False)
    df['clean_sseqid'] = df['sseqid'].str.extract(r"\|(.*?)\|")[0]
    df = df[df['clean_sseqid'].isin(target_uniprots)]
    
    cond0 = df['off_by'] == 0
    cond1 = (df['off_by'] == 1) & (df['snp_explains'] == 1)
    cond2 = (df['off_by'] == 2) & (df['snp_explains'] == 2)
    df = df[cond0 | cond1 | cond2]
    
    for _, row in df.iterrows():
        uid = row['clean_sseqid']
        seq = row['whole_query']
        relevant_sequences.add(seq)
        
        protein_alignments[uid].append({
            'seq': seq, 
            'start': int(row['sstart']), 
            'end': int(row['send'])
        })

print(f"   Total unique relevant peptides: {len(relevant_sequences)}")

print(">>> 4. Streaming PSM Files...")
peptide_stats = defaultdict(lambda: {
    'counts': {'mq': 0, 'dnps': 0, 'both': 0},
    'psite_counts': {'mq': 0, 'dnps': 0},
    'psms': []
})

psm_files = glob.glob(os.path.join(PSM_DIR, "*.tsv"))
total_psms_loaded = 0

for i, psm_file in enumerate(psm_files):
    batch_name = os.path.basename(psm_file)
    chunk_iter = pd.read_csv(psm_file, sep='\t', chunksize=100000, low_memory=False)
    
    for chunk in chunk_iter:
        cols_to_str = ['final_dnps_unmod', 'final_mq_unmod', 'DNPS_pred', 'MQ_pred']
        for c in cols_to_str:
            if c in chunk.columns: chunk[c] = chunk[c].fillna('').astype(str)
        
        mask = chunk['final_dnps_unmod'].isin(relevant_sequences) | chunk['final_mq_unmod'].isin(relevant_sequences)
        filtered = chunk[mask]
        
        if filtered.empty: continue
        total_psms_loaded += len(filtered)
        
        for _, row in filtered.iterrows():
            dnps_seq = row.get('final_dnps_unmod', '')
            mq_seq = row.get('final_mq_unmod', '')
            
            is_dnps_relevant = dnps_seq in relevant_sequences
            is_mq_relevant = mq_seq in relevant_sequences
            
            target_seq, row_type = "", ""
            if is_dnps_relevant and dnps_seq == mq_seq:
                target_seq, row_type = dnps_seq, "both"
            elif is_dnps_relevant and not is_mq_relevant:
                target_seq, row_type = dnps_seq, "dnps"
            elif is_mq_relevant and not is_dnps_relevant:
                target_seq, row_type = mq_seq, "mq"
            elif is_dnps_relevant and is_mq_relevant:
                target_seq, row_type = dnps_seq, "dnps"
            else:
                continue

            stats = peptide_stats[target_seq]
            stats['counts'][row_type] += 1
            
            mod_details = get_mod_status_detailed(row)
            if mod_details['dnps_p']: stats['psite_counts']['dnps'] += 1
            if mod_details['mq_p']: stats['psite_counts']['mq'] += 1
            
            stats['psms'].append({
                'type': row_type, 
                'has_p_row': mod_details['any_p'],
                'scan': row.get('scan_id'), 
                'raw': row.get('rawfile'),
                'dnps_pred': mod_details['dnps_pred_str'],
                'mq_pred': mod_details['mq_pred_str'],
                'dnps_score': row.get('DNPS_score'),
                'mq_prot': row.get('MQ_proteins')
            })
    print(f"   Processed {i+1}/{len(psm_files)} (Total matched PSMs: {total_psms_loaded})", end='\r')

print("\n>>> 5. Building Data Structure for Lollipop Report...")
final_data = {}
for uid, gene in TARGET_MAP.items():
    if uid not in ref_seqs: continue
    
    processed_peptides = []
    max_count_in_protein = 0
    
    for al in protein_alignments.get(uid, []):
        seq = al['seq']
        st = peptide_stats.get(seq)
        
        if not st: continue
        
        counts = st['counts']
        total_count = counts['mq'] + counts['dnps'] + counts['both']
        if total_count == 0: continue
        
        if total_count > max_count_in_protein: max_count_in_protein = total_count
        
        psms = st['psms']
        p_psms = [p for p in psms if p['has_p_row']]
        
        p_labels = []
        if p_psms:
            # Use new robust extraction
            p_labels = extract_psite_info(seq, ref_seqs[uid], al['start'], p_psms)

        processed_peptides.append({
            'seq': seq,
            'start': al['start'],
            'end': al['end'],
            'center': (al['start'] + al['end']) / 2,
            'counts': counts,
            'total_count': total_count,
            'p_labels': p_labels,
            'psms': psms
        })
        
    final_data[uid] = {
        'gene': gene, 
        'sequence': ref_seqs[uid], 
        'peptides': processed_peptides,
        'max_y': max_count_in_protein
    }

print(">>> 6. Generating HTML with Grid & Labels...")

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Lollipop Proteomics</title>
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; background: #f9f9f9; }
        header { background: #2c3e50; color: #fff; padding: 15px; display: flex; justify-content: space-between; align-items: center; }
        h1 { margin: 0; font-size: 22px; }
        .tab-container { background: #fff; padding: 10px 20px; border-bottom: 1px solid #ddd; display: flex; gap: 10px; }
        .tab-btn { padding: 10px 25px; cursor: pointer; border: none; background: #e0e0e0; border-radius: 4px; font-weight: bold; color: #555; transition: 0.3s; }
        .tab-btn.active { background: #3498db; color: white; }
        .tab-btn:hover:not(.active) { background: #d0d0d0; }
        
        .protein-view { display: none; padding: 20px; }
        .protein-view.active { display: block; }
        
        /* LAYOUT FOR CHART */
        .chart-container {
            display: flex;
            background: white;
            border: 1px solid #ccc;
            box-shadow: 0 2px 5px rgba(0,0,0,0.05);
            height: 600px;
            position: relative;
            margin-top: 20px;
        }

        /* Y-AXIS (FIXED) */
        .y-axis-area {
            width: 60px;
            background: #fff;
            border-right: 1px solid #eee;
            position: relative;
            z-index: 100;
            flex-shrink: 0;
        }
        .y-tick-label {
            position: absolute;
            right: 10px;
            font-size: 11px;
            color: #666;
            transform: translateY(50%); /* Center vertically on the line */
        }
        .y-axis-title {
            position: absolute;
            left: 5px;
            top: 50%;
            transform: rotate(-90deg) translateX(50%);
            transform-origin: left center;
            font-size: 12px;
            font-weight: bold;
            color: #333;
            white-space: nowrap;
        }

        /* SCROLLABLE CANVAS */
        .chart-scroll-wrapper {
            flex-grow: 1;
            overflow-x: auto;
            position: relative;
        }
        
        .canvas-area {
            position: relative;
            height: 100%;
            /* Width set dynamically */
        }

        /* GRID LINES */
        .grid-line {
            position: absolute;
            left: 0;
            width: 100%;
            border-bottom: 1px dashed #eee;
            pointer-events: none;
            z-index: 0;
        }

        /* SEQUENCE STRIP */
        .sequence-strip {
            position: absolute;
            bottom: 40px;
            left: 0;
            font-family: 'Courier New', monospace;
            font-size: 14px;
            white-space: nowrap;
            border-top: 2px solid #333;
            padding-top: 5px;
            z-index: 5;
        }
        .aa-char {
            display: inline-block;
            text-align: center;
            color: #ccc;
        }
        .aa-char:nth-child(10n) { color: #333; font-weight: bold; }

        /* LOLLIPOP ELEMENTS */
        .stick {
            position: absolute;
            bottom: 65px; /* Must align with grid zero */
            width: 2px;
            background-color: #888;
            pointer-events: none;
            z-index: 10;
        }
        
        .pop-head {
            position: absolute;
            width: 20px;
            height: 20px;
            border-radius: 50%;
            cursor: pointer;
            box-shadow: 1px 1px 3px rgba(0,0,0,0.3);
            transition: transform 0.2s;
            z-index: 20;
        }
        .pop-head:hover { transform: scale(1.3); z-index: 30; border: 2px solid #333; }
        
        .p-label {
            position: absolute;
            font-size: 10px;
            background: #ffecb3;
            border: 1px solid #ffc107;
            padding: 1px 3px;
            border-radius: 3px;
            white-space: nowrap;
            transform: translate(15px, -10px);
            pointer-events: none;
            color: #d35400;
            font-weight: bold;
            z-index: 25;
        }

        /* LEGEND */
        .legend {
            position: fixed; top: 80px; right: 40px;
            background: rgba(255,255,255,0.9); border: 1px solid #ccc;
            padding: 10px; border-radius: 5px; z-index: 1000; font-size: 12px;
        }
        .legend-item { display: flex; align-items: center; margin-bottom: 5px; }
        .color-box { width: 12px; height: 12px; margin-right: 8px; border-radius: 2px; }
        
        /* MODAL */
        .modal { display: none; position: fixed; z-index: 2000; left: 0; top: 0; width: 100%; height: 100%; background-color: rgba(0,0,0,0.5); }
        .modal-content { background-color: #fff; margin: 5% auto; padding: 20px; border: 1px solid #888; width: 80%; max-height: 80vh; overflow-y: auto; border-radius: 5px; box-shadow: 0 4px 15px rgba(0,0,0,0.3); }
        .close { color: #aaa; float: right; font-size: 28px; font-weight: bold; cursor: pointer; }
        table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 12px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; position: sticky; top: 0; }
        .row-dnps { background-color: #e8f5e9; } .row-mq { background-color: #fdeaea; } .row-both { background-color: #e3f2fd; }

    </style>
</head>
<body>

<header>
    <h1>Lollipop Proteomics (Grid View)</h1>
    <div id="loading">Loading Data...</div>
</header>

<div class="legend">
    <div class="legend-item"><div class="color-box" style="background: #d9534f;"></div>MQ Only</div>
    <div class="legend-item"><div class="color-box" style="background: #5cb85c;"></div>DNPS Only</div>
    <div class="legend-item"><div class="color-box" style="background: #0275d8;"></div>Both</div>
    <div class="legend-item"><div class="color-box" style="background: #ffecb3; border:1px solid orange;"></div>Phospho Label</div>
</div>

<div class="tab-container" id="tab-container"></div>
<div id="content-area"></div>

<div id="myModal" class="modal">
  <div class="modal-content">
    <span class="close" onclick="closeModal()">&times;</span>
    <h2 id="modal-title">Peptide Details</h2>
    <div id="modal-summary" style="margin-bottom: 15px; font-size: 14px;"></div>
    <table>
        <thead><tr><th>Type</th><th>Has P?</th><th>Scan</th><th>Rawfile</th><th>DNPS Pred</th><th>Score</th><th>MQ Pred</th></tr></thead>
        <tbody id="modal-body"></tbody>
    </table>
  </div>
</div>

<script>
    const data = __JSON_DATA__;
    const AA_WIDTH = 12; // Slightly wider for readability
    const CHART_HEIGHT = 450; // Total container height inside scroll wrapper
    const PADDING_BOTTOM = 65; // Matches .stick bottom css
    const PLOT_AREA_HEIGHT = CHART_HEIGHT - PADDING_BOTTOM - 20; // Available height for sticks

    function init() {
        const tabs = document.getElementById('tab-container');
        const content = document.getElementById('content-area');
        let first = true;
        
        for (const [uid, info] of Object.entries(data)) {
            const btn = document.createElement('button');
            btn.className = 'tab-btn' + (first ? ' active' : '');
            btn.innerText = info.gene + ' (' + uid + ')';
            btn.onclick = () => openTab(uid, btn);
            tabs.appendChild(btn);
            
            const div = document.createElement('div');
            div.id = 'view-' + uid;
            div.className = 'protein-view' + (first ? ' active' : '');
            div.innerHTML = renderFullChart(uid, info);
            content.appendChild(div);
            first = false;
        }
        document.getElementById('loading').style.display = 'none';
    }
    
    function openTab(uid, btn) {
        document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        document.querySelectorAll('.protein-view').forEach(v => v.classList.remove('active'));
        document.getElementById('view-' + uid).classList.add('active');
    }

    function calculateTicks(maxY) {
        // Simple logic to find nice round ticks
        let step = 1;
        if (maxY > 10) step = 5;
        if (maxY > 50) step = 10;
        if (maxY > 100) step = 20;
        if (maxY > 500) step = 50;
        if (maxY > 1000) step = 100;
        
        const ticks = [];
        // Go slightly above maxY
        const top = Math.ceil(maxY / step) * step;
        for (let i = 0; i <= top; i += step) {
            ticks.push(i);
        }
        return ticks;
    }

    function renderFullChart(uid, info) {
        const seq = info.sequence;
        const totalWidth = (seq.length * AA_WIDTH) + 100; 
        
        // 1. Calculate Grid
        const ticks = calculateTicks(Math.max(5, info.max_y));
        const maxTick = ticks[ticks.length - 1];
        const pxPerCount = PLOT_AREA_HEIGHT / maxTick;

        // Generate Y-Axis HTML and Grid Lines HTML
        let yAxisHtml = '<div class="y-axis-title">Number of PSMs</div>';
        let gridHtml = '';
        
        ticks.forEach(val => {
            const h = (val * pxPerCount) + PADDING_BOTTOM;
            // Y-Axis Label
            yAxisHtml += `<div class="y-tick-label" style="bottom:${h}px;">${val}</div>`;
            // Grid Line (Background)
            gridHtml += `<div class="grid-line" style="bottom:${h}px;"></div>`;
        });

        // 2. Generate Sequence Strip
        let seqHtml = '<div class="sequence-strip">';
        for(let i=0; i<seq.length; i++) {
            seqHtml += `<span class="aa-char" style="width:${AA_WIDTH}px;">${seq[i]}</span>`;
        }
        seqHtml += '</div>';

        // 3. Generate Lollipops
        let lollipopsHtml = '';
        info.peptides.forEach((pep, idx) => {
            const leftPos = (pep.center - 0.5) * AA_WIDTH;
            const height = pep.total_count * pxPerCount;
            
            // Gradients
            const c = pep.counts;
            const t = pep.total_count;
            const degMQ = (c.mq / t) * 360;
            const degDNPS = (c.dnps / t) * 360;
            
            const gradient = `conic-gradient(
                #d9534f 0deg ${degMQ}deg, 
                #5cb85c ${degMQ}deg ${degMQ + degDNPS}deg, 
                #0275d8 ${degMQ + degDNPS}deg 360deg
            )`;

            lollipopsHtml += `<div class="stick" style="left:${leftPos}px; height:${height}px;"></div>`;
            
            lollipopsHtml += `<div class="pop-head" 
                style="left:${leftPos - 10}px; bottom:${height + PADDING_BOTTOM - 10}px; background:${gradient};"
                onclick="showPopup('${uid}', ${idx})"
                title="Seq: ${pep.seq}\nCount: ${pep.total_count}"></div>`;
                
            if(pep.p_labels && pep.p_labels.length > 0) {
                const labelText = pep.p_labels.join(', ');
                lollipopsHtml += `<div class="p-label" style="left:${leftPos - 10}px; bottom:${height + PADDING_BOTTOM + 15}px;">${labelText}</div>`;
            }
        });

        return `
            <div class="chart-container">
                <div class="y-axis-area" style="height:${CHART_HEIGHT}px;">
                    ${yAxisHtml}
                </div>
                <div class="chart-scroll-wrapper" style="height:${CHART_HEIGHT}px;">
                    <div class="canvas-area" style="width:${totalWidth}px; height:${CHART_HEIGHT}px;">
                        ${gridHtml}
                        ${seqHtml}
                        ${lollipopsHtml}
                    </div>
                </div>
            </div>
        `;
    }

    function showPopup(uid, pepIdx) {
        const pep = data[uid].peptides[pepIdx];
        document.getElementById('modal-title').innerText = 'Details: ' + pep.seq;
        
        let summary = `<b>Total Cases:</b> ${pep.total_count}<br>`;
        summary += `<b>Breakdown:</b> <span style="color:#d9534f">MQ: ${pep.counts.mq}</span> | <span style="color:#5cb85c">DNPS: ${pep.counts.dnps}</span> | <span style="color:#0275d8">Both: ${pep.counts.both}</span><br>`;
        summary += `<b>Position:</b> ${pep.start} - ${pep.end}<br>`;
        if(pep.p_labels.length > 0) summary += `<b>P-Sites:</b> ${pep.p_labels.join(', ')}`;
        
        document.getElementById('modal-summary').innerHTML = summary;
        
        const tbody = document.getElementById('modal-body');
        tbody.innerHTML = '';
        pep.psms.forEach(row => {
            let tr = document.createElement('tr');
            tr.className = 'row-' + row.type;
            let inner = `<td>${row.type}</td>
                        <td>${row.has_p_row ? '<b>Yes</b>' : 'No'}</td>
                        <td>${row.scan}</td>
                        <td>${row.raw}</td>
                        <td>${row.dnps_pred}</td>
                        <td>${parseFloat(row.dnps_score || 0).toFixed(2)}</td>
                        <td>${row.mq_pred}</td>`;
            tr.innerHTML = inner;
            tbody.appendChild(tr);
        });
        
        document.getElementById('myModal').style.display = "block";
    }
    
    function closeModal() { document.getElementById('myModal').style.display = "none"; }
    window.onclick = (event) => { if (event.target == document.getElementById('myModal')) closeModal(); };
    
    if(typeof data !== 'undefined') init();
</script>
</body>
</html>
"""

# Inject Data
json_data = json.dumps(final_data)
html_output = HTML_TEMPLATE.replace('__JSON_DATA__', json_data)

with open(OUTPUT_HTML, "w", encoding='utf-8') as f:
    f.write(html_output)

print(f"\n>>> Done! Lollipop Report generated at {OUTPUT_HTML}")
EOF_PY
