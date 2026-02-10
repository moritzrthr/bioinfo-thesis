#!/usr/bin/env bash
#SBATCH --job-name=PSM_Rank_Monomethyl_Analysis
#SBATCH --partition=noninterruptive
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=%x_%j.log

# Load Conda Environment
source /data/nasif12/modules_if12/SL7/i12g/miniforge/24.9.0-0/etc/profile.d/conda.sh
conda activate denovo_align_blastpDiamond

python3 - << 'EOF_PY'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import re

# --- CONFIGURATION ---
PSM_DIR = "/s/project/denovo-prosit/data/TMT-cancer-kusterlab20251119/dnps_glioma_predictions/ckpt1PP_fullrun/ckpt1mzTabMQmerged/extended_tsvs/beam5_rerun/"
ALIGN_BASE = "/s/project/denovo-prosit/MoritzRothhaar/cancerAppliKuster/alignment/fullAlign/alignFastas12.12/alignflow06.01/step4enrich/run2fixedIdenti/4runsFilteredenriched_v2/"

ALIGN_FILES = [
    os.path.join(ALIGN_BASE, "dnps_unique_COMPLETE_enriched.tsv"),
    os.path.join(ALIGN_BASE, "mq_unique_qcov80max1top100_filtered_enriched.tsv"),
    os.path.join(ALIGN_BASE, "shared_peptides_qcov80max1top100_filtered_enriched.tsv"),
    "/data/nasif12/home_if12/rotm/Moe/cancerAppliKuster/runCkpt1Beam/run11.01/align/enrich/4runsFilteredenriched_v2/dnps_beam5_new_unique_qcov80max1top100_filtered_enriched.tsv"
]

# --- 1. LOAD PSM DATA ---
print("Loading PSM TSV files...")
psm_files = glob.glob(os.path.join(PSM_DIR, "*.tsv"))
df_psms_list = []

for f in psm_files:
    tmp = pd.read_csv(f, sep='\t', usecols=['DNPS_pred', 'DNPS_score', 'final_dnps_unmod'])
    df_psms_list.append(tmp)

df_psms = pd.concat(df_psms_list, ignore_index=True)

# PTM Classification Logik inkl. Monomethylierung
# TMT/Standard Mods: [+15.995]=Ox, [+42.011]=Ac, [+229.163]=TMT6, [+271.174]=TMTpro
unmod_ox_ac = {"[+15.995]", "[+42.011]", "[+229.163]", "[+271.174]"}

def get_ptm_cat(seq_str):
    if pd.isna(seq_str): return "Unknown"
    # Monomethylierung Check (+14.016)
    if "[+14.016]" in seq_str:
        return "Monomethylation"
    # Phosphorylierung Check (+79.966)
    if "[+79.966]" in seq_str:
        return "Phosphorylation"
    
    mods = re.findall(r"\[\+.*?\]", seq_str)
    if not mods or all(m in unmod_ox_ac for m in mods):
        return "Unmodified or Ox/Ac"
    return "Other PTM type"

print("Classifying PTM types...")
df_psms['ptm_cat'] = df_psms['DNPS_pred'].apply(get_ptm_cat)

# --- 2. LOAD ALIGNMENT DATA ---
print("Loading Alignment files...")
align_list = []
for f in ALIGN_FILES:
    if os.path.exists(f):
        tmp_align = pd.read_csv(f, sep='\t', usecols=['whole_query', 'off_by'])
        align_list.append(tmp_align)

df_align = pd.concat(align_list, ignore_index=True)
df_align = df_align.sort_values('off_by').drop_duplicates('whole_query')

# --- 3. MERGE ---
print("Merging PSMs with Alignments...")
df_master = pd.merge(df_psms, df_align, left_on='final_dnps_unmod', right_on='whole_query', how='left')

# --- 4. PLOTTING (Off-by-0 Only) ---
def create_rank_plot_refined(df_input):
    df = df_input.sort_values('DNPS_score', ascending=False).copy()
    df['rank'] = np.arange(1, len(df) + 1)
    
    # Stratifizierung
    categories = ["Monomethylation", "Phosphorylation", "Unmodified or Ox/Ac", "Other PTM type"]
    colors = {
        "Monomethylation": "#1f77b4", # Blau
        "Phosphorylation": "#d62728", # Rot
        "Unmodified or Ox/Ac": "#2ca02c", # Grün
        "Other PTM type": "#9467bd", # Lila
    }
    
    thresholds = [0.80, 0.90, 0.95]
    markers = {0.80: 'D', 0.90: '*', 0.95: 'h'}
    
    fig, ax = plt.subplots(figsize=(10, 7))

    for cat in categories:
        temp_df = df[df['ptm_cat'] == cat].copy()
        
        if temp_df.empty: continue

        # Off-by-0 Logic
        temp_df['is_correct'] = (temp_df['off_by'] == 0).astype(int)
        temp_df['running_prop'] = temp_df['is_correct'].expanding().mean()
        
        # Plot Main Line
        line, = ax.plot(temp_df['rank'], temp_df['running_prop'], 
                        label=cat, color=colors[cat], linewidth=2)
        
        # Marker for Scores
        for t in thresholds:
            t_subset = temp_df[temp_df['DNPS_score'] >= t]
            if not t_subset.empty:
                last_point = t_subset.iloc[-1]
                ax.scatter(last_point['rank'], last_point['running_prop'], 
                           marker=markers[t], color=colors[cat], s=100, zorder=5)

    # Formatting
    ax.set_xscale('log')
    ax.set_xlabel("PSM Rank", fontsize=12)
    ax.set_ylabel("Proportion of PSMs with perfect blastp alignment", fontsize=12)
    ax.grid(True, which="both", ls="-", alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_title("De novo Performance: Perfect blastp Alignment (off_by=0)", fontweight='bold')

    # Custom Legend Handling
    handles, labels = ax.get_legend_handles_labels()
    
    # Add Score Threshold Icons to Legend
    from matplotlib.lines import Line2D
    legend_elements = handles # Start with category lines
    legend_elements.append(Line2D([0], [0], color='w', label='')) # Spacer
    legend_elements.append(Line2D([0], [0], marker=markers[0.80], color='gray', label='Score ≥ 0.80', markersize=10, linestyle='None'))
    legend_elements.append(Line2D([0], [0], marker=markers[0.90], color='gray', label='Score ≥ 0.90', markersize=10, linestyle='None'))
    legend_elements.append(Line2D([0], [0], marker=markers[0.95], color='gray', label='Score ≥ 0.95', markersize=10, linestyle='None'))

    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
    
    plt.tight_layout()
    plt.savefig("rank_precision_offby0_stratified.png", dpi=300, bbox_inches='tight')
    plt.close()

print("Generating Refined Plot...")
create_rank_plot_refined(df_master)
print("Analysis Complete.")
EOF_PY
