import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

csv_path = 'korea_plant_viruses_core_perfect.csv'
if not os.path.exists(csv_path):
    for root, dirs, files in os.walk('.'):
        if 'korea_plant_viruses_core_perfect.csv' in files:
            csv_path = os.path.join(root, 'korea_plant_viruses_core_perfect.csv')
            break

df = pd.read_csv(csv_path)

print("=" * 70)
print(" [Figure 8] Streamlined 3-Panel (A, B, C) High-Resolution Export Script")
print("=" * 70)

# 1행 3열 패널 레이아웃 구성
fig, axes = plt.subplots(1, 3, figsize=(24, 7), dpi=600)
fig.patch.set_facecolor('#ffffff')

# ----------------------------------------------------
# Panel A: Cross-Tabulation Heatmap: Major Hosts vs. Genome Type
# ----------------------------------------------------
top_hosts = df['Standardized_Host'].value_counts().head(5)
top_host_list = top_hosts.index.tolist()
filtered_host_df = df[df['Standardized_Host'].isin(top_host_list)]
host_genome_cross = pd.crosstab(filtered_host_df['Standardized_Host'], filtered_host_df['Genome_Type'])

sns.heatmap(
    host_genome_cross,
    annot=True,
    fmt='d',
    cmap='YlGnBu',
    cbar=True,
    ax=axes[0],
    linewidths=1.0,
    linecolor='#cbd5e1',
    annot_kws={'size': 11, 'weight': 'bold'}
)
axes[0].set_title('Major Hosts vs. Genome Type Heatmap', fontsize=13, fontweight='bold', pad=15)
axes[0].set_xlabel('Genome Type', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Standardized Host', fontsize=12, fontweight='bold')
axes[0].tick_params(axis='x', rotation=25, labelsize=11)
axes[0].tick_params(axis='y', labelsize=11)
for label in axes[0].get_xticklabels():
    label.set_fontweight('bold')
for label in axes[0].get_yticklabels():
    label.set_fontweight('bold')
axes[0].text(-0.22, -0.16, 'A', transform=axes[0].transAxes, fontsize=18, fontweight='bold', va='bottom', ha='left', color='#0f172a')

# ----------------------------------------------------
# Panel B: Association Between Genome Segmentation and Vector Transmission
# ----------------------------------------------------
seg_vector = pd.crosstab(df['Vector_Mode'], df['Genome_Segmentation'])
seg_vector.plot(kind='barh', ax=axes[1], colormap='viridis', alpha=0.85, edgecolor='none', width=0.75)
axes[1].set_title('Genome Segmentation vs. Vector Transmission', fontsize=13, fontweight='bold', pad=15)
axes[1].set_xlabel('Frequency', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Vector Transmission Mode', fontsize=12, fontweight='bold')
axes[1].tick_params(axis='y', labelsize=10)
axes[1].tick_params(axis='x', labelsize=11)
for label in axes[1].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[1].get_xticklabels():
    label.set_fontweight('bold')
axes[1].legend(title='Segmentation', fontsize=8, title_fontsize=9, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', loc='lower right')
axes[1].grid(True, linestyle=':', alpha=0.6, linewidth=1.2)
axes[1].text(-0.28, -0.16, 'B', transform=axes[1].transAxes, fontsize=18, fontweight='bold', va='bottom', ha='left', color='#0f172a')

# ----------------------------------------------------
# Panel C: Refined Primary Symptom Phenotype Distribution (Mutual Exclusive)
# ----------------------------------------------------
def categorize_symptom(val):
    if pd.isna(val):
        return 'Unspecified / General'
    val_lower = str(val).lower()
    if 'yellow' in val_lower or 'chlorosis' in val_lower or 'chlorotic' in val_lower:
        return 'Chlorosis / Yellowing'
    elif 'mosaic' in val_lower or 'mottling' in val_lower:
        return 'Mosaic / Mottling'
    elif 'stunting' in val_lower or 'distortion' in val_lower or 'curling' in val_lower:
        return 'Stunting / Distortion'
    elif 'necrosis' in val_lower or 'ring spot' in val_lower or 'ringspot' in val_lower:
        return 'Necrosis / Ringspot'
    elif 'latent' in val_lower or 'asymptomatic' in val_lower:
        return 'Latent / Asymptomatic'
    else:
        return 'Unspecified / General'

df['Primary_Symptom'] = df['Pathogenicity_Symptoms'].apply(categorize_symptom)
primary_symptoms = df['Primary_Symptom'].value_counts()

axes[2].barh(primary_symptoms.index, primary_symptoms.values, color='#10b981', edgecolor='#047857', alpha=0.85, height=0.6)
axes[2].set_title('Primary Symptom Phenotype Distribution', fontsize=13, fontweight='bold', pad=15)
axes[2].set_xlabel('Detection Count', fontsize=12, fontweight='bold')
axes[2].set_ylabel('Primary Symptom Category', fontsize=12, fontweight='bold')
axes[2].tick_params(axis='y', labelsize=11)
axes[2].tick_params(axis='x', labelsize=11)
for label in axes[2].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[2].get_xticklabels():
    label.set_fontweight('bold')
axes[2].grid(True, linestyle=':', alpha=0.6, linewidth=1.2)
axes[2].text(-0.26, -0.16, 'C', transform=axes[2].transAxes, fontsize=18, fontweight='bold', va='bottom', ha='left', color='#0f172a')

plt.subplots_adjust(wspace=0.58)

# 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장
base_filename = "Figure_8_Streamlined_3Panel"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)
print("=" * 70)