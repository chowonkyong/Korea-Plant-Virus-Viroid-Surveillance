import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 1. 데이터 로드 및 전처리
csv_path = 'korea_plant_viruses_core_perfect.csv'
if not os.path.exists(csv_path):
    for root, dirs, files in os.walk('.'):
        if 'korea_plant_viruses_core_perfect.csv' in files:
            csv_path = os.path.join(root, 'korea_plant_viruses_core_perfect.csv')
            break

df = pd.read_csv(csv_path)

print("=" * 70)
print(" [Graphical Abstract - Elsevier Compliant] 531 x 1328 Proportion (13.28 x 5.31 in, 600 DPI)")
print("=" * 70)

# 스타일 및 글씨 두께 설정
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Elsevier 규격에 맞춘 가로세로 비율 (Width: 13.28 inches, Height: 5.31 inches -> 2.5:1 비율)
fig, axes = plt.subplots(1, 3, figsize=(13.28, 5.31), dpi=600)
fig.patch.set_facecolor('#ffffff')

# --- Panel 1: National Virome Dataset (Genome Types) ---
g_counts = df['Genome_Type'].value_counts().head(5)
sns.barplot(x=g_counts.values, y=g_counts.index, ax=axes[0], hue=g_counts.index, palette='Blues_r', legend=False)
axes[0].set_title('1. National Virome Dataset\n(6,538 Records, 1990–2026)', fontsize=10, fontweight='bold', pad=10)
axes[0].set_xlabel('Number of Records', fontsize=9, fontweight='bold')
axes[0].set_ylabel('Genome Type', fontsize=9, fontweight='bold')
axes[0].tick_params(axis='both', labelsize=8)
for tick in axes[0].get_yticklabels():
    tick.set_fontweight('bold')
for tick in axes[0].get_xticklabels():
    tick.set_fontweight('bold')

total_n = len(df)
for i, v in enumerate(g_counts.values):
    pct = (v / total_n) * 100
    axes[0].text(v + 50, i, f'{pct:.1f}%', va='center', fontsize=8, fontweight='bold')

# --- Panel 2: Diagnostic Paradigm Shift (Pre-2015 vs HTS Era) ---
yearly_df = df.groupby(['Collection_Year', df['Collection_Year'].apply(lambda x: 'Phase I (<2015)' if x < 2015 else 'Phase II (HTS ≥2015)')]).size().unstack(fill_value=0)
yearly_df.plot(kind='bar', stacked=True, ax=axes[1], color=['#94a3b8', '#3b82f6'], width=0.8, edgecolor='none')
axes[1].set_title('2. Diagnostic Paradigm Shift\n(Pre-2015 vs. HTS Era)', fontsize=10, fontweight='bold', pad=10)
axes[1].set_xlabel('Collection Year', fontsize=9, fontweight='bold')
axes[1].set_ylabel('Detection Frequency', fontsize=9, fontweight='bold')
axes[1].tick_params(axis='both', labelsize=8)
plt.sca(axes[1])
plt.xticks(ticks=range(0, len(yearly_df.index), 5), labels=yearly_df.index[::5], rotation=45, fontweight='bold')
for tick in axes[1].get_yticklabels():
    tick.set_fontweight('bold')
axes[1].legend(fontsize=8, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', loc='upper left')
axes[1].grid(axis='y', linestyle=':', alpha=0.6)

# --- Panel 3: Primary Horticultural Hotspots ---
top_hosts = df[df['Standardized_Host'] != 'Unknown']['Standardized_Host'].value_counts().head(5)
sns.barplot(x=top_hosts.values, y=top_hosts.index, ax=axes[2], hue=top_hosts.index, palette='Purples_r', legend=False)
axes[2].set_title('3. Primary Horticultural Hotspots\n& Vector Pathways', fontsize=10, fontweight='bold', pad=10)
axes[2].set_xlabel('Record Count', fontsize=9, fontweight='bold')
axes[2].set_ylabel('Host Plant Species', fontsize=9, fontweight='bold')
axes[2].tick_params(axis='both', labelsize=8)
for tick in axes[2].get_yticklabels():
    tick.set_fontweight('bold')
    tick.set_style('italic')
for tick in axes[2].get_xticklabels():
    tick.set_fontweight('bold')

for i, v in enumerate(top_hosts.values):
    axes[2].text(v + 15, i, str(v), va='center', fontsize=8, fontweight='bold')

plt.tight_layout()

# 600 DPI 고해상도로 TIFF, PDF, SVG 포맷 동시 저장 (Elsevier 권장 포맷 반영)
base_filename = "Graphical_Abstract_Elsevier"
for ext in ['tiff', 'pdf', 'svg', 'png']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI, Widescreen 2.5:1 Ratio)")

plt.close(fig)
print("=" * 70)