import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
print(" [Figure 6] High-Resolution (600 DPI) Multi-Format Export Script")
print("=" * 70)

# 2. 연도별 진단 건수 및 유전체 유형별 데이터 집계
df['Diagnostic_Phase'] = df['Collection_Year'].apply(lambda x: 'Phase I (Traditional)' if x < 2015 else 'Phase II (HTS Metagenomics)')

if 'Genome_Type' not in df.columns:
    def assign_genome_type(name):
        name_lower = str(name).lower()
        if 'viroid' in name_lower:
            return 'Viroid (ncRNA)'
        elif 'orthotospovirus' in name_lower or 'potyvirus' in name_lower or 'mosaic' in name_lower:
            return 'ssRNA(+) Virus'
        elif 'motto' in name_lower or 'latent' in name_lower:
            return 'ssRNA(-) / dsRNA'
        else:
            return 'DNA Virus'
    df['Genome_Type'] = df['Organism'].apply(assign_genome_type)

yearly_genome = df.groupby(['Collection_Year', 'Genome_Type']).size().unstack(fill_value=0)

# 3. Figure 6: 캔버스 크기 확장 및 글자 크기/두께 강화 (600 DPI 적용)
fig, ax = plt.subplots(figsize=(18, 10), dpi=600)
fig.patch.set_facecolor('#ffffff')

years = yearly_genome.index
columns = yearly_genome.columns
palette = ['#3b82f6', '#10b981', '#f59e0b', '#ef4444']

ax.stackplot(years, yearly_genome.T.values, labels=columns, colors=palette[:len(columns)], alpha=0.85, edgecolor='none')

ax.axvline(x=2015, color='#7c3aed', linestyle='--', linewidth=3.0, label='HTS Mass Surveillance Introduction (2015)')

# 상단 타이틀은 완전히 제거하고 축 레이블과 눈금 폰트 크기/두께 대폭 강화
ax.set_xlabel('Collection Year', fontsize=14, fontweight='bold', labelpad=10)
ax.set_ylabel('Detection Count (Frequency)', fontsize=14, fontweight='bold', labelpad=10)

ax.tick_params(axis='both', labelsize=12, width=2)
for label in ax.get_yticklabels():
    label.set_fontweight('bold')
for label in ax.get_xticklabels():
    label.set_fontweight('bold')

ax.legend(loc='upper left', fontsize=12, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', borderpad=1.2)
ax.grid(True, linestyle=':', alpha=0.6, linewidth=1.2)

plt.tight_layout()

# 4. 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장
base_filename = "Figure_6_Stacked_Area_Chart"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)