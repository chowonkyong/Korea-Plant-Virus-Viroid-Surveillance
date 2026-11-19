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
print(" [Figure 7] High-Resolution (600 DPI) Multi-Format Export Script")
print("=" * 70)

# 유전체 유형 컬럼 정의
if 'Genome_Type' not in df.columns:
    def assign_genome_type(name):
        name_lower = str(name).lower()
        if 'viroid' in name_lower:
            return 'viroids'
        elif 'orthotospovirus' in name_lower or 'potyvirus' in name_lower or 'mosaic' in name_lower:
            return 'ssRNA (positive)'
        elif 'motto' in name_lower or 'latent' in name_lower:
            return 'ssRNA (negative)'
        else:
            return 'dsDNA'
    df['Genome_Type'] = df['Organism'].apply(assign_genome_type)

pivot_df = pd.crosstab(df['Genome_Type'], df['Collection_Year'])
all_years = range(int(df['Collection_Year'].min()), int(df['Collection_Year'].max()) + 1)
pivot_df = pivot_df.reindex(columns=all_years, fill_value=0)

# 2. Figure 7 히트맵 시각화 렌더링 (600 DPI 고해상도 및 확장된 캔버스 적용)
fig, ax = plt.subplots(figsize=(18, 9), dpi=600)
fig.patch.set_facecolor('#ffffff')

sns.heatmap(pivot_df, cmap='YlGnBu', annot=True, fmt='d', linewidths=1.0, 
            linecolor='#cbd5e1', cbar_kws={'label': 'Detection Frequency'}, ax=ax,
            annot_kws={'size': 11, 'weight': 'bold'})

# 상단 타이틀은 완전히 삭제하고 축 레이블과 눈금 폰트 크기/두께 대폭 강화
ax.set_xlabel('Collection Year', fontsize=14, fontweight='bold', labelpad=12)
ax.set_ylabel('Genome Type', fontsize=14, fontweight='bold', labelpad=12)

plt.xticks(rotation=45, ha='right', fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')

# 컬러바 레이블 폰트 두께 조정
cbar = ax.collections[0].colorbar
cbar.ax.yaxis.label.set_fontweight('bold')
cbar.ax.yaxis.label.set_fontsize(12)
for label in cbar.ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()

# 3. 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장
base_filename = "Figure_7_Temporal_Heatmap"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)