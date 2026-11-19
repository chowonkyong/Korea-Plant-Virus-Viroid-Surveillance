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
print(" [Figure 10] Top Genera and Species Temporal Stacked Bar Charts (600 DPI)")
print("=" * 70)

# 시대(Era) 분류 함수 정의
def assign_era(yr):
    if yr < 2000:
        return '1990s'
    elif yr < 2010:
        return '2000s'
    elif yr < 2020:
        return '2010s'
    else:
        return '2020s'

df['Era'] = df['Collection_Year'].apply(assign_era)
era_palette = {'1990s': '#2563eb', '2000s': '#d97706', '2010s': '#059669', '2020s': '#dc2626'}

# 2. 패널 A: 상위 10대 속(Genus) 및 시대별 누적 데이터 집계
top_genera = df['Genus'].value_counts().head(10).index
df_genera = df[df['Genus'].isin(top_genera)]
cross_genera = pd.crosstab(df_genera['Genus'], df_genera['Era'])[['1990s', '2000s', '2010s', '2020s']]
cross_genera['Total'] = cross_genera.sum(axis=1)
cross_genera = cross_genera.sort_values(by='Total', ascending=True).drop(columns=['Total'])

# 3. 패널 B: 상위 10대 종(Species/Organism) 및 시대별 누적 데이터 집계
top_species = df['Organism'].value_counts().head(10).index
df_species = df[df['Organism'].isin(top_species)]
cross_species = pd.crosstab(df_species['Organism'], df_species['Era'])[['1990s', '2000s', '2010s', '2020s']]
cross_species['Total'] = cross_species.sum(axis=1)
cross_species = cross_species.sort_values(by='Total', ascending=True).drop(columns=['Total'])

# 4. 2단 패널 캔버스 설정 (600 DPI 적용)
fig, axes = plt.subplots(1, 2, figsize=(22, 11), dpi=600)
fig.patch.set_facecolor('#ffffff')

# 패널 (A): Top 10 Genera Stacked Horizontal Bar Chart
cross_genera.plot(kind='barh', stacked=True, ax=axes[0], color=[era_palette[col] for col in cross_genera.columns], edgecolor='none', alpha=0.85, width=0.7)
axes[0].set_title('A. Temporal Distribution of Top 10 Viral & Viroid Genera', fontsize=13, fontweight='bold', pad=15)
axes[0].set_xlabel('Detection Frequency (Count)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Viral & Viroid Genus', fontsize=12, fontweight='bold')
axes[0].tick_params(axis='both', labelsize=11)
for label in axes[0].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[0].get_xticklabels():
    label.set_fontweight('bold')
axes[0].legend(title='Era', fontsize=10, title_fontsize=11, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', loc='lower right')
axes[0].grid(axis='x', linestyle=':', alpha=0.6, linewidth=1.2)

# 패널 (B): Top 10 Species Stacked Horizontal Bar Chart
cross_species.plot(kind='barh', stacked=True, ax=axes[1], color=[era_palette[col] for col in cross_species.columns], edgecolor='none', alpha=0.85, width=0.7)
axes[1].set_title('B. Temporal Distribution of Top 10 Viral & Viroid Species', fontsize=13, fontweight='bold', pad=15)
axes[1].set_xlabel('Detection Frequency (Count)', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Viral & Viroid Species (Organism)', fontsize=12, fontweight='bold')
axes[1].tick_params(axis='both', labelsize=11)
for label in axes[1].get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontstyle('italic') # 수정된 부분: set_fontweight와 set_fontstyle 분리 적용
for label in axes[1].get_xticklabels():
    label.set_fontweight('bold')
axes[1].legend(title='Era', fontsize=10, title_fontsize=11, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', loc='lower right')
axes[1].grid(axis='x', linestyle=':', alpha=0.6, linewidth=1.2)

plt.tight_layout()

# 5. 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장
base_filename = "Figure_10_Genera_Species_Temporal"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)
print("=" * 70)