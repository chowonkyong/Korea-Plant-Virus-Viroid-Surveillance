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
print(" [Figure 8] High-Resolution (600 DPI) Multi-Format Export Script")
print("=" * 70)

# 캔버스 크기를 대폭 키우고 600 DPI 적용
fig, axes = plt.subplots(2, 2, figsize=(20, 16), dpi=600)
fig.patch.set_facecolor('#ffffff')

# Panel A: 주요 기주 식물별 검출 빈도 상위 5개 바차트
top_hosts = df['Standardized_Host'].value_counts().head(5)
bars_a = axes[0, 0].bar(top_hosts.index, top_hosts.values, color='#3b82f6', edgecolor='#1d4ed8', alpha=0.85, width=0.6)
axes[0, 0].set_title('A. Pathogen Detection Frequency by Major Host Plants', fontsize=13, fontweight='bold', pad=15)
axes[0, 0].set_ylabel('Detection Count', fontsize=12, fontweight='bold')
axes[0, 0].tick_params(axis='x', rotation=15, labelsize=11)
axes[0, 0].tick_params(axis='y', labelsize=11)
for label in axes[0, 0].get_xticklabels():
    label.set_fontweight('bold')
for label in axes[0, 0].get_yticklabels():
    label.set_fontweight('bold')
axes[0, 0].grid(True, linestyle=':', alpha=0.6, linewidth=1.2)

# Panel B: 기주 식물(상위 5개)과 유전체 타입 간의 교차탭 히트맵
top_host_list = top_hosts.index.tolist()
filtered_host_df = df[df['Standardized_Host'].isin(top_host_list)]
host_genome_cross = pd.crosstab(filtered_host_df['Standardized_Host'], filtered_host_df['Genome_Type'])
sns.heatmap(host_genome_cross, annot=True, fmt='d', cmap='YlGnBu', cbar=True, ax=axes[0, 1], linewidths=1.0, linecolor='#cbd5e1', annot_kws={'size': 11, 'weight': 'bold'})
axes[0, 1].set_title('B. Cross-Tabulation Heatmap: Major Hosts vs. Genome Type', fontsize=13, fontweight='bold', pad=15)
axes[0, 1].set_xlabel('Genome Type', fontsize=12, fontweight='bold')
axes[0, 1].set_ylabel('Standardized Host', fontsize=12, fontweight='bold')
axes[0, 1].tick_params(axis='x', rotation=25, labelsize=11)
axes[0, 1].tick_params(axis='y', labelsize=11)
for label in axes[0, 1].get_xticklabels():
    label.set_fontweight('bold')
for label in axes[0, 1].get_yticklabels():
    label.set_fontweight('bold')

# Panel C: 유전체 분절 구조와 벡터 전파 모드 연관성 분석 (가로형 바차트)
seg_vector = pd.crosstab(df['Vector_Mode'], df['Genome_Segmentation'])
seg_vector.plot(kind='barh', ax=axes[1, 0], colormap='viridis', alpha=0.85, edgecolor='none', width=0.75)
axes[1, 0].set_title('C. Association Between Genome Segmentation and Vector Transmission', fontsize=13, fontweight='bold', pad=15)
axes[1, 0].set_xlabel('Frequency', fontsize=12, fontweight='bold')
axes[1, 0].set_ylabel('Vector Transmission Mode', fontsize=12, fontweight='bold')
axes[1, 0].tick_params(axis='y', labelsize=10)
axes[1, 0].tick_params(axis='x', labelsize=11)
for label in axes[1, 0].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[1, 0].get_xticklabels():
    label.set_fontweight('bold')
axes[1, 0].legend(title='Segmentation', fontsize=9, title_fontsize=10, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', loc='lower right')
axes[1, 0].grid(True, linestyle=':', alpha=0.6, linewidth=1.2)

# Panel D: 관측된 주요 병원체 증상 빈도 상위 5개 바차트
top_symptoms = df['Pathogenicity_Symptoms'].value_counts().head(5)
axes[1, 1].barh(top_symptoms.index, top_symptoms.values, color='#10b981', edgecolor='#047857', alpha=0.85, height=0.6)
axes[1, 1].set_title('D. Symptom Phenotype Frequency Distribution in Surveillance', fontsize=13, fontweight='bold', pad=15)
axes[1, 1].set_xlabel('Detection Count', fontsize=12, fontweight='bold')
axes[1, 1].set_ylabel('Observed Symptom Phenotype', fontsize=12, fontweight='bold')
axes[1, 1].tick_params(axis='y', labelsize=11)
axes[1, 1].tick_params(axis='x', labelsize=11)
for label in axes[1, 1].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[1, 1].get_xticklabels():
    label.set_fontweight('bold')
axes[1, 1].grid(True, linestyle=':', alpha=0.6, linewidth=1.2)

# 상단 Figure 8 전체 타이틀은 완전히 삭제하고 패널 간격 최적화
plt.tight_layout()

# 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장
base_filename = "Figure_8_True_MultiPanel_Composite"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)
print("=" * 70)
