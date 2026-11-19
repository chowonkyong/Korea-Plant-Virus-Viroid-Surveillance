import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 스타일 및 글씨 두께 설정
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# 실제 마스터 데이터셋 로드
df = pd.read_csv('korea_plant_viruses_core_perfect.csv')

# Chrysanthemum sp.를 Chrysanthemum morifolium으로 변경
df['Standardized_Host'] = df['Standardized_Host'].replace({'Chrysanthemum sp.': 'Chrysanthemum morifolium'})

total_n = len(df)

# (A) 유전체 유형별 분포 비율 바차트
g_counts = df['Genome_Type'].value_counts()
g_labels = g_counts.index.tolist()
g_vals = g_counts.values.tolist()

sns.barplot(x=g_vals, y=g_labels, ax=axes[0], hue=g_labels, palette='Blues_r', legend=False)
axes[0].set_title('A. Distribution of Viral Genome Types', fontsize=15, fontweight='bold')
axes[0].set_xlabel('Number of Records', fontsize=13, fontweight='bold')
axes[0].set_ylabel('Genome Type', fontsize=13, fontweight='bold')
axes[0].tick_params(axis='both', which='major', labelsize=11)
for tick in axes[0].get_yticklabels():
    tick.set_fontweight('bold')
for tick in axes[0].get_xticklabels():
    tick.set_fontweight('bold')

for i, v in enumerate(g_vals):
    pct = (v / total_n) * 100
    axes[0].text(v + 50, i, f'{v} ({pct:.1f}%)', va='center', fontsize=11, fontweight='bold')

# (B) 주요 기주 식물(Host) 순위 바차트
host_counts = df[df['Standardized_Host'] != 'Unknown']['Standardized_Host'].value_counts().head(6)
sns.barplot(x=host_counts.values, y=host_counts.index, ax=axes[1], hue=host_counts.index, palette='Purples_r', legend=False)
axes[1].set_title('B. Top 6 Host Plants', fontsize=15, fontweight='bold')
axes[1].set_xlabel('Number of Records', fontsize=13, fontweight='bold')
axes[1].set_ylabel('Host Plant Species', fontsize=13, fontweight='bold')
axes[1].tick_params(axis='both', which='major', labelsize=11)
for tick in axes[1].get_yticklabels():
    tick.set_fontweight('bold')
    tick.set_style('italic')
for tick in axes[1].get_xticklabels():
    tick.set_fontweight('bold')

for i, v in enumerate(host_counts.values):
    axes[1].text(v + 10, i, str(v), va='center', fontsize=11, fontweight='bold')

plt.tight_layout()

# 고해상도 600 DPI PNG 및 논문 출판용 벡터 포맷(PDF, SVG)으로 저장
plt.savefig('Figure_1_AB_Overview_Metadata_Distribution.png', dpi=600, bbox_inches='tight')
plt.savefig('Figure_1_AB_Overview_Metadata_Distribution.pdf', format='pdf', bbox_inches='tight')
plt.savefig('Figure_1_AB_Overview_Metadata_Distribution.svg', format='svg', bbox_inches='tight')
print('Figure A and B generated and saved successfully as 600 DPI PNG, PDF, and SVG formats with bold text and updated host nomenclature.')