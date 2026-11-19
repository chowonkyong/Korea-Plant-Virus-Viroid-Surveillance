import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 스타일 및 글씨 두께 설정
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

fig, axes = plt.subplots(1, 2, figsize=(20, 8))

# 실제 마스터 데이터셋 로드
df = pd.read_csv('korea_plant_viruses_core_perfect.csv')

# --- Panel A: Taxonomic Distribution ---
family_counts = df['Family'].value_counts().head(10)
y_pos = range(len(family_counts))
colors = sns.color_palette("husl", len(family_counts))

bars = axes[0].barh(y_pos, family_counts.values, color=colors, edgecolor='black', alpha=0.85)
axes[0].set_yticks(y_pos)
axes[0].set_yticklabels(family_counts.index, fontsize=12, fontweight='bold', fontstyle='italic')
axes[0].invert_yaxis()

axes[0].set_title('A. Taxonomic Distribution of Major Viral Families', fontsize=15, fontweight='bold', pad=15)
axes[0].set_xlabel('Number of Records (n)', fontsize=13, fontweight='bold')
axes[0].set_ylabel('Viral Family', fontsize=13, fontweight='bold')
axes[0].tick_params(axis='both', which='major', labelsize=11)
for tick in axes[0].get_xticklabels():
    tick.set_fontweight('bold')

for bar in bars:
    width = bar.get_width()
    axes[0].text(width + 15, bar.get_y() + bar.get_height()/2, f'n = {width}', 
                 va='center', ha='left', fontsize=11, fontweight='bold')

# --- Panel B: Sequence Length Distribution by Genome Type ---
g_order = ['ssRNA (positive)', 'ssRNA (negative)', 'dsDNA', 'dsRNA', 'ssDNA', 'viroids']
available_orders = [g for g in g_order if g in df['Genome_Type'].unique()]

sns.violinplot(data=df, x='Length', y='Genome_Type', hue='Genome_Type', order=available_orders, 
               palette='Set2', ax=axes[1], inner='quartile', cut=0, legend=False)

axes[1].set_title('B. Sequence Length Distribution by Genome Type', fontsize=15, fontweight='bold', pad=15)
axes[1].set_xlabel('Sequence Length (bp)', fontsize=13, fontweight='bold')
axes[1].set_ylabel('Genome Type', fontsize=13, fontweight='bold')
axes[1].tick_params(axis='both', which='major', labelsize=11)
for tick in axes[1].get_yticklabels():
    tick.set_fontweight('bold')
for tick in axes[1].get_xticklabels():
    tick.set_fontweight('bold')

plt.tight_layout()

# 고해상도 600 DPI PNG 및 논문 출판용 벡터 포맷(PDF, SVG)으로 저장
plt.savefig('Figure_2_Combined_Taxonomy_and_Length.png', dpi=600, bbox_inches='tight')
plt.savefig('Figure_2_Combined_Taxonomy_and_Length.pdf', format='pdf', bbox_inches='tight')
plt.savefig('Figure_2_Combined_Taxonomy_and_Length.svg', format='svg', bbox_inches='tight')
print('Figure 2 generated and saved successfully as 600 DPI PNG, PDF, and SVG formats with italicized family names and bold text.')
