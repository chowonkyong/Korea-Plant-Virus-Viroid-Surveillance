import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 1. 최종 정제된 데이터 로드
df = pd.read_csv('korea_plant_viruses_final_clean.csv')

# 2. 상위 10개 지역 및 기관 집계 (Unspecified 제외 혹은 포함 선택 가능)
top_locations = df[~df['Location_Clean'].isin(['Unspecified', 'missing: third party data'])]['Location_Clean'].value_counts().head(10).sort_values(ascending=True)
top_institutions = df[~df['Institution_Clean'].isin(['Unspecified', 'General GenBank Submission'])]['Institution_Clean'].value_counts().head(10).sort_values(ascending=True)

# 3. Figure 11 컴포지트 그래프 생성 (600 DPI)
fig, axes = plt.subplots(1, 2, figsize=(22, 10), dpi=600)
fig.patch.set_facecolor('#ffffff')

# Panel A: Cities & Regions
top_locations.plot(kind='barh', ax=axes[0], color='#0284c7', edgecolor='none', width=0.7)
axes[0].set_title('A. Top South Korean Cities & Regions in Plant Virus Surveillance', fontsize=13, fontweight='bold', pad=15)
axes[0].set_xlabel('Detection Frequency (Count)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Region / City', fontsize=12, fontweight='bold')
axes[0].tick_params(axis='both', labelsize=11)
for label in axes[0].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[0].get_xticklabels():
    label.set_fontweight('bold')
axes[0].grid(axis='x', linestyle=':', alpha=0.6, linewidth=1.2)

# Panel B: Institutions & Universities
top_institutions.plot(kind='barh', ax=axes[1], color='#0d9488', edgecolor='none', width=0.7)
axes[1].set_title('B. Top Research Institutions & Universities Contributing Sequences', fontsize=13, fontweight='bold', pad=15)
axes[1].set_xlabel('Sequence Submission / Detection Count', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Research Institution / University', fontsize=12, fontweight='bold')
axes[1].tick_params(axis='both', labelsize=11)
for label in axes[1].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[1].get_xticklabels():
    label.set_fontweight('bold')
axes[1].grid(axis='x', linestyle=':', alpha=0.6, linewidth=1.2)

plt.tight_layout()

# 저장
for ext in ['png', 'pdf', 'svg']:
    plt.savefig(f'Figure_11_Composite_Updated.{ext}', dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')

plt.close(fig)
print("Successfully regenerated Figure 11 using standardized clean data!")
