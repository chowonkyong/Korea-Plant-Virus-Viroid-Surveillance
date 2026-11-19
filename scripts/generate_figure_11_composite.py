import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# 1. 데이터 정의 (도시 및 기관 상위 목록)
city_counts = pd.Series({
    'Choong-Nam': 122,
    'Suwon': 32,
    'Jincheon': 23,
    'Andong': 23,
    'Gwangju': 17,
    'Hongsan': 16,
    'Anseong': 16,
    'Nonsan': 15,
    'Gyeongju': 14,
    'Daegu': 14
}).sort_values(ascending=True)

inst_counts = pd.Series({
    'Chonnam Nat. Univ.': 532,
    'Seoul Nat. Univ.': 451,
    'Kyungpook Nat. Univ.': 152,
    'Andong Nat. Univ.': 148,
    'NIHHS (RDA)': 140,
    'NAAS (RDA)': 135,
    'Chung-Ang Univ.': 122,
    'Crop Science (RDA)': 55,
    'Jeju Nat. Univ.': 46,
    'Kangwon Nat. Univ.': 35
}).sort_values(ascending=True)

# 2. 2단 패널 컴포지트 캔버스 설정 (600 DPI)
fig, axes = plt.subplots(1, 2, figsize=(22, 10), dpi=600)
fig.patch.set_facecolor('#ffffff')

# 패널 A: 도시 및 지역별 검출 빈도
city_counts.plot(kind='barh', ax=axes[0], color='#0284c7', edgecolor='none', width=0.7)
axes[0].set_title('A. Top South Korean Cities & Regions in Plant Virus Surveillance', fontsize=13, fontweight='bold', pad=15)
axes[0].set_xlabel('Detection Frequency (Count)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Region / City', fontsize=12, fontweight='bold')
axes[0].tick_params(axis='both', labelsize=11)
for label in axes[0].get_yticklabels():
    label.set_fontweight('bold')
for label in axes[0].get_xticklabels():
    label.set_fontweight('bold')
axes[0].grid(axis='x', linestyle=':', alpha=0.6, linewidth=1.2)

# 패널 B: 연구기관 및 대학별 기탁 빈도
inst_counts.plot(kind='barh', ax=axes[1], color='#0d9488', edgecolor='none', width=0.7)
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

# 3. 600 DPI 고해상도 PNG, SVG, PDF 포맷 동시 저장
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"Figure_11_Composite.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)
