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
print(" [Figure 9] High-Resolution (600 DPI) Multi-Format Export Script")
print("=" * 70)

# 2. 연도별 유전체 완전성(Completeness_Flag) 백분율 집계
completeness_cross = pd.crosstab(df['Collection_Year'], df['Completeness_Flag'], normalize='index') * 100
all_years = range(int(df['Collection_Year'].min()), int(df['Collection_Year'].max()) + 1)
completeness_cross = completeness_cross.reindex(all_years, fill_value=0.0)

# 3. Figure 9: 연도별 완전 유전체 vs 부분 유전자 조각 비율 100% 누적 바차트 렌더링
fig, ax = plt.subplots(figsize=(18, 9), dpi=600)
fig.patch.set_facecolor('#ffffff')

years = completeness_cross.index
complete_pct = completeness_cross['Complete Genome'] if 'Complete Genome' in completeness_cross.columns else pd.Series(0, index=years)
partial_pct = completeness_cross['Partial Gene Fragment'] if 'Partial Gene Fragment' in completeness_cross.columns else pd.Series(0, index=years)

ax.bar(years, complete_pct, label='Complete Genome', color='#3b82f6', edgecolor='none', alpha=0.85, width=0.8)
ax.bar(years, partial_pct, bottom=complete_pct, label='Partial Gene Fragment', color='#f59e0b', edgecolor='none', alpha=0.85, width=0.8)

ax.axvline(x=2015, color='#7c3aed', linestyle='--', linewidth=3.0, label='HTS Mass Surveillance Introduction (2015)')

# 상단 타이틀은 완전히 제거하고 축 레이블과 눈금 폰트 크기/두께 강화
ax.set_xlabel('Collection Year', fontsize=14, fontweight='bold', labelpad=12)
ax.set_ylabel('Proportion (%)', fontsize=14, fontweight='bold', labelpad=12)
ax.set_ylim(0, 105)

ax.tick_params(axis='both', labelsize=12, width=2)
for label in ax.get_yticklabels():
    label.set_fontweight('bold')
for label in ax.get_xticklabels():
    label.set_fontweight('bold')

ax.legend(loc='upper right', fontsize=12, frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1', borderpad=1.2)
ax.grid(True, linestyle=':', alpha=0.6, linewidth=1.2)

plt.tight_layout()

# 4. 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장 (파일명에 figure9 포함)
base_filename = "Figure_9_Completeness_Proportion"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)