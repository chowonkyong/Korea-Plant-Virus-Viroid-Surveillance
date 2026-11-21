import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 1. 실제 핵심 데이터셋 로드 및 무결성 검사
csv_path = 'korea_plant_viruses_core_perfect.csv'
if not os.path.exists(csv_path):
    for root, dirs, files in os.walk('.'):
        if 'korea_plant_viruses_core_perfect.csv' in files:
            csv_path = os.path.join(root, 'korea_plant_viruses_core_perfect.csv')
            break

if not os.path.exists(csv_path):
    raise FileNotFoundError("Error: 'korea_plant_viruses_core_perfect.csv' file could not be found in the directory.")

df = pd.read_csv(csv_path)

print("=" * 70)
print(" [Core Dataset Exploratory Data Analysis & Statistical Summary]")
print("=" * 70)
print(f"Total Records (Rows): {df.shape[0]}")
print(f"Total Features (Columns): {df.shape[1]}")
print("\n[Column Information & Data Types]")
print(df.dtypes)
print("\n[Missing Values Check]")
print(df.isnull().sum())
print("=" * 70)

# 2. 연도별 검출 건수 통계 및 추이 분석
if 'Collection_Year' in df.columns:
    yearly_counts = df['Collection_Year'].value_counts().sort_index()
    print("\n[Annual Detection Counts Summary]")
    print(yearly_counts.to_string())
    print("=" * 70)

# 3. 주요 바이러스/병원체(Organism) 상위 빈도 분석
if 'Organism' in df.columns:
    top_organisms = df['Organism'].value_counts().head(15)
    print("\n[Top 15 Most Frequently Detected Organisms]")
    print(top_organisms.to_string())
    print("=" * 70)

# 4. 유전체 유형(Genome Type) 분포 분석 (컬럼이 없을 경우 자동 분류 파생)
if 'Genome_Type' not in df.columns and 'Organism' in df.columns:
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

if 'Genome_Type' in df.columns:
    genome_counts = df['Genome_Type'].value_counts()
    print("\n[Genome Type Distribution Summary]")
    print(genome_counts.to_string())
    print("=" * 70)

# 5. 종합 데이터 분석 요약 시각화 패널 생성 (EDA Summary Figure)
fig, axes = plt.subplots(1, 2, figsize=(16, 6), dpi=300)
fig.patch.set_facecolor('#ffffff')

# Chart 1: 연도별 검출 건수 트렌드
if 'Collection_Year' in df.columns:
    yearly_counts.plot(kind='bar', ax=axes[0], color='#2563eb', edgecolor='#1d4ed8', alpha=0.85)
    axes[0].set_title('Annual Plant Virus Detection Trends', fontsize=13, fontweight='bold', pad=12)
    axes[0].set_xlabel('Collection Year', fontsize=11, fontweight='bold')
    axes[0].set_ylabel('Detection Frequency', fontsize=11, fontweight='bold')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].grid(True, linestyle=':', alpha=0.5)

# Chart 2: 유전체 유형별 비율 파이 차트 또는 바 차트
if 'Genome_Type' in df.columns:
    genome_counts.plot(kind='barh', ax=axes[1], color='#059669', edgecolor='#047857', alpha=0.85)
    axes[1].set_title('Detection Frequency by Genome Type', fontsize=13, fontweight='bold', pad=12)
    axes[1].set_xlabel('Detection Count', fontsize=11, fontweight='bold')
    axes[1].set_ylabel('Genome Type', fontsize=11, fontweight='bold')
    axes[1].grid(True, linestyle=':', alpha=0.5)

plt.tight_layout()
output_chart = "Dataset_EDA_Summary_Report.png"
plt.savefig(output_chart, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
plt.close(fig)
print(f"\nSuccessfully generated EDA summary visualization: {output_chart}")
print("=" * 70)
