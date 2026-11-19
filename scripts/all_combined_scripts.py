

# ==========================================
# FILE: add_completeness.py
# ==========================================

import pandas as pd
from Bio import SeqIO
import os

gbk_path = "../data/korea_plant_viruses_clean.gbk"
master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(gbk_path):
    print(f"[Error] GenBank 파일을 찾을 수 없습니다: {gbk_path}")
    exit()

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("GenBank 파일에서 서열 설명(Description) 및 completeness 정보 추출 중...")
completeness_dict = {}

for record in SeqIO.parse(gbk_path, "genbank"):
    acc_full = record.id
    acc_clean = acc_full.split('.')[0]
    description = record.description.lower()
    
    if "complete genome" in description or "complete sequence" in description:
        status = "Complete Genome"
    elif "partial" in description or "fragment" in description or "gene" in description:
        status = "Partial / Gene Fragment"
    else:
        status = "Unspecified / Other"
        
    completeness_dict[acc_full] = status
    completeness_dict[acc_clean] = status

df = pd.read_csv(master_csv)

def get_completeness(row):
    acc = str(row['Accession'])
    acc_clean = acc.split('.')[0]
    if acc in completeness_dict:
        return completeness_dict[acc]
    elif acc_clean in completeness_dict:
        return completeness_dict[acc_clean]
    return "Unknown"

print("마스터 메타데이터에 Complete_Genome_Status 병합 중...")
df['Complete_Genome_Status'] = df.apply(get_completeness, axis=1)

df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 전장 유전체 여부 정보가 마스터 파일({master_csv})에 성공적으로 통합되었습니다!")
print("\n[분포 현황 요약]")
print(df['Complete_Genome_Status'].value_counts())

print("\n[상위 5개 결과 미리보기]")
cols_to_show = [c for c in ['Accession', 'Viral_Species', 'Sequence_Length', 'Complete_Genome_Status'] if c in df.columns]
print(df[cols_to_show].head())


# ==========================================
# FILE: add_genome_segmentation.py
# ==========================================

import pandas as pd
import os

csv_path = 'plant_viruses_metadata_enriched.csv'
if not os.path.exists(csv_path):
    csv_path = 'plant_viruses_metadata_standardized.csv'

df = pd.read_csv(csv_path)

def assign_segmentation(row):
    tax = str(row['Taxonomy_Lineage']).lower()
    species = str(row['Viral_Species']).lower()
    
    if 'viroid' in tax or 'viroid' in species:
        return 'Not Applicable (Viroid)'
    elif 'germiniviridae' in tax:
        if 'begomovirus' in tax:
            return 'Bipartite / Monopartite (Geminiviridae)'
        else:
            return 'Monopartite (Geminiviridae)'
    elif any(fam in tax for fam in ['bromoviridae', 'cucumoviridae', 'secoviridae', 'partitiviridae', 'comoviridae']):
        return 'Multipartite / Bipartite'
    elif any(fam in tax for fam in ['potyviridae', 'tombusviridae', 'betaflexiviridae', 'alphaflexiviridae', 'virgaviridae', 'closteroviridae', 'tymoviridae', 'capillovirus', 'trichovirus', 'foveavirus']):
        return 'Monopartite'
    elif any(fam in tax for fam in ['rhabdoviridae', 'tospoviridae', 'bunyaviridae', 'phenuiviridae']):
        return 'Multipartite / Negative-sense'
    else:
        return 'Monopartite (Default / Other RNA)'

df['Genome_Segmentation'] = df.apply(assign_segmentation, axis=1)

print("Updated Segmentation Counts:")
print(df['Genome_Segmentation'].value_counts())

# 마스터 파일 덮어쓰기 저장
df.to_csv('plant_viruses_metadata_enriched.csv', index=False, encoding='utf-8-sig')
print("Master file updated with Genome_Segmentation successfully.")


# ==========================================
# FILE: add_genome_type.py
# ==========================================

import pandas as pd
import os
import glob

# 파일이 있는 경로를 자동으로 탐색
possible_paths = [
    "plant_viruses_metadata_standardized.csv",
    "../plant_viruses_metadata_standardized.csv",
    "../data/plant_viruses_metadata_standardized.csv",
    "/mnt/c/virus_korea_project/data/plant_viruses_metadata_standardized.csv"
]

input_csv = None
for path in possible_paths:
    if os.path.exists(path):
        input_csv = path
        break

# 만약 위 경로들에 없다면 전체 시스템 검색
if not input_csv:
    found_files = glob.glob('**/plant_viruses_metadata_standardized.csv', recursive=True)
    if found_files:
        input_csv = found_files[0]

if not input_csv or not os.path.exists(input_csv):
    print("[Error] 'plant_viruses_metadata_standardized.csv' 파일을 찾을 수 없습니다. 파일 위치를 확인해 주세요.")
    exit()

print(f"파일을 성공적으로 찾았습니다: {os.path.abspath(input_csv)}")
df = pd.read_csv(input_csv)

# 유전체 유형(Genome Type) 자동 분류 함수
def assign_genome_type(row):
    tax = str(row['Taxonomy_Lineage']).lower()
    species = str(row['Viral_Species']).lower()
    
    if 'viroid' in tax or 'viroid' in species:
        return 'Viroid'
    elif 'germiniviridae' in tax or 'nanoviridae' in tax:
        return 'ssDNA'
    elif 'caulimoviridae' in tax:
        return 'dsDNA'
    elif 'partitiviridae' in tax:
        return 'dsRNA'
    elif 'rhabdoviridae' in tax or 'tospoviridae' in tax:
        return '-ssRNA'
    elif any(fam in tax for fam in ['potyviridae', 'bromoviridae', 'betaflexiviridae', 'alphaflexiviridae', 'secoviridae', 'virgaviridae', 'tombusviridae', 'closteroviridae', 'tymoviridae']):
        return '+ssRNA'
    else:
        return 'Other / Unclassified'

print("유전체 유형(Genome Type) 매핑 중...")
df['Genome_Type'] = df.apply(assign_genome_type, axis=1)

# 결과 확인 출력
print("\n[ 유전체 유형별 분포 현황 ]")
print(df['Genome_Type'].value_counts())

output_csv = "plant_viruses_metadata_enriched.csv"
df.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"\n유전체 유형이 추가된 파일이 성공적으로 저장되었습니다: {os.path.abspath(output_csv)}")


# ==========================================
# FILE: add_geographic_region.py
# ==========================================

import pandas as pd
from Bio import SeqIO
import os

gbk_path = "../data/korea_plant_viruses_clean.gbk"
master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(gbk_path):
    print(f"[Error] GenBank 파일을 찾을 수 없습니다: {gbk_path}")
    exit()

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("GenBank 파일에서 지리적 출처(Country / Province) 정보 추출 중...")
region_dict = {}

for record in SeqIO.parse(gbk_path, "genbank"):
    acc_full = record.id
    acc_clean = acc_full.split('.')[0]
    
    country_info = "South Korea (Unspecified Province)"
    
    for feature in record.features:
        if feature.type == "source":
            if "country" in feature.qualifiers:
                country_info = feature.qualifiers["country"][0]
                
    region_dict[acc_full] = country_info
    region_dict[acc_clean] = country_info

df = pd.read_csv(master_csv)

def get_region(row):
    acc = str(row['Accession'])
    acc_clean = acc.split('.')[0]
    if acc in region_dict:
        return region_dict[acc]
    elif acc_clean in region_dict:
        return region_dict[acc_clean]
    return "South Korea (Unspecified)"

print("마스터 메타데이터에 Geographic_Region 컬럼 병합 중...")
df['Geographic_Region'] = df.apply(get_region, axis=1)

df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 지리적 세부 출원/채집 지역 정보가 마스터 파일({master_csv})에 성공적으로 통합되었습니다!")
print("\n[주요 지역별 분포 현황]")
print(df['Geographic_Region'].value_counts().head(10))

print("\n[상위 5개 결과 미리보기]")
cols = [c for c in ['Accession', 'Viral_Species', 'Geographic_Region'] if c in df.columns]
print(df[cols].head())


# ==========================================
# FILE: add_pathogenicity_and_analyze.py
# ==========================================

import pandas as pd
import os

master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("마스터 메타데이터 로드 중...")
df = pd.read_csv(master_csv)

# 5. 증상 및 병원성 유형(Pathogenicity / Symptoms) 자동 매핑 함수
def assign_symptoms(row):
    tax = str(row['Taxonomy_Lineage']).lower()
    species = str(row['Viral_Species']).lower()
    
    if 'viroid' in tax or 'viroid' in species:
        return 'Stunting / Chlorosis / Asymptomatic'
    elif 'potyviridae' in tax or 'potyvirus' in tax:
        return 'Mosaic / Mottling / Leaf Distortion'
    elif 'germiniviridae' in tax or 'begomovirus' in tax:
        return 'Yellow Leaf Curl / Stunting / Leaf Curling'
    elif 'tospoviridae' in tax or 'orthotospovirus' in tax:
        return 'Necrosis / Ringspot / Wilt'
    elif 'cucumoviridae' in tax or 'bromoviridae' in tax:
        return 'Mosaic / Stunting / Ring spot'
    elif 'rhabdoviridae' in tax:
        return 'Stunting / Enation / Chlorosis'
    elif 'secoviridae' in tax or 'nepovirus' in tax:
        return 'Ring spot / Mosaic / Stunting'
    elif 'partitiviridae' in tax or 'endornaviridae' in tax:
        return 'Latent / Asymptomatic (Cryptic)'
    elif 'tobamovirus' in species or 'virgaviridae' in tax:
        return 'Mosaic / Necrosis / Leaf Blister'
    elif 'alphaflexiviridae' in tax or 'betaflexiviridae' in tax:
        return 'Mosaic / Chlorotic Spot / Latent'
    else:
        return 'Mosaic / General Chlorosis / Stunting'

print("바이러스 과/속/종 분류군을 바탕으로 병원성 및 증상(Pathogenicity_Symptoms) 정보 매핑 중...")
df['Pathogenicity_Symptoms'] = df.apply(assign_symptoms, axis=1)

# 마스터 파일 덮어쓰기 저장 (총 13개 컬럼 완성)
df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 증상 및 병원성 유형 정보가 마스터 파일({master_csv})에 성공적으로 통합되었습니다!")
print(f"최종 마스터 데이터 규격: 총 {df.shape[0]}행, {df.shape[1]}개 컬럼")

print("\n[보유 중인 13개 핵심 컬럼 리스트]")
print(df.columns.tolist())

print("\n[병원성 및 증상 유형별 분포 현황]")
print(df['Pathogenicity_Symptoms'].value_counts())


# ==========================================
# FILE: add_sequence_length.py
# ==========================================

import pandas as pd
from Bio import SeqIO
import os

# 파일 경로 설정 (환경에 맞게 수정 가능)
gbk_path = "../data/korea_plant_viruses_clean.gbk"  # GenBank 파일 경로
master_csv = "plant_viruses_metadata_enriched.csv"     # 마스터 메타데이터 CSV

if not os.path.exists(gbk_path):
    print(f"[Error] GenBank 파일을 찾을 수 없습니다: {gbk_path}")
    exit()

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("GenBank 파일에서 각 Accession별 유전체 크기(Sequence Length) 추출 중...")
length_dict = {}
for record in SeqIO.parse(gbk_path, "genbank"):
    # 버전 포함(Accession.version)과 순수 Accession 모두 매핑
    acc_full = record.id
    acc_clean = acc_full.split('.')[0]
    seq_len = len(record.seq)
    
    length_dict[acc_full] = seq_len
    length_dict[acc_clean] = seq_len

print(f"총 {len(length_dict)}개의 서열 길이를 파싱했습니다.")

# 마스터 메타데이터 로드
df = pd.read_csv(master_csv)

# Sequence_Length 매핑 함수
def get_length(row):
    acc = str(row['Accession'])
    acc_clean = acc.split('.')[0]
    if acc in length_dict:
        return length_dict[acc]
    elif acc_clean in length_dict:
        return length_dict[acc_clean]
    return None

print("마스터 메타데이터에 유전체 사이즈(Sequence_Length) 병합 중...")
df['Sequence_Length'] = df.apply(get_length, axis=1)

# 마스터 파일 덮어쓰기 저장 (통합 유지)
df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 유전체 사이즈 정보가 마스터 파일({master_csv})에 성공적으로 통합되었습니다!")
print(f"누락된 사이즈 데이터 수: {df['Sequence_Length'].isnull().sum()}건")
print("\n[상위 5개 결과 미리보기]")
print(df[['Accession', 'Viral_Species', 'Genome_Type', 'Sequence_Length']].head())


# ==========================================
# FILE: add_vector_mode.py
# ==========================================

import pandas as pd
import os

master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("마스터 메타데이터 로드 중...")
df = pd.read_csv(master_csv)

def assign_vector_mode(row):
    tax = str(row['Taxonomy_Lineage']).lower()
    species = str(row['Viral_Species']).lower()
    
    if 'viroid' in tax or 'viroid' in species:
        return 'Not Applicable (Viroid / Mechanical/Seed)'
    elif 'potyviridae' in tax or 'potyvirus' in tax:
        return 'Insect: Aphid (Non-persistent)'
    elif 'germiniviridae' in tax or 'begomovirus' in tax:
        return 'Insect: Whitefly (Bemisia tabaci)'
    elif 'tospoviridae' in tax or 'orthotospovirus' in tax:
        return 'Insect: Thrips (Frankliniella / Thrips)'
    elif 'rhabdoviridae' in tax:
        return 'Insect: Leafhopper / Planthopper'
    elif 'cucumoviridae' in tax or 'bromoviridae' in tax:
        return 'Insect: Aphid / Mechanical'
    elif 'secoviridae' in tax or 'nepovirus' in tax:
        return 'Nematode / Mechanical'
    elif 'partitiviridae' in tax or 'endornaviridae' in tax:
        return 'Fungal / Seed / Pollen (Vertical)'
    elif 'tobamovirus' in species or 'virgaviridae' in tax:
        return 'Mechanical / Seed (Contact transmission)'
    elif 'alphaflexiviridae' in tax or 'betaflexiviridae' in tax:
        return 'Mechanical / Aphid / Grafting'
    else:
        return 'Mechanical / Seed / Other Insect'

print("바이러스 과/속/종 분류군을 기반으로 매개체(Vector & Transmission Mode) 정보 매핑 중...")
df['Vector_Transmission_Mode'] = df.apply(assign_vector_mode, axis=1)

# 마스터 파일 덮어쓰기 저장
df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 매개체 및 전파 방식 정보가 마스터 파일({master_csv})에 성공적으로 통합되었습니다!")
print("\n[매개체 유형별 분포 현황]")
print(df['Vector_Transmission_Mode'].value_counts())

print("\n[상위 5개 결과 미리보기]")
cols = [c for c in ['Accession', 'Viral_Species', 'Genome_Type', 'Vector_Transmission_Mode'] if c in df.columns]
print(df[cols].head())


# ==========================================
# FILE: analyze_architecture.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*50)
print(" 🧬 [분석 1] 유전체 구조별 서열 길이 통계")
print("="*50)
print(df.groupby(['Genome_Type', 'Genome_Segmentation'])['Sequence_Length'].agg(
    Count='count', Mean='mean', Median='median', Min='min', Max='max'
).round(1))

# 박스플롯 시각화 생성 (이상치 제외 가독성 극대화)
plt.figure(figsize=(12, 6))
sns.boxplot(x='Genome_Type', y='Sequence_Length', hue='Genome_Segmentation', data=df, palette='Set2', showfliers=False)
plt.title('Sequence Length Distribution by Genome Type and Architecture', fontsize=14, fontweight='bold')
plt.xlabel('Genome Type', fontsize=12)
plt.ylabel('Sequence Length (bp)', fontsize=12)
plt.xticks(rotation=30)
plt.legend(title='Genome Segmentation', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('genomic_architecture_box.png', dpi=300)
plt.close()
print("\n[완료] 시각화 저장: genomic_architecture_box.png")


# ==========================================
# FILE: analyze_detailed_viruses_and_viroids.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for detailed virus shift and viroid analysis...")
df = pd.read_csv(csv_file)

# 연도 추출 로직
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        y = int(match.group(1))
        if 2000 <= y <= 2026:
            return y
    return None

df['Year'] = df['Collection_Date'].apply(extract_year)
df_clean = df.dropna(subset=['Year', 'Organism']).copy()

# 1. 'Other' 비중을 줄이고 Top 10 바이러스 종으로 세분화
top_10_viruses = df_clean['Organism'].value_counts().head(10).index.tolist()

def group_virus_detailed(org):
    org_str = str(org)
    if "viroid" in org_str.lower():
        return "Viroids (Group)"
    elif org_str in top_10_viruses:
        return org_str
    else:
        return "Other Minor Viruses"

df_clean['Detailed_Group'] = df_clean['Organism'].apply(group_virus_detailed)

# 연도별 x 상세 바이러스 그룹 교차 테이블 (비율 계산)
yearly_detailed_xtab = pd.crosstab(df_clean['Year'], df_clean['Detailed_Group'], normalize='index') * 100

# CSV 저장 (세분화된 연도별 변동)
shift_path = os.path.join(results_dir, "korea_detailed_virus_yearly_shift.csv")
yearly_detailed_xtab.to_csv(shift_path)
print(f"-> Saved detailed yearly shift CSV: {shift_path}")

# 시각화 1: 세분화된 우점종 역학 변동 스택 바차트 (Top 10 + Viroids)
plt.figure(figsize=(14, 8))
ax = yearly_detailed_xtab.plot(
    kind='bar', 
    stacked=True, 
    figsize=(14, 8), 
    colormap='tab20', 
    edgecolor='none'
)
plt.title('Detailed Shift in Dominant Plant Viruses & Viroids Over Time in Korea (2000–2026)', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Proportion of Reported Sequences (%)', fontsize=12)
plt.legend(title='Virus Species / Viroids', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()

fig_path_1 = os.path.join(fig_dir, "korea_detailed_dominant_virus_shift.png")
plt.savefig(fig_path_1, dpi=300, bbox_inches='tight')
plt.close()
print(f"-> Saved high-res detailed shift figure: {fig_path_1}")

# 2. 비로이드(Viroids) 전용 심층 분석
df_viroids = df_clean[df_clean['Organism'].str.contains("viroid", case=False, na=False)].copy()

if len(df_viroids) > 0:
    viroid_counts = df_viroids['Organism'].value_counts().reset_index()
    viroid_counts.columns = ['Viroid_Species', 'Count']
    viroid_summary_path = os.path.join(results_dir, "korea_viroids_summary.csv")
    viroid_counts.to_csv(viroid_summary_path, index=False)
    print(f"-> Saved viroids summary CSV: {viroid_summary_path}")
    print("\nViroids Found in Korea Preview:")
    print(viroid_counts.head(10))

    # 비로이드 전용 바차트 생성
    plt.figure(figsize=(10, 5))
    sns.barplot(
        data=viroid_counts.head(8), 
        x='Count', 
        y='Viroid_Species', 
        hue='Viroid_Species',
        palette='YlOrRd_r',
        legend=False
    )
    plt.title('Detected Viroid Species in South Korea', fontsize=13, fontweight='bold', pad=15)
    plt.xlabel('Number of Sequences', fontsize=12)
    plt.ylabel('Viroid Species', fontsize=12)
    plt.tight_layout()

    fig_path_2 = os.path.join(fig_dir, "korea_viroids_top_species.png")
    plt.savefig(fig_path_2, dpi=300)
    plt.close()
    print(f"-> Saved viroids top species figure: {fig_path_2}")
else:
    print("[Notice] No specific viroid records found matching the keyword.")

print("All detailed virus and viroid analyses completed successfully!")


# ==========================================
# FILE: analyze_diagnostic_cohort.py
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

# 1. 수집 연도(Collection_Year) 기준 진단 기술 코호트 구분
df['Collection_Year'] = pd.to_numeric(df['Collection_Year'], errors='coerce')
df['Diagnostic_Cohort'] = pd.cut(
    df['Collection_Year'],
    bins=[1980, 2014, 2026],
    labels=['Traditional Era (~2014)', 'NGS Metagenomics Era (2015~)']
)

# 2. 전장 유전체 확보 여부 이진형 변환
df['Is_Complete_Genome'] = df['Complete_Genome_Status'].apply(lambda x: 1 if str(x).strip() == 'Complete Genome' else 0)
cohort_df = df.dropna(subset=['Diagnostic_Cohort', 'Is_Complete_Genome']).copy()

# 3. 코호트별 성공률 집계
cohort_summary = cohort_df.groupby('Diagnostic_Cohort', observed=False)['Is_Complete_Genome'].agg(['mean', 'count']).reset_index()
cohort_summary['Success_Rate_%'] = cohort_summary['mean'] * 100

print("\n--- [코호트별 전장 유전체 확보율 요약] ---")
print(cohort_summary[['Diagnostic_Cohort', 'count', 'Success_Rate_%']])

# 4. 시각화 (Seaborn 경고 해결을 위해 hue 및 legend 옵션 추가)
plt.figure(figsize=(10, 6))

ax = sns.barplot(
    data=cohort_summary, x='Diagnostic_Cohort', y='Success_Rate_%', 
    hue='Diagnostic_Cohort', palette=['#4682B4', '#2E8B57'], edgecolor='black', legend=False
)

plt.title('Cohort Analysis: Complete Genome Recovery Rate by Diagnostic Era', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Diagnostic Technology Cohort', fontsize=12)
plt.ylabel('Complete Genome Recovery Rate (%)', fontsize=12)
plt.ylim(0, 105)

# 바 상단 데이터 레이블 부착
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(f"{height:.1f}%",
                    (p.get_x() + p.get_width() / 2., height),
                    ha='center', va='bottom', fontsize=11, fontweight='bold',
                    xytext=(0, 3), textcoords='offset points')

plt.tight_layout()
output_fig = "diagnostic_cohort_analysis.png"
plt.savefig(output_fig, dpi=300)
plt.close()

print(f"\n[완료] 정제된 코호트 분석 그래프 저장 성공: {output_fig}")


# ==========================================
# FILE: analyze_dominant_virus_shift.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for dominant virus shift over time analysis...")
df = pd.read_csv(csv_file)

# 연도 추출 로직
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        y = int(match.group(1))
        if 2000 <= y <= 2026: # 최근 역학 변동이 집중된 2000년 이후 대상
            return y
    return None

df['Year'] = df['Collection_Date'].apply(extract_year)
df_clean = df.dropna(subset=['Year', 'Organism']).copy()

# 상위 5~6개 주요 바이러스 종 선정 후 나머지는 'Other'로 묶기
top_viruses = df_clean['Organism'].value_counts().head(6).index.tolist()
def group_virus(org):
    if org in top_viruses:
        return org
    else:
        return "Other Viruses"

df_clean['Virus_Group'] = df_clean['Organism'].apply(group_virus)

# 연도별 x 바이러스 그룹 교차 테이블 (비율 계산)
yearly_virus_xtab = pd.crosstab(df_clean['Year'], df_clean['Virus_Group'], normalize='index') * 100

# CSV 저장
shift_path = os.path.join(results_dir, "korea_dominant_virus_yearly_shift.csv")
yearly_virus_xtab.to_csv(shift_path)
print(f"-> Saved yearly shift summary CSV: {shift_path}")

print("\nYearly Dominant Virus Shift Preview (Percentage):")
print(yearly_virus_xtab.tail(10))

# 시각화: 스택드 바차트 (Stacked Bar Chart)
plt.figure(figsize=(14, 7))
ax = yearly_virus_xtab.plot(
    kind='bar', 
    stacked=True, 
    figsize=(14, 7), 
    colormap='tab10', 
    edgecolor='none'
)
plt.title('Shift in Dominant Plant Virus Species Over Time in Korea (2000–2026)', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Proportion of Reported Sequences (%)', fontsize=12)
plt.legend(title='Virus Species / Group', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.xticks(rotation=45)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_dominant_virus_shift.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"-> Saved high-res yearly shift figure: {fig_path}")


# ==========================================
# FILE: analyze_enriched_dataset.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

input_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(input_csv):
    print(f"[Error] File not found: {input_csv}")
    exit()

print("="*60)
print(" LOADING AND ANALYZING PLANT VIRUS METADATA")
print("="*60)

df = pd.read_csv(input_csv)

# 1. 기본 정보 요약
total_records = len(df)
print(f"Total Records: {total_records}")

# 2. 유전체 유형(Genome Type) 통계
genome_counts = df['Genome_Type'].value_counts()
print("\n[1] Genome Type Distribution:")
for g_type, count in genome_counts.items():
    pct = (count / total_records) * 100
    print(f"  - {g_type:<20}: {count}건 ({pct:.2f}%)")

# 3. 주요 기주(Host Plant) 통계 ('Not Specified' 제외)
df_known_hosts = df[df['Host_Plant'] != 'Not Specified']
host_counts = df_known_hosts['Host_Plant'].value_counts().head(10)
print("\n[2] Top 10 Plant Hosts (Standardized Scientific Names):")
for host, count in host_counts.items():
    print(f"  - {host:<30}: {count}건")

# 4. 주요 바이러스 종(Viral Species) 통계
species_counts = df['Viral_Species'].value_counts().head(10)
print("\n[3] Top 10 Viral & Viroid Species:")
for species, count in species_counts.items():
    print(f"  - {species:<40}: {count}건")

# 5. 시각화 생성 (유전체 유형별 바차트)
plt.figure(figsize=(10, 6), dpi=300)
sns.countplot(data=df, y='Genome_Type', order=genome_counts.index, palette='viridis')
plt.title('Distribution of Plant Pathogens by Genome Type', fontsize=14, weight='bold')
plt.xlabel('Number of Records', fontsize=12)
plt.ylabel('Genome Type', fontsize=12)
plt.grid(axis='x', linestyle='--', alpha=0.7)

output_img = 'genome_type_distribution.png'
plt.tight_layout()
plt.savefig(output_img, dpi=300)
plt.close()

print("="*60)
print(f"분석 완료! 시각화 차트가 '{os.path.abspath(output_img)}'에 저장되었습니다.")
print("="*60)


# ==========================================
# FILE: analyze_families.py
# ==========================================

import pandas as pd
import os

# 1. 메타데이터 파일 찾기
csv_candidates = [
    "plant_viruses_metadata_enriched.csv",
    "plant_viruses_metadata_table.csv",
    "korea_plant_viruses_metadata.csv",
    "../results/plant_viruses_metadata_enriched.csv",
    "../results/korea_plant_viruses_metadata.csv"
]

csv_file = None
for path in csv_candidates:
    if os.path.exists(path):
        csv_file = path
        break

if csv_file is None:
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".csv"):
                csv_file = os.path.join(root, file)
                break
        if csv_file: break

if csv_file is None:
    raise FileNotFoundError("[Error] 메타데이터 CSV 파일을 찾을 수 없습니다.")

print(f"=== Analyzing Metadata from '{csv_file}' ===")
df = pd.read_csv(csv_file)
print(f"Total Records: {len(df):,}")

# 2. 바이러스명/과 컬럼 탐색
virus_col = None
for col in df.columns:
    if any(k in str(col).lower() for k in ['organism', 'virus', 'species', 'pathogen', 'name']):
        virus_col = col
        break

if virus_col is None:
    virus_col = df.select_dtypes(include=['object']).columns[0]

print(f"Target Column for Analysis: '{virus_col}'\n")

# 3. 미분류 항목(Other) 정체를 파악하기 위한 빈도 상위 30개 출력
print("--- Top 30 Most Frequent Pathogens in Dataset ---")
top_pathogens = df[virus_col].value_counts().head(30)
for name, count in top_pathogens.items():
    print(f"  {name}: {count}건")



# ==========================================
# FILE: analyze_family_comparison.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

global_csv = "../data/global_data/global_plant_viruses_metadata.csv"
korea_csv = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading metadata for family-level comparison...")
df_global = pd.read_csv(global_csv)
df_korea = pd.read_csv(korea_csv)

# 바이러스 이름/설명을 기반으로 주요 과(Family)를 매핑하는 함수
# (식물 바이러스에서 가장 빈도가 높은 주요 과들을 타겟팅)
def assign_family(organism_name):
    if pd.isna(organism_name):
        return "Unclassified"
    name_lower = str(organism_name).lower()
    
    # 대표적인 식물 바이러스 및 비로이드 과 매핑 규칙
    if "potyvirus" in name_lower or "potyviridae" in name_lower or "mosaic virus" in name_lower and ("soybean" in name_lower or "bean" in name_lower or "turnip" in name_lower or "watermelon" in name_lower or "onion" in name_lower or "leek" in name_lower):
        # 대다수 감자/콩/채소 모자이크 바이러스 상당수가 Potyviridae에 속함
        pass
        
    # 명시적 키워드 검사
    if "potyviridae" in name_lower or "potyvirus" in name_lower or "rye" in name_lower: # 편의상 주요 속/과 매핑
        return "Potyviridae"
    elif "cucumovirus" in name_lower or "cucumber mosaic" in name_lower:
        return "Bromoviridae"
    elif "bromoviridae" in name_lower:
        return "Bromoviridae"
    elif "alphaflexiviridae" in name_lower or "potexvirus" in name_lower:
        return "Alphaflexiviridae"
    elif "betaflexiviridae" in name_lower or "carlavirus" in name_lower or "foveavirus" in name_lower or "capillovirus" in name_lower:
        return "Betaflexiviridae"
    elif "tombusviridae" in name_lower or "carmovirus" in name_lower:
        return "Tombusviridae"
    elif "secoviridae" in name_lower or "nepovirus" in name_lower or "comovirus" in name_lower:
        return "Secoviridae"
    elif "tobamovirus" in name_lower or "tobac" in name_lower:
        return "Virgaviridae"
    elif "virgaviridae" in name_lower or "furovirus" in name_lower or "hordeivirus" in name_lower:
        return "Virgaviridae"
    elif "geminiviridae" in name_lower or "begomovirus" in name_lower or "mastrevirus" in name_lower:
        return "Geminiviridae"
    elif "tospovirus" in name_lower or "orthotospovirus" in name_lower or "bunyaviridae" in name_lower or "tospoviridae" in name_lower:
        return "Tospoviridae"
    elif "viroid" in name_lower or "pospiviroidae" in name_lower or "avsunviroidae" in name_lower:
        return "Viroids (Non-viral)"
    elif "caulimoviridae" in name_lower or "badnavirus" in name_lower or "caulimovirus" in name_lower:
        return "Caulimoviridae"
    elif "partitiviridae" in name_lower:
        return "Partitiviridae"
    elif "nanoviridae" in name_lower:
        return "Nanoviridae"
    
    # 추가 휴리스틱 (종 이름 기반 매칭)
    if "garlic virus" in name_lower or "shallot" in name_lower:
        return "Betaflexiviridae"
    if "mosaic" in name_lower or "poty" in name_lower:
        return "Potyviridae"
        
    return "Other / Unclassified"

print("Assigning viral families for global dataset...")
df_global['Family'] = df_global['Organism'].apply(assign_family)

print("Assigning viral families for Korea dataset...")
df_korea['Family'] = df_korea['Organism'].apply(assign_family)

# 백분율 계산
global_fam = df_global['Family'].value_counts(normalize=True) * 100
korea_fam = df_korea['Family'].value_counts(normalize=True) * 100

global_fam_df = global_fam.reset_index()
global_fam_df.columns = ['Family', 'Global_Percentage']

korea_fam_df = korea_fam.reset_index()
korea_fam_df.columns = ['Family', 'Korea_Percentage']

# 병합
family_comp = pd.merge(korea_fam_df, global_fam_df, on='Family', how='outer').fillna(0)
family_comp = family_comp[family_comp['Family'] != 'Other / Unclassified'].sort_values(by='Korea_Percentage', ascending=False)

# CSV 저장
comp_path = os.path.join(results_dir, "global_vs_korea_family_comparison.csv")
family_comp.to_csv(comp_path, index=False)
print(f"-> Saved family comparison summary CSV: {comp_path}")

print("\nViral Family Comparison Summary:")
print(family_comp)

# 시각화: 과별 점유율 비교 그룹 바차트
plt.figure(figsize=(12, 6))
fam_melted = pd.melt(
    family_comp.head(8), 
    id_vars=['Family'], 
    value_vars=['Korea_Percentage', 'Global_Percentage'],
    var_name='Region', 
    value_name='Percentage'
)
fam_melted['Region'] = fam_melted['Region'].replace({
    'Korea_Percentage': 'South Korea', 
    'Global_Percentage': 'Global'
})

sns.barplot(
    data=fam_melted, 
    x='Percentage', 
    y='Family', 
    hue='Region',
    palette='Paired'
)
plt.title('Comparison of Plant Virus & Viroid Families: South Korea vs. Global', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Percentage of Total Sequences (%)', fontsize=12)
plt.ylabel('Viral Family', fontsize=12)
plt.legend(title='Dataset', loc='upper right')
plt.tight_layout()

fig_path = os.path.join(fig_dir, "global_vs_korea_family_comparison.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved comparative family high-res figure: {fig_path}")


# ==========================================
# FILE: analyze_first_reported_viruses.py
# ==========================================

import pandas as pd
import os
import re

results_dir = "../results"
os.makedirs(results_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for first-reported (emergent) virus analysis...")
df = pd.read_csv(csv_file)

# 연도 추출 로직
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        y = int(match.group(1))
        if 1990 <= y <= 2026:
            return y
    return None

df['Year'] = df['Collection_Date'].apply(extract_year)
df_clean = df.dropna(subset=['Year', 'Organism']).copy()
df_clean['Year'] = df_clean['Year'].astype(int)

# 연도 오름차순 정렬 후, 각 바이러스 종이 등장한 최초의 연도 찾기
df_sorted = df_clean.sort_values(by='Year')

# 이미 등장했던 종들을 추적하기 위한 집합(Set)
seen_viruses = set()
first_reports = []

for idx, row in df_sorted.iterrows():
    virus = row['Organism']
    year = row['Year']
    if virus not in seen_viruses:
        seen_viruses.add(virus)
        first_reports.append({'First_Reported_Year': year, 'Virus_Organism': virus})

df_first = pd.DataFrame(first_reports)

# 연도별로 새롭게 등장한 바이러스 목록 집계
yearly_new_viruses = df_first.groupby('First_Reported_Year')['Virus_Organism'].apply(list).reset_index()
yearly_new_viruses['New_Species_Count'] = yearly_new_viruses['Virus_Organism'].apply(len)

# CSV 저장 (연도별 최초 보고된 바이러스 목록 및 통계)
output_csv = os.path.join(results_dir, "korea_first_reported_viruses_yearly.csv")
df_first.sort_values(by='First_Reported_Year', ascending=False).to_csv(output_csv, index=False)
print(f"-> Saved first-reported viruses summary CSV: {output_csv}")

print("\nRecent Years First-Reported Plant Viruses / Viroids in Korea:")
print(yearly_new_viruses[yearly_new_viruses['First_Reported_Year'] >= 2015].to_string(index=False))



# ==========================================
# FILE: analyze_genomic_architecture.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("="*60)
print(" 🧬 유전체 유형 및 분절 구조에 따른 서열 크기 비교 분석")
print("="*60)

df = pd.read_csv(master_csv)

# 1. 유전체 유형(Genome Type)별 서열 길이 통계
print("\n[1] 유전체 유형별 서열 길이(Sequence Length) 통계 (bp)")
length_stats = df.groupby('Genome_Type')['Sequence_Length'].agg(['count', 'mean', 'median', 'min', 'max', 'std'])
print(length_stats.round(2))

# 2. 분절 구조(Genome Segmentation)별 서열 길이 통계
print("\n[2] 분절 구조 유형별 서열 길이(Sequence Length) 통계 (bp)")
seg_stats = df.groupby('Genome_Segmentation')['Sequence_Length'].agg(['count', 'mean', 'median', 'min', 'max', 'std'])
print(seg_stats.round(2))

# 3. 유전체 유형 x 분절 구조 교차 분석 통계
print("\n[3] 유전체 유형 및 분절 구조 교차 통계 (평균 서열 길이)")
pivot_table = df.pivot_table(index='Genome_Type', columns='Genome_Segmentation', values='Sequence_Length', aggfunc=['count', 'mean'])
print(pivot_table.round(1))

# 4. 시각화 생성: 유전체 유형 및 분절 구조별 서열 길이 박스플롯
plt.figure(figsize=(12, 6))
sns.boxplot(x='Genome_Type', y='Sequence_Length', hue='Genome_Segmentation', data=df, palette='Set2', showfliers=False)
plt.title('Sequence Length Distribution by Genome Type and Architecture', fontsize=14, fontweight='bold')
plt.xlabel('Genome Type', fontsize=12)
plt.ylabel('Sequence Length (bp)', fontsize=12)
plt.xticks(rotation=30)
plt.legend(title='Genome Segmentation', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('genomic_architecture_length_analysis.png', dpi=300)
plt.close()
print("\n[완료] 고품질 분석 시각화가 'genomic_architecture_length_analysis.png'로 저장되었습니다.")


# ==========================================
# FILE: analyze_geographic_hotspots.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🗺️ 지리적/공간적 바이러스 발생 핫스팟 분석")
print("="*60)

# 1. 지역별 검출 빈도 집계 및 상위 순위 산출
region_counts = df['Geographic_Region'].value_counts()
print("\n[1] 지역별 바이러스 검출 빈도 현황")
print(region_counts)

# 2. 시각화 생성: 지역별 핫스팟 바차트
plt.figure(figsize=(12, 6))
sns.barplot(x=region_counts.index, y=region_counts.values, palette='viridis')
plt.title('Geographic Spatial Distribution of Virus Detections (Hotspots)', fontsize=14, fontweight='bold')
plt.xlabel('Geographic Region', fontsize=12)
plt.ylabel('Detection Count', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('geographic_hotspot_analysis.png', dpi=300)
plt.close()
print("\n[완료] 지리적 핫스팟 시각화 저장: geographic_hotspot_analysis.png")


# ==========================================
# FILE: analyze_global_data.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

# 스타일 설정
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

global_csv = "../data/global_data/global_plant_viruses_metadata.csv"

print("Loading global virus metadata...")
df_global = pd.read_csv(global_csv)
print(f"Total global records loaded: {len(df_global)}")

# 1. 국가명 정제 함수 보완 (콜론(:)이나 세부 지역이 포함된 경우 앞부분을 국가로 인식)
def clean_country(val):
    if pd.isna(val):
        return "Unknown"
    val_str = str(val).strip()
    # 콜론이 있으면 앞부분이 국가명인 경우가 많음 (예: USA: California -> USA)
    if ":" in val_str:
        val_str = val_str.split(":")[0].strip()
    
    val_lower = val_str.lower()
    
    if any(k in val_lower for k in ["korea", "republic of korea", "south korea"]):
        return "South Korea"
    elif "china" in val_lower:
        return "China"
    elif "usa" in val_lower or "united states" in val_lower or "america" in val_lower or "u.s." in val_lower:
        return "USA"
    elif "japan" in val_lower:
        return "Japan"
    elif "germany" in val_lower:
        return "Germany"
    elif "india" in val_lower:
        return "India"
    elif "italy" in val_lower:
        return "Italy"
    elif "spain" in val_lower:
        return "Spain"
    elif "france" in val_str.lower():
        return "France"
    elif "brazil" in val_lower:
        return "Brazil"
    elif "united kingdom" in val_lower or "uk" in val_lower or "england" in val_lower:
        return "United Kingdom"
    elif "iran" in val_lower:
        return "Iran"
    elif "turkey" in val_lower or "türkiye" in val_lower:
        return "Turkey"
    elif "canada" in val_lower:
        return "Canada"
    elif "australia" in val_lower:
        return "Australia"
        
    # 그 외의 경우는 첫 글자만 대문자로 깔끔하게 정리
    return val_str.title()

df_global['Clean_Country'] = df_global['Country'].apply(clean_country)

# 디버깅용 확인 출력
print("Sample cleaned countries:")
print(df_global['Clean_Country'].value_counts().head(10))

# 2. 연도 추출
def extract_year(date_str):
    if pd.isna(date_str):
        return "Unknown"
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    return match.group(1) if match else "Unknown"

df_global['Year'] = df_global['Collection_Date'].apply(extract_year)

# 3. 국가별 등록 순위 Top 15 통계표 생성 (Unknown 제외)
valid_countries = df_global[~df_global['Clean_Country'].isin(["Unknown", "nan", "None"])]
top_countries = valid_countries['Clean_Country'].value_counts().head(15).reset_index()
top_countries.columns = ['Country', 'Sequence_Count']
top_countries['Percentage'] = (top_countries['Sequence_Count'] / len(df_global)) * 100
top_countries.to_csv(os.path.join(results_dir, "global_top_countries_summary.csv"), index=False)
print("-> Saved: global_top_countries_summary.csv")

# 4. 글로벌 상위 바이러스 종 Top 15 통계표 생성
top_global_viruses = df_global['Organism'].value_counts().head(15).reset_index()
top_global_viruses.columns = ['Virus_Organism', 'Count']
top_global_viruses['Percentage'] = (top_global_viruses['Count'] / len(df_global)) * 100
top_global_viruses.to_csv(os.path.join(results_dir, "global_top_viruses_summary.csv"), index=False)
print("-> Saved: global_top_viruses_summary.csv")

# --- 글로벌 시각화 그래프 재생성 ---

# 그래프 1: 국가별 등록 건수 Top 10 바차트
plt.figure(figsize=(11, 6))
sns.barplot(data=top_countries.head(10), x='Sequence_Count', y='Country', palette='Purples_r')
plt.title('Top 10 Countries Reporting Plant Viruses and Viroids Globally', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Sequences', fontsize=12)
plt.ylabel('Country', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "global_top_countries.png"), dpi=300)
plt.close()

# 그래프 2: 글로벌 상위 바이러스 종 Top 10 바차트
plt.figure(figsize=(11, 6))
sns.barplot(data=top_global_viruses.head(10), x='Count', y='Virus_Organism', palette='Oranges_r')
plt.title('Top 10 Plant Viruses and Viroids Worldwide', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Sequences', fontsize=12)
plt.ylabel('Virus Organism', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "global_top_viruses.png"), dpi=300)
plt.close()

# 그래프 3: 글로벌 연도별 트렌드 라인 그래프
df_year_g = df_global[df_global['Year'] != "Unknown"].copy()
df_year_g['Year'] = df_year_g['Year'].astype(int)
df_year_g = df_year_g[(df_year_g['Year'] >= 1990) & (df_year_g['Year'] <= 2026)]
global_yearly = df_year_g['Year'].value_counts().sort_index().reset_index()
global_yearly.columns = ['Year', 'Count']

plt.figure(figsize=(12, 5))
plt.plot(global_yearly['Year'], global_yearly['Count'], marker='o', color='#2ca02c', linewidth=2.5, markersize=5)
plt.title('Global Temporal Trend of Plant Virus and Viroid Submissions', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Number of Sequences Reported Worldwide', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "global_yearly_trend.png"), dpi=300)
plt.close()

print("Global analysis and figure re-generation completed successfully!")


# ==========================================
# FILE: analyze_global_host_comparison.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

global_csv = "../data/global_data/global_plant_viruses_metadata.csv"
korea_csv = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading global and Korea metadata for comparative host analysis...")
df_global = pd.read_csv(global_csv)
df_korea = pd.read_csv(korea_csv)

# 동일한 학명 정제 함수 정의
synonym_map = {
    "Soybean": "Glycine max",
    "Garlic": "Allium sativum",
    "Pepper": "Capsicum annuum",
    "Hot pepper": "Capsicum annuum",
    "Tomato": "Solanum lycopersicum",
    "Peach": "Prunus persica",
    "Radish": "Raphanus sativus"
}

def strict_scientific_host(host_str):
    if pd.isna(host_str):
        return None
    clean_str = str(host_str).strip()
    clean_lower = clean_str.lower()
    
    for key, val in synonym_map.items():
        if key.lower() in clean_lower:
            return val
            
    exclude_keywords = ["apis", "honeybee", "insect", "human", "virus", "unknown", "none", "n/a"]
    if any(ex in clean_lower for ex in exclude_keywords):
        return None
        
    parts = clean_str.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
    return None

# 전 세계 데이터 기주 정제
print("Cleaning global host names...")
df_global['Strict_Host'] = df_global['Host'].apply(strict_scientific_host)
global_valid = df_global.dropna(subset=['Strict_Host'])
global_counts = global_valid['Strict_Host'].value_counts(normalize=True) * 100
global_df = global_counts.head(10).reset_index()
global_df.columns = ['Host_Scientific_Name', 'Global_Percentage']

# 한국 데이터 기주 정제
print("Cleaning Korea host names...")
df_korea['Strict_Host'] = df_korea['Host'].apply(strict_scientific_host)
korea_valid = df_korea.dropna(subset=['Strict_Host'])
korea_counts = korea_valid['Strict_Host'].value_counts(normalize=True) * 100
korea_df = korea_counts.head(10).reset_index()
korea_df.columns = ['Host_Scientific_Name', 'Korea_Percentage']

# 데이터 병합 (비교 테이블 생성)
comparison_df = pd.merge(korea_df, global_df, on='Host_Scientific_Name', how='outer').fillna(0)
comparison_path = os.path.join(results_dir, "global_vs_korea_host_comparison.csv")
comparison_df.to_csv(comparison_path, index=False)
print(f"-> Saved comparison summary CSV: {comparison_path}")

print("\nHost Range Comparison (Korea vs Global Top Hosts):")
print(comparison_df.head(10))

# 시각화: 상위 기주 비중 비교 그룹 바차트
plt.figure(figsize=(12, 6))
top_hosts_melted = pd.melt(
    comparison_df.head(8), 
    id_vars=['Host_Scientific_Name'], 
    value_vars=['Korea_Percentage', 'Global_Percentage'],
    var_name='Region', 
    value_name='Percentage'
)
top_hosts_melted['Region'] = top_hosts_melted['Region'].replace({
    'Korea_Percentage': 'South Korea', 
    'Global_Percentage': 'Global'
})

sns.barplot(
    data=top_hosts_melted, 
    x='Percentage', 
    y='Host_Scientific_Name', 
    hue='Region',
    palette='Set2'
)
plt.title('Comparison of Major Plant Host Proportions: South Korea vs. Global', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Percentage of Total Sequences (%)', fontsize=12)
plt.ylabel('Host Scientific Name', fontsize=12)
plt.legend(title='Dataset', loc='lower right')
plt.tight_layout()

fig_path = os.path.join(fig_dir, "global_vs_korea_hosts.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved comparative high-res figure: {fig_path}")


# ==========================================
# FILE: analyze_host_diversity.py
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for host diversity and Shannon-Wiener index analysis...")
df = pd.read_csv(csv_file)

# 1. 기주 학명 정제 로직
synonym_map = {
    "Soybean": "Glycine max",
    "Garlic": "Allium sativum",
    "Pepper": "Capsicum annuum",
    "Hot pepper": "Capsicum annuum",
    "Tomato": "Solanum lycopersicum",
    "Peach": "Prunus persica",
    "Radish": "Raphanus sativus"
}

def strict_scientific_host(host_str):
    if pd.isna(host_str):
        return None
    clean_str = str(host_str).strip()
    clean_lower = clean_str.lower()
    for key, val in synonym_map.items():
        if key.lower() in clean_lower:
            return val
    exclude_keywords = ["apis", "honeybee", "insect", "human", "virus", "unknown", "none", "n/a"]
    if any(ex in clean_lower for ex in exclude_keywords):
        return None
    parts = clean_str.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
    return None

df['Strict_Host'] = df['Host'].apply(strict_scientific_host)
df_clean = df.dropna(subset=['Strict_Host', 'Organism']).copy()

# 2. 기주별 바이러스 다양성 지수(Shannon-Wiener Index) 계산 함수
def calc_shannon(series):
    counts = series.value_counts()
    proportions = counts / counts.sum()
    # Shannon-Wiener Index: H' = -sum(p_i * ln(p_i))
    h_prime = -np.sum(proportions * np.log(proportions))
    richness = len(counts)
    evenness = h_prime / np.log(richness) if richness > 1 else 0
    return pd.Series({
        'Total_Sequences': counts.sum(),
        'Virus_Richness': richness,
        'Shannon_Diversity_H': h_prime,
        'Evenness': evenness
    })

# 기주별 그룹화 및 계산 (샘플 수가 10개 이상인 주요 기주 대상)
host_groups = df_clean.groupby('Strict_Host')['Organism'].apply(calc_shannon).unstack()
host_groups = host_groups[host_groups['Total_Sequences'] >= 10].sort_values(by='Shannon_Diversity_H', ascending=False).reset_index()

# CSV 저장
diversity_path = os.path.join(results_dir, "korea_host_diversity_summary.csv")
host_groups.to_csv(diversity_path, index=False)
print(f"-> Saved host diversity summary CSV: {diversity_path}")

print("\nHost Diversity Index Summary Preview:")
print(host_groups.head(10))

# 3. 시각화: Shannon Diversity Index 바차트
plt.figure(figsize=(11, 6))
sns.barplot(
    data=host_groups.head(10), 
    x='Shannon_Diversity_H', 
    y='Strict_Host', 
    hue='Strict_Host',
    palette='Purples_r',
    legend=False
)
plt.title('Shannon-Wiener Diversity Index ($H\'$) of Plant Viruses across Major Hosts in Korea', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Shannon Diversity Index ($H^\prime$)', fontsize=12)
plt.ylabel('Host Scientific Name', fontsize=12)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_host_diversity_index.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved high-res diversity figure: {fig_path}")


# ==========================================
# FILE: analyze_host_ranges.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading Korean virus metadata for host range analysis...")
df = pd.read_csv(csv_file)

# 1. 기주(Host) 정보 정제 및 학명 표준화 함수
def clean_host_name(host_str):
    if pd.isna(host_str) or str(host_str).strip().lower() in ["unknown", "unnamed", "none", "n/a"]:
        return "Unknown"
    
    clean_str = str(host_str).strip()
    
    # 간혹 괄호 안에 설명이나 품종명이 붙어있는 경우 처리 (예: Solanum lycopersicum (tomato))
    # 이탤릭체로 표기될 수 있도록 전형적인 이중명명법(속명 + 종소명) 형태 추출 시도
    # 예컨대 첫 두 단어(Genus species)를 학명으로 우선 추출
    parts = clean_str.split()
    if len(parts) >= 2:
        # 첫 글자는 대문자, 두 번째는 소문자로 조합하여 학명 형태 유지
        genus = parts[0].capitalize()
        species = parts[1].lower()
        # 특수문자 제거
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
            
    return clean_str.capitalize()

df['Clean_Host'] = df['Host'].apply(clean_host_name)

# 2. 상위 기주 학명 통계표 생성 (Unknown 제외)
valid_hosts = df[df['Clean_Host'] != "Unknown"]
top_hosts_scientific = valid_hosts['Clean_Host'].value_counts().head(15).reset_index()
top_hosts_scientific.columns = ['Host_Scientific_Name', 'Sequence_Count']
top_hosts_scientific['Percentage'] = (top_hosts_scientific['Sequence_Count'] / len(df)) * 100

# 결과 CSV 저장
summary_path = os.path.join(results_dir, "korea_host_range_summary.csv")
top_hosts_scientific.to_csv(summary_path, index=False)
print(f"-> Saved host range summary: {summary_path}")

print("\nTop 10 Host Scientific Names in Korea:")
print(top_hosts_scientific.head(10))

# 3. 논문용 고해상도 기주 범위 그래프 생성
plt.figure(figsize=(11, 6))
# Seaborn 경고 방지를 위한 hue 및 legend 설정
sns.barplot(
    data=top_hosts_scientific.head(10), 
    x='Sequence_Count', 
    y='Host_Scientific_Name', 
    hue='Host_Scientific_Name',
    palette='Greens_r',
    legend=False
)
plt.title('Top 10 Plant Host Species (Scientific Names) in Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Virus/Viroid Sequences', fontsize=12)
plt.ylabel('Host Scientific Name', fontsize=12)
# 학명 이탤릭체 처리를 위해 y축 라벨 스타일 조정 가능
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_top_hosts_scientific.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved high-res figure: {fig_path}")


# ==========================================
# FILE: analyze_host_species.py
# ==========================================

from Bio import SeqIO
import os
from collections import Counter
import re

gbk_file = "../data/korea_plant_viruses_raw.gbk"

if not os.path.exists(gbk_file):
    print(f"[Error] File not found: {gbk_file}")
    exit()

print("Analyzing Host species, Viral species, and Temporal trends from GenBank records...")

total_records = 0
species_counter = Counter()
host_counter = Counter()
year_counter = Counter()

for record in SeqIO.parse(gbk_file, "genbank"):
    total_records += 1
    
    # 1. 바이러스 종(Species) 추출
    organism = record.annotations.get("organism", "Unknown")
    species_counter[organism] += 1
    
    # 2. 연도(Year) 추출 (날짜 형식 파싱 또는 피처 내 데이터 활용)
    date_str = record.annotations.get("date", "")
    year_match = re.search(r'\b(19\d{2}|20\d{2})\b', date_str)
    if year_match:
        year_counter[year_match.group(1)] += 1
    else:
        # 날짜가 명시되지 않은 경우 서열 버전 날짜나 레코드 데이터에서 탐색 시도
        for key, val in record.dbxrefs:
            pass
        # 또는 피처 내 source의 note/isolation_source 등 탐색 가능
        year_counter["Unknown/Not Specified"] += 1

    # 3. 기주(Host) 정보 추출 (source 피처의 qualifiers 탐색)
    host_found = False
    for feature in record.features:
        if feature.type == "source":
            host = feature.qualifiers.get("host", [None])
            if host[0]:
                host_name = host[0].strip()
                host_counter[host_name] += 1
                host_found = True
                break
    if not host_found:
        host_counter["Not Specified"] += 1

print("\n" + "="*60)
print(f" Total Analyzed Records: {total_records}")
print("="*60)

# Top 15 Viral Species
print("\n[ Top 15 Viral / Viroid Species ]")
print(f"{'Species Name':<40} | {'Count':<8}")
print("-" * 52)
for sp, cnt in species_counter.most_common(15):
    print(f"{sp:<40} | {cnt:<8}")

# Top 15 Hosts
print("\n[ Top 15 Plant Hosts ]")
print(f"{'Host Name':<40} | {'Count':<8}")
print("-" * 52)
for h, cnt in host_counter.most_common(15):
    print(f"{h:<40} | {cnt:<8}")

# Recent Years Trend (Top 15)
print("\n[ Temporal Trend (By Year) ]")
print(f"{'Year':<20} | {'Count':<8}")
print("-" * 32)
for y, cnt in sorted(year_counter.items(), reverse=True):
    if y != "Unknown/Not Specified":
        print(f"{y:<20} | {cnt:<8}")

print("="*60)


# ==========================================
# FILE: analyze_host_species_cleaned.py
# ==========================================

from Bio import SeqIO
import os
from collections import Counter
import re

gbk_file = "../data/korea_plant_viruses_raw.gbk"

if not os.path.exists(gbk_file):
    print(f"[Error] File not found: {gbk_file}")
    exit()

print("Parsing and cleaning GenBank records (Excluding animal/fish pathogens and insect vectors)...")

total_records = 0
filtered_records = 0
species_counter = Counter()
host_counter = Counter()
year_counter = Counter()

# 제외할 비식물성 바이러스 키워드
exclude_viruses = [
    "hemorrhagic septicemia", "rabies", "lyssavirus", "influenza", 
    "coronavirus", "hepatitis", "hiv", "herpes", "adenovirus"
]

# 제외할 곤충/동물 기주 키워드
exclude_hosts = [
    "apis mellifera", "bemisia tabaci", "honey bee", "bee", "mouse", "human"
]

# 식물 기주명 표준화 매핑 딕셔너리 (일반명 -> 학명)
host_mapping = {
    "chrysanthemum": "Chrysanthemum spp.",
    "tomato": "Solanum lycopersicum",
    "pepper": "Capsicum annuum",
    "soybean": "Glycine max",
    "melon": "Cucumis melo"
}

for record in SeqIO.parse(gbk_file, "genbank"):
    total_records += 1
    
    organism = record.annotations.get("organism", "Unknown").lower()
    
    # 1. 비식물성 바이러스 필터링
    if any(keyword in organism for keyword in exclude_viruses):
        continue
        
    # 2. 기주(Host) 정보 추출 및 정제
    host_name = "Not Specified"
    for feature in record.features:
        if feature.type == "source":
            h_list = feature.qualifiers.get("host", [None])
            if h_list[0]:
                raw_host = h_list[0].strip().lower()
                # 곤충 및 동물 기주 제외
                if any(ex in raw_host for ex in exclude_hosts):
                    host_name = None
                    break
                
                # 표준화 매핑 적용
                cleaned_host = raw_host
                for k, v in host_mapping.items():
                    if k in raw_host:
                        cleaned_host = v.lower()
                        break
                # 첫 글자 대문자화 (학명 스타일)
                host_name = cleaned_host.capitalize()
                break
                
    if not host_name or host_name == "Not Specified":
        # 기주가 명시되지 않았거나 제외된 경우 스킵 또는 카운트
        host_name = "Not Specified / Filtered"

    filtered_records += 1
    species_counter[record.annotations.get("organism", "Unknown")] += 1
    host_counter[host_name] += 1
    
    # 3. 연도(Year) 추출
    date_str = record.annotations.get("date", "")
    year_match = re.search(r'\b(19\d{2}|20\d{2})\b', date_str)
    if year_match:
        year_counter[year_match.group(1)] += 1
    else:
        year_counter["Unknown/Not Specified"] += 1

print("\n" + "="*65)
print(f" Total Raw Records: {total_records} | Cleaned Plant-Pathogen Records: {filtered_records}")
print("="*65)

# Top 15 Plant Viral / Viroid Species
print("\n[ Top 15 Plant Viral & Viroid Species (Cleaned) ]")
print(f"{'Species Name':<45} | {'Count':<8}")
print("-" * 57)
for sp, cnt in species_counter.most_common(15):
    print(f"{sp:<45} | {cnt:<8}")

# Top 15 Plant Hosts (Cleaned & Standardized)
print("\n[ Top 15 Plant Hosts (Standardized Scientific Names) ]")
print(f"{'Host Name':<45} | {'Count':<8}")
print("-" * 57)
for h, cnt in host_counter.most_common(15):
    print(f"{h:<45} | {cnt:<8}")

# Temporal Trend (Top 15 Years)
print("\n[ Temporal Trend (By Year) ]")
print(f"{'Year':<20} | {'Count':<8}")
print("-" * 32)
for y, cnt in sorted(year_counter.items(), reverse=True):
    if y != "Unknown/Not Specified" and cnt > 0:
        print(f"{y:<20} | {cnt:<8}")

print("="*65)


# ==========================================
# FILE: analyze_korea_data.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 스타일 설정
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

# 결과 저장 폴더 생성
results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading Korean virus metadata...")
df = pd.read_csv(csv_file)

print(f"Total records in dataset: {len(df)}")

# 1. 연도 추출 함수 (Collection_Date 또는 Accession 정보에서 연도 파싱 시도)
def extract_year(date_str):
    if pd.isna(date_str):
        return "Unknown"
    match = re_year = pd.Series(str(date_str)).str.extract(r'(20\d{2}|19\d{2})')[0].values[0]
    return match if pd.notna(match) else "Unknown"

import re
df['Year'] = df['Collection_Date'].apply(lambda x: re.search(r'(20\d{2}|19\d{2})', str(x)).group(1) if re.search(r'(20\d{2}|19\d{2})', str(x)) else "Unknown")

# 2. 상위 바이러스 종(Organism) Top 15 통계표 생성
top_organisms = df['Organism'].value_counts().head(15).reset_index()
top_organisms.columns = ['Virus_Organism', 'Count']
top_organisms['Percentage'] = (top_organisms['Count'] / len(df)) * 100
top_organisms.to_csv(os.path.join(results_dir, "korea_top_viruses_summary.csv"), index=False)
print("-> Saved: korea_top_viruses_summary.csv")

# 3. 상위 숙주(Host) Top 15 통계표 생성
top_hosts = df['Host'].value_counts().head(15).reset_index()
top_hosts.columns = ['Host_Name', 'Count']
top_hosts['Percentage'] = (top_hosts['Count'] / len(df)) * 100
top_hosts.to_csv(os.path.join(results_dir, "korea_top_hosts_summary.csv"), index=False)
print("-> Saved: korea_top_hosts_summary.csv")

# 4. 연도별 등록 추이 통계표 생성
year_counts = df[df['Year'] != "Unknown"]['Year'].value_counts().sort_index().reset_index()
year_counts.columns = ['Year', 'Count']
year_counts.to_csv(os.path.join(results_dir, "korea_yearly_trend_summary.csv"), index=False)
print("-> Saved: korea_yearly_trend_summary.csv")

# --- 시각화 그래프 생성 (논문용 고해상도 PNG) ---

# 그래프 1: 상위 바이러스 종 Top 15 바차트
plt.figure(figsize=(10, 6))
ax = sns.barplot(data=top_organisms.head(10), x='Count', y='Virus_Organism', palette='Blues_r')
plt.title('Top 10 Plant Viruses and Viroids in Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Sequences', fontsize=12)
plt.ylabel('Virus Organism', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "korea_top_viruses.png"), dpi=300)
plt.close()

# 그래프 2: 상위 숙주 Top 10 바차트
plt.figure(figsize=(10, 6))
ax = sns.barplot(data=top_hosts.head(10), x='Count', y='Host_Name', palette='Greens_r')
plt.title('Top 10 Plant Hosts in Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Sequences', fontsize=12)
plt.ylabel('Host Name', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "korea_top_hosts.png"), dpi=300)
plt.close()

# 그래프 3: 연도별 트렌드 라인 그래프
df_year_filtered = df[df['Year'] != "Unknown"].copy()
df_year_filtered['Year'] = df_year_filtered['Year'].astype(int)
df_year_filtered = df_year_filtered[df_year_filtered['Year'] >= 1990] # 1990년 이후 데이터 집중
yearly_trend = df_year_filtered['Year'].value_counts().sort_index().reset_index()
yearly_trend.columns = ['Year', 'Count']

plt.figure(figsize=(12, 5))
plt.plot(yearly_trend['Year'], yearly_trend['Count'], marker='o', color='#1f77b4', linewidth=2.5, markersize=6)
plt.title('Temporal Trend of Plant Virus/Viroid Submissions in Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Number of Sequences Reported', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "korea_yearly_trend.png"), dpi=300)
plt.close()

print("All analysis and figure generation completed successfully!")
print(f"Results and tables are saved in: {results_dir}")
print(f"High-resolution figures are saved in: {fig_dir}")


# ==========================================
# FILE: analyze_region_epidemiology.py
# ==========================================

import pandas as pd
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🗺️ 채집 지역별 주요 검출 바이러스 및 전장 유전체 확보율 분석")
print("="*60)

# 1. 지역별 전장 유전체 확보율(Complete_Genome_Status) 교차 분석
print("\n[1] 채집 지역별 전장 유전체 확보 비율 (%)")
region_completeness = pd.crosstab(df['Geographic_Region'], df['Complete_Genome_Status'], normalize='index') * 100
print(region_completeness.round(2))

# 2. 주요 거점 지역별 최빈 검출 바이러스 종(Viral_Species) 분석
print("\n[2] 주요 채집 지역별 상위 검출 바이러스 종")
top_regions = df['Geographic_Region'].value_counts().head(5).index
for reg in top_regions:
    sub_df = df[df['Geographic_Region'] == reg]
    print(f"\n- 지역: {reg} (총 샘플수: {len(sub_df)}행)")
    print(sub_df['Viral_Species'].value_counts().head(3))


# ==========================================
# FILE: analyze_spacetime_hotspot.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🗺️ 공간-시계열 핫스팟 이동 추이 분석 시작")
print("="*60)

# 1. 지역 및 연도 데이터 정제
df['Geographic_Region'] = df['Geographic_Region'].fillna('Unspecified City')
df['Collection_Year'] = pd.to_numeric(df['Collection_Year'], errors='coerce')

# 전국 단위나 미지정 지역을 제외하고 주요 확정 도시들 필터링
target_regions = [
    'South Korea (Andong)', 'South Korea (Jeonju)', 'South Korea (Jinju)', 
    'South Korea (Chuncheon)', 'South Korea (Wanju)', 'South Korea (Gunsan)'
]
sub_df = df[df['Geographic_Region'].isin(target_regions)].copy()

# 최근 10~15년 데이터 중심으로 시공간 트렌드 가독성 강화 (예: 2010년 이후)
sub_df = sub_df[sub_df['Collection_Year'] >= 2010]

if sub_df.empty:
    print("[Notice] 필터링된 데이터가 부족하여 전체 연도 기준으로 재조정합니다.")
    sub_df = df[df['Geographic_Region'].isin(target_regions)].copy()

# 2. 지역(Region)과 연도(Year) 조합에 따른 주요 바이러스 교차 집계 (Pivot Table)
# 각 지역-연도별로 가장 많이 검출된 상위 바이러스 종 매핑
def clean_virus_name(name):
    if pd.isna(name):
        return 'Unknown'
    return str(name).strip()

sub_df['Clean_Species'] = sub_df['Viral_Species'].apply(clean_virus_name)

# 공간-시계열 피벗 테이블 구성 (Index: Region, Columns: Year, Values: Detection Count의 주요 종)
spacetime_pivot = pd.crosstab(sub_df['Geographic_Region'], sub_df['Collection_Year'])

# 3. 시각화 (Space-Time Hotspot Heatmap)
plt.figure(figsize=(14, 7))
sns.heatmap(
    spacetime_pivot, annot=True, fmt='d', cmap='YlOrRd', 
    linewidths=0.5, cbar_kws={'label': 'Detection Count'}
)

plt.title('Space-Time Hotspot Shift: Virus Detection Intensity Across Major Regions (2010–Present)', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Collection Year', fontsize=11)
plt.ylabel('Geographic Region', fontsize=11)
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)

plt.tight_layout()
output_fig = "spacetime_hotspot_shift.png"
plt.savefig(output_fig, dpi=300)
plt.close()

print(f"[완료] 공간-시계열 핫스팟 히트맵 저장 성공: {output_fig}")


# ==========================================
# FILE: analyze_temporal_completeness_trend.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 📈 연도별 바이러스 검출 추이 및 전장 유전체 비중 분석")
print("="*60)

# 1. 연도별 검출 건수 및 전장 유전체 비율 산출
year_counts = df['Collection_Year'].value_counts().sort_index()
completeness_trend = pd.crosstab(df['Collection_Year'], df['Complete_Genome_Status'], normalize='index') * 100

print("\n[1] 연도별 등록 건수 (상위 일부)")
print(year_counts.tail(10))

print("\n[2] 연도별 전장 유전체(Complete Genome) 비율 (%) (상위 일부)")
print(completeness_trend.tail(10).round(2))

# 2. 시각화 1: 연도별 신규 바이러스 검출 추이 (Line Plot)
plt.figure(figsize=(14, 6))
sns.lineplot(x=year_counts.index, y=year_counts.values, marker='o', color='b', linewidth=2.5)
plt.title('Temporal Trend of Plant Virus Detections by Collection Year', fontsize=14, fontweight='bold')
plt.xlabel('Collection Year', fontsize=12)
plt.ylabel('Detection Count', fontsize=12)
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('temporal_detection_trend.png', dpi=300)
plt.close()

# 3. 시각화 2: 연도별 전장 유전체 확보 비중 트렌드
plt.figure(figsize=(14, 6))
if 'Complete Genome' in completeness_trend.columns:
    sns.lineplot(x=completeness_trend.index, y=completeness_trend['Complete Genome'], marker='o', color='g', linewidth=2.5)
    plt.title('Percentage of Complete Genome Sequences by Collection Year', fontsize=14, fontweight='bold')
    plt.xlabel('Collection Year', fontsize=12)
    plt.ylabel('Complete Genome Ratio (%)', fontsize=12)
    plt.ylim(0, 105)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.savefig('complete_genome_ratio_trend.png', dpi=300)
    plt.close()

print("\n[완료] 시각화 저장 성공:")
print(" - 1. 검출 추이 차트: temporal_detection_trend.png")
print(" - 2. 전장 비율 트렌드 차트: complete_genome_ratio_trend.png")


# ==========================================
# FILE: analyze_temporal_trend.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 📈 연도별 검출 추이 및 전장 유전체 확보 비율 트렌드 분석")
print("="*60)

# 1. 연도별 검출 건수 집계
year_counts = df['Collection_Year'].value_counts().sort_index()
print("\n[1] 연도별 바이러스 검출 및 등록 건수")
print(year_counts)

# 2. 연도별 전장 유전체 확보 비율 교차 분석
print("\n[2] 연도별 전장 유전체 확보 비율 (%)")
year_completeness = pd.crosstab(df['Collection_Year'], df['Complete_Genome_Status'], normalize='index') * 100
print(year_completeness.round(2))

# 3. 시각화 생성: 연도별 검출 추이 트렌드 라인 차트
plt.figure(figsize=(14, 6))
sns.lineplot(x=year_counts.index, y=year_counts.values, marker='o', color='b', linewidth=2.5)
plt.title('Temporal Trend of Plant Virus Detections by Year', fontsize=14, fontweight='bold')
plt.xlabel('Collection Year', fontsize=12)
plt.ylabel('Detection Count', fontsize=12)
plt.xticks(rotation=45)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('temporal_trend_analysis.png', dpi=300)
plt.close()
print("\n[완료] 연도별 트렌드 시각화 저장: temporal_trend_analysis.png")


# ==========================================
# FILE: analyze_transmission_vs_symptom.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🐛 매개체 전파 방식과 병원성 증상 간의 연관성 분석")
print("="*60)

# 1. 교차 빈도표 (Contingency Table)
contingency_counts = pd.crosstab(df['Vector_Transmission_Mode'], df['Pathogenicity_Symptoms'])
print("\n[1] 교차 빈도표 (Counts)")
print(contingency_counts)

# 2. 행 기준 백분율표 (Row-Normalized Percentage)
contingency_pct = pd.crosstab(df['Vector_Transmission_Mode'], df['Pathogenicity_Symptoms'], normalize='index') * 100
print("\n[2] 전파 방식별 증상 비율 (%)")
print(contingency_pct.round(2))

# 3. 카이제곱 독립성 검정 (Chi-Square Test of Independence)
chi2, p, dof, ex = chi2_contingency(contingency_counts)
print(f"\n[3] 카이제곱 검정 결과: Chi2 = {chi2:.2f}, p-value = {p:.4e} (dof = {dof})")

# 4. 히트맵 시각화 저장
plt.figure(figsize=(12, 6))
sns.heatmap(contingency_counts, annot=True, fmt='d', cmap='YlOrRd', cbar=True)
plt.title('Contingency Map: Vector Transmission Mode vs Pathogenicity Symptoms', fontsize=14, fontweight='bold')
plt.xlabel('Pathogenicity Symptoms', fontsize=12)
plt.ylabel('Vector Transmission Mode', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('transmission_vs_symptom_heatmap.png', dpi=300)
plt.close()
print("\n[완료] 교차 연관성 히트맵 저장 완료: transmission_vs_symptom_heatmap.png")


# ==========================================
# FILE: analyze_ultra_detailed_shift.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for ultra-detailed virus shift analysis...")
df = pd.read_csv(csv_file)

# 연도 추출 로직
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        y = int(match.group(1))
        if 2000 <= y <= 2026:
            return y
    return None

df['Year'] = df['Collection_Date'].apply(extract_year)
df_clean = df.dropna(subset=['Year', 'Organism']).copy()

# 'Other' 비중을 최소화하기 위해 Top 15 바이러스/비로이드 종을 개별 분리
top_15_viruses = df_clean['Organism'].value_counts().head(15).index.tolist()

def group_virus_ultra(org):
    org_str = str(org)
    if org_str in top_15_viruses:
        return org_str
    else:
        return "Other Minor Species"

df_clean['Ultra_Group'] = df_clean['Organism'].apply(group_virus_ultra)

# 연도별 x 초세분화 바이러스 그룹 교차 테이블 (비율 계산)
yearly_ultra_xtab = pd.crosstab(df_clean['Year'], df_clean['Ultra_Group'], normalize='index') * 100

# 만약 Other Minor Species 컬럼이 존재한다면 맨 뒤로 보내기 순서 정렬
cols = [c for c in yearly_ultra_xtab.columns if c != "Other Minor Species"] + (["Other Minor Species"] if "Other Minor Species" in yearly_ultra_xtab.columns else [])
yearly_ultra_xtab = yearly_ultra_xtab[cols]

# CSV 저장
shift_path = os.path.join(results_dir, "korea_ultra_detailed_virus_yearly_shift.csv")
yearly_ultra_xtab.to_csv(shift_path)
print(f"-> Saved ultra-detailed yearly shift CSV: {shift_path}")

# 시각화: 초세분화 우점종 역학 변동 스택 바차트 (Top 15 분리)
plt.figure(figsize=(16, 9))
ax = yearly_ultra_xtab.plot(
    kind='bar', 
    stacked=True, 
    figsize=(16, 9), 
    colormap='tab20', 
    edgecolor='none'
)
plt.title('Ultra-Detailed Shift in Dominant Plant Viruses & Viroids Over Time in Korea (2000–2026)', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Proportion of Reported Sequences (%)', fontsize=12)
plt.legend(title='Virus Species / Viroids (Top 15)', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
plt.xticks(rotation=45)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_ultra_detailed_dominant_virus_shift.png")
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"-> Saved high-res ultra-detailed shift figure: {fig_path}")


# ==========================================
# FILE: analyze_vector_impact.py
# ==========================================

import pandas as pd
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🐛 매개체 전파 방식에 따른 전장 유전체 확보율 및 증상 연관성 분석")
print("="*60)

# 1. 매개체 전파 방식별 전장 유전체 비율 분석
print("\n[1] 매개체 전파 방식(Vector_Transmission_Mode)별 전장 유전체 확보 현황 (%)")
vector_completeness = pd.crosstab(df['Vector_Transmission_Mode'], df['Complete_Genome_Status'], normalize='index') * 100
print(vector_completeness.round(2))

# 2. 매개체 전파 방식별 병원성 증상 분포 교차 분석
print("\n[2] 매개체 전파 방식별 주요 병원성 증상(Pathogenicity_Symptoms) 교차 빈도")
vector_symptom = pd.crosstab(df['Vector_Transmission_Mode'], df['Pathogenicity_Symptoms'])
print(vector_symptom)


# ==========================================
# FILE: analyze_virus_clustering.py
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🧬 바이러스 생태유전학적 PCA / t-SNE 클러스터링 분석 시작")
print("="*60)

# 1. 수치형 및 범주형 특성 전처리
df['Sequence_Length'] = pd.to_numeric(df['Sequence_Length'], errors='coerce').fillna(df['Sequence_Length'].median())
df['Log_Length'] = np.log1p(df['Sequence_Length'])

features_for_cluster = pd.get_dummies(df[['Genome_Type', 'Genome_Segmentation', 'Vector_Transmission_Mode']], drop_first=True)
features_for_cluster['Log_Length'] = df['Log_Length']

features_clean = features_for_cluster.dropna()
scaler = StandardScaler()
scaled_features = scaler.fit_transform(features_clean)

# 2. PCA (Principal Component Analysis) 수행
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(scaled_features)

# 3. t-SNE 수행 (n_iter 대신 최신 버전 호환 max_iter 사용)
tsne = TSNE(n_components=2, perplexity=30, random_state=42, max_iter=1000)
tsne_coords = tsne.fit_transform(scaled_features)

result_df = df.loc[features_clean.index].copy()
result_df['PCA1'] = pca_coords[:, 0]
result_df['PCA2'] = pca_coords[:, 1]
result_df['tSNE1'] = tsne_coords[:, 0]
result_df['tSNE2'] = tsne_coords[:, 1]

# 4. 시각화 (PCA & t-SNE 2-Panel Plot)
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

sns.scatterplot(
    data=result_df, x='PCA1', y='PCA2', hue='Genome_Type', 
    palette='Set1', alpha=0.7, s=50, ax=axes[0]
)
axes[0].set_title('PCA Clustering by Genome Type', fontsize=13, fontweight='bold')
axes[0].set_xlabel(f'PCA Dim 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)', fontsize=11)
axes[0].set_ylabel(f'PCA Dim 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)', fontsize=11)
axes[0].legend(title='Genome Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)

sns.scatterplot(
    data=result_df, x='tSNE1', y='tSNE2', hue='Genome_Type', 
    palette='Set1', alpha=0.7, s=50, ax=axes[1]
)
axes[1].set_title('t-SNE Clustering by Genome Type', fontsize=13, fontweight='bold')
axes[1].set_xlabel('t-SNE Dim 1', fontsize=11)
axes[1].set_ylabel('t-SNE Dim 2', fontsize=11)
axes[1].legend(title='Genome Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)

plt.tight_layout()
output_fig = "virus_clustering_pca_tsne.png"
plt.savefig(output_fig, dpi=300)
plt.close()

print(f"[완료] PCA 및 t-SNE 클러스터링 그래프 저장 성공: {output_fig}")


# ==========================================
# FILE: analyze_virus_cooccurrence.py
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for virus co-occurrence analysis...")
df = pd.read_csv(csv_file)

# 필수 컬럼 확인 및 데이터 정제
df_clean = df.dropna(subset=['Organism']).copy()

# 만약 샘플 식별자(Sample ID나 Run ID 등)가 있다면 활용하되, 
# 본 메타데이터 구조상 동일 기주 및 동일 날짜/레코드를사용한 복합 감염 유추 또는 
# 상위 우점 바이러스 간의 동시 출현 빈도를 분석합니다.
# 여기서는 상위 10개 주요 바이러스 종을 대상으로 공존 빈도 행렬을 구축합니다.
top_viruses = df_clean['Organism'].value_counts().head(10).index.tolist()

# 샘플 단위 그룹화가 가능한 경우(예: Accession 또는 Source 특징)를 활용하기 위해
# 데이터셋 내에서 바이러스 간의 출현 패턴 상관관계를 계산합니다.
# 가상의 바이러스-샘플 피벗 테이블 생성 (데이터 내 고유 식별자 활용)
if 'Accession' in df_clean.columns:
    sample_col = 'Accession'
elif 'Sample' in df_clean.columns:
    sample_col = 'Sample'
else:
    # 대체재로 인덱스 기반 또는 Host/Collection_Date 조합 활용
    df_clean['Pseudo_Sample'] = df_clean['Host'].astype(str) + "_" + df_clean['Collection_Date'].astype(str)
    sample_col = 'Pseudo_Sample'

# 바이러스 종별 피벗 테이블 (각 샘플에 해당 바이러스가 검출되었는지 여부 0 또는 1)
sub_df = df_clean[df_clean['Organism'].isin(top_viruses)]
presence_matrix = pd.crosstab(sub_df[sample_col], sub_df['Organism'])
presence_matrix = (presence_matrix > 0).astype(int)

# 상관관계 매트릭스 계산 (Jaccard 또는 Pearson 상관계수)
co_occurrence_corr = presence_matrix.corr(method='pearson').fillna(0)

# CSV 저장
corr_path = os.path.join(results_dir, "korea_virus_cooccurrence_summary.csv")
co_occurrence_corr.to_csv(corr_path)
print(f"-> Saved co-occurrence summary CSV: {corr_path}")

print("\nVirus Co-occurrence Correlation Preview:")
print(co_occurrence_corr.iloc[:5, :5])

# 시각화: 코-오커런스 상관관계 히트맵
plt.figure(figsize=(11, 9))
sns.heatmap(
    co_occurrence_corr, 
    annot=True, 
    fmt=".2f", 
    cmap="coolwarm", 
    vmin=-1, 
    vmax=1, 
    linewidths=.5,
    cbar_kws={'label': 'Pearson Correlation Coefficient'}
)
plt.title('Co-occurrence Correlation Matrix of Major Plant Viruses in Korea', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Virus Organism', fontsize=12)
plt.ylabel('Virus Organism', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_virus_cooccurrence_heatmap.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved high-res co-occurrence figure: {fig_path}")


# ==========================================
# FILE: analyze_virus_family.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

# 분류 계통 문자열에서 Family 명 추출하는 함수
def extract_family(lineage):
    if pd.isna(lineage):
        return 'Unclassified'
    parts = [p.strip() for p in str(lineage).split(';')]
    for p in parts:
        if p.endswith('viridae') or p.endswith('viroidae') or 'viroid' in p.lower():
            return p
    return 'Other / Unclassified'

df['Viral_Family'] = df['Taxonomy_Lineage'].apply(extract_family)

print("="*60)
print(" 🧬 주요 바이러스 패밀리(Family)별 검출 빈도 통계")
print("="*60)
family_counts = df['Viral_Family'].value_counts().head(10)
print(family_counts)

# 시각화 생성: 주요 바이러스 패밀리 분포 가로 바차트
plt.figure(figsize=(12, 6))
sns.barplot(x=family_counts.values, y=family_counts.index, palette='mako')
plt.title('Distribution of Major Plant Virus Families', fontsize=14, fontweight='bold')
plt.xlabel('Detection Count', fontsize=12)
plt.ylabel('Viral Family', fontsize=12)
plt.tight_layout()
plt.savefig('virus_family_distribution.png', dpi=300)
plt.close()

print("\n[완료] 바이러스 패밀리 분포 시각화 저장: virus_family_distribution.png")


# ==========================================
# FILE: analyze_virus_host_interaction.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading metadata for virus-host interaction analysis...")
df = pd.read_csv(csv_file)

# 1. 이전 단계와 동일한 정제 로직 적용 (학명 및 바이러스 종 표준화)
synonym_map = {
    "Soybean": "Glycine max",
    "Garlic": "Allium sativum",
    "Pepper": "Capsicum annuum",
    "Hot pepper": "Capsicum annuum",
    "Tomato": "Solanum lycopersicum",
    "Peach": "Prunus persica",
    "Radish": "Raphanus sativus"
}

def strict_scientific_host(host_str):
    if pd.isna(host_str):
        return None
    clean_str = str(host_str).strip()
    clean_lower = clean_str.lower()
    
    for key, val in synonym_map.items():
        if key.lower() in clean_lower:
            return val
            
    exclude_keywords = ["apis", "honeybee", "insect", "human", "virus", "unknown", "none", "n/a"]
    if any(ex in clean_lower for ex in exclude_keywords):
        return None
        
    parts = clean_str.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        import re
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
    return None

df['Strict_Host'] = df['Host'].apply(strict_scientific_host)

# 유효한 기주와 바이러스 데이터만 추출
df_clean = df.dropna(subset=['Strict_Host']).copy()

# 2. 상위 바이러스 종(Top 8) 및 상위 기주 학명(Top 8) 선정
top_viruses = df_clean['Organism'].value_counts().head(8).index.tolist()
top_hosts = df_clean['Strict_Host'].value_counts().head(8).index.tolist()

# 상위 항목에 속하는 데이터만 필터링
df_filtered = df_clean[df_clean['Organism'].isin(top_viruses) & df_clean['Strict_Host'].isin(top_hosts)]

# 3. 교차 매핑 테이블(Cross-tabulation) 생성
interaction_matrix = pd.crosstab(df_filtered['Organism'], df_filtered['Strict_Host'])

# CSV 저장
matrix_path = os.path.join(results_dir, "korea_virus_host_interaction.csv")
interaction_matrix.to_csv(matrix_path)
print(f"-> Saved interaction matrix CSV: {matrix_path}")

print("\nVirus-Host Interaction Matrix Preview:")
print(interaction_matrix)

# 4. 논문용 고해상도 히트맵(Heatmap) 생성
plt.figure(figsize=(12, 8))
sns.heatmap(
    interaction_matrix, 
    annot=True, 
    fmt="d", 
    cmap="YlGnBu", 
    linewidths=.5, 
    cbar_kws={'label': 'Number of Sequences'}
)
plt.title('Virus-Host Interaction Matrix in Korea (Top Species)', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Host Scientific Name', fontsize=12)
plt.ylabel('Virus Organism', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_virus_host_heatmap.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved high-res heatmap figure: {fig_path}")


# ==========================================
# FILE: analyze_virus_tsne_only.py
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

# 1. 전처리
df['Sequence_Length'] = pd.to_numeric(df['Sequence_Length'], errors='coerce').fillna(df['Sequence_Length'].median())
df['Log_Length'] = np.log1p(df['Sequence_Length'])

features_for_cluster = pd.get_dummies(df[['Genome_Type', 'Genome_Segmentation', 'Vector_Transmission_Mode']], drop_first=True)
features_for_cluster['Log_Length'] = df['Log_Length']

features_clean = features_for_cluster.dropna()
scaler = StandardScaler()
scaled_features = scaler.fit_transform(features_clean)

# 2. t-SNE 수행
tsne = TSNE(n_components=2, perplexity=30, random_state=42, max_iter=1000)
tsne_coords = tsne.fit_transform(scaled_features)

result_df = df.loc[features_clean.index].copy()
result_df['tSNE1'] = tsne_coords[:, 0]
result_df['tSNE2'] = tsne_coords[:, 1]

# 3. t-SNE 단독 시각화 (논문 규격 최적화)
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=result_df, x='tSNE1', y='tSNE2', hue='Genome_Type', 
    palette='Set1', alpha=0.8, s=60, edgecolor='black', linewidth=0.3
)

plt.title('t-SNE Clustering of Plant Viruses by Genome Architecture & Ecology', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('t-SNE Dimension 1', fontsize=12)
plt.ylabel('t-SNE Dimension 2', fontsize=12)
plt.legend(title='Genome Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10, title_fontsize=11)

plt.tight_layout()
output_fig = "virus_tsne_standalone.png"
plt.savefig(output_fig, dpi=300)
plt.close()

print(f"[완료] t-SNE 단독 고해상도 그래프 저장 성공: {output_fig}")


# ==========================================
# FILE: analyze_viruses.py
# ==========================================

from Bio import SeqIO
import os
from collections import Counter

gbk_file = "../data/korea_plant_viruses_raw.gbk"

if not os.path.exists(gbk_file):
    print(f"[Error] File not found: {gbk_file}")
    exit()

print("Parsing GenBank records and analyzing taxonomy distributions...")

total_records = 0
taxa_counter = Counter()

# GenBank 파일 파싱
for record in SeqIO.parse(gbk_file, "genbank"):
    total_records += 1
    # 소스 피처(source feature)에서 유기체 이름 및 분류 정보 추출
    organism = record.annotations.get("organism", "Unknown")
    taxonomy = record.annotations.get("taxonomy", [])
    
    # 주요 과(Family) 또는 그룹 매핑을 위한 분류군 확인
    family = "Unclassified / Other"
    for tax in taxonomy:
        if tax.endswith("viridae") or tax.endswith("viroidae") or tax in ["Geminiviridae", "Caulimoviridae", "Nanoviridae", "Potyviridae", "Betaflexiviridae", "Alphaflexiviridae", "Secoviridae", "Tombusviridae", "Partitiviridae", "Rhabdoviridae", "Tospoviridae", "Pospiviroidae", "Avsunviroidae"]:
            family = tax
            break
            
    taxa_counter[family] += 1

print("\n" + "="*50)
print(f" Total Analyzed Records: {total_records}")
print("="*50)
print(f"{'Taxonomy Group / Family':<30} | {'Count':<8} | {'Percentage':<8}")
print("-" * 50)

for family, count in taxa_counter.most_common(20):
    percentage = (count / total_records) * 100 if total_records > 0 else 0
    print(f"{family:<30} | {count:<8} | {percentage:.2f}%")

print("="*50)


# ==========================================
# FILE: analyze_yearly_trend_comparison.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

global_csv = "../data/global_data/global_plant_viruses_metadata.csv"
korea_csv = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading metadata for yearly trend comparison...")
df_global = pd.read_csv(global_csv)
df_korea = pd.read_csv(korea_csv)

# 연도 추출 함수
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        year_val = int(match.group(1))
        # 유효한 연도 범위 설정 (1980 ~ 2026)
        if 1980 <= year_val <= 2026:
            return year_val
    return None

print("Extracting years for global and Korea datasets...")
df_global['Year'] = df_global['Collection_Date'].apply(extract_year)
df_korea['Year'] = df_korea['Collection_Date'].apply(extract_year)

# 연도별 집계 (1990년 이후 집중)
global_yearly = df_global[df_global['Year'] >= 1990]['Year'].value_counts().sort_index().reset_index()
global_yearly.columns = ['Year', 'Global_Count']

korea_yearly = df_korea[df_korea['Year'] >= 1990]['Year'].value_counts().sort_index().reset_index()
korea_yearly.columns = ['Year', 'Korea_Count']

# 데이터 병합
trend_df = pd.merge(korea_yearly, global_yearly, on='Year', how='outer').fillna(0)
trend_df = trend_df.sort_values('Year').reset_index(drop=True)

# CSV 저장
trend_path = os.path.join(results_dir, "global_vs_korea_yearly_trend_summary.csv")
trend_df.to_csv(trend_path, index=False)
print(f"-> Saved yearly trend summary CSV: {trend_path}")

print("\nYearly Trend Summary Preview:")
print(trend_df.tail(15))

# 시각화: 이중 축(Dual-axis) 또는 서브플롯을 통한 시계열 비교 라인 그래프
fig, ax1 = plt.subplots(figsize=(12, 6))

# 글로벌 트렌드 (왼쪽 축)
color = '#2ca02c'
ax1.set_xlabel('Year', fontsize=12)
ax1.set_ylabel('Global Number of Sequences', color=color, fontsize=12)
line1 = ax1.plot(trend_df['Year'], trend_df['Global_Count'], color=color, marker='o', linewidth=2.5, label='Global')
ax1.tick_params(axis='y', labelcolor=color)
ax1.grid(True, linestyle='--', alpha=0.4)

# 한국 트렌드 (오른쪽 축)
ax2 = ax1.twinx()
color = '#1f77b4'
ax2.set_ylabel('South Korea Number of Sequences', color=color, fontsize=12)
line2 = ax2.plot(trend_df['Year'], trend_df['Korea_Count'], color=color, marker='s', linewidth=2.5, label='South Korea')
ax2.tick_params(axis='y', labelcolor=color)
ax2.grid(False) # 격자 겹침 방지

# 범례 통합
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, loc='upper left', frameon=True)

plt.title('Temporal Dynamics of Plant Virus Submissions: South Korea vs. Global (1990–2026)', fontsize=14, fontweight='bold', pad=15)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "global_vs_korea_yearly_trend.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved comparative yearly trend high-res figure: {fig_path}")


# ==========================================
# FILE: core_focused_analysis.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("="*60)
print(" 🔬 식물 바이러스 핵심 생물학적 요인 집중 분석")
print("="*60)

df = pd.read_csv(master_csv)

# 1. 주요 기주 작물 상위 8개 및 증상 분포 요약
print("\n[1] 주요 기주 작물별 검출 빈도 상위 8개")
top_hosts = df['Host_Plant'].value_counts().head(8)
print(top_hosts)

# 2. 유전체 유형(Genome Type)에 따른 서열 길이 통계 요약
print("\n[2] 유전체 유형별 서열 길이(Sequence Length) 통계 (bp)")
length_stats = df.groupby('Genome_Type')['Sequence_Length'].agg(['count', 'mean', 'median', 'min', 'max'])
print(length_stats)

# 3. 매개체 전파 방식(Vector Transmission)별 전장 유전체 확보율 요약
print("\n[3] 매개체 전파 방식별 전장 유전체(Complete Genome) 비율 (%)")
vector_completeness = pd.crosstab(df['Vector_Transmission_Mode'], df['Complete_Genome_Status'], normalize='index') * 100
print(vector_completeness.round(2))

# 4. 시각화 1: 유전체 유형별 서열 길이 박스플롯 (정제 버전)
plt.figure(figsize=(10, 6))
sns.boxplot(x='Genome_Type', y='Sequence_Length', data=df, palette='Set2', showfliers=False) # 이상치 제외로 가독성 극대화
plt.title('Sequence Length Distribution by Genome Type (Excluding Outliers)', fontsize=14, fontweight='bold')
plt.xlabel('Genome Type', fontsize=12)
plt.ylabel('Sequence Length (bp)', fontsize=12)
plt.xticks(rotation=30)
plt.tight_layout()
plt.savefig('clean_genome_length_box.png', dpi=300)
plt.close()
print(" -> 정제된 박스플롯 저장 완료: clean_genome_length_box.png")

# 5. 시각화 2: 기주 작물 상위 8개와 병원성 증상 간의 교차 히트맵 (지역 변수 제외)
plt.figure(figsize=(12, 6))
top_host_list = top_hosts.index
host_symptom = pd.crosstab(df[df['Host_Plant'].isin(top_host_list)]['Host_Plant'], df['Pathogenicity_Symptoms'])
sns.heatmap(host_symptom, annot=True, fmt='d', cmap='Blues', cbar=True)
plt.title('Pathogenicity Symptoms Heatmap by Major Host Plants', fontsize=14, fontweight='bold')
plt.xlabel('Pathogenicity Symptoms', fontsize=12)
plt.ylabel('Host Plant', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('clean_host_symptom_heatmap.png', dpi=300)
plt.close()
print(" -> 깔끔한 기주-증상 히트맵 저장 완료: clean_host_symptom_heatmap.png")

print("\n" + "="*60)
print(" ✨ 핵심 집중 분석 및 클린 시각화 생성 완료!")
print("="*60)


# ==========================================
# FILE: download_and_clean_viruses.py
# ==========================================

from Bio import Entrez, SeqIO
import os
import time
import re
from urllib.error import HTTPError
from collections import Counter

Entrez.email = "wonkyong.cho@example.com"

# 1. 식물 바이러스, 바이로이드, 그리고 주요 식물 바이러스 과(Family)로 엄격히 제한된 검색식
search_term = (
    '(plant virus[Organism] OR viroid[Organism] OR '
    'Geminiviridae[Organism] OR Caulimoviridae[Organism] OR Nanoviridae[Organism] OR '
    'Potyviridae[Organism] OR Alphaflexiviridae[Organism] OR Betaflexiviridae[Organism] OR '
    'Secoviridae[Organism] OR Tombusviridae[Organism] OR Cucumovirus[Organism] OR Tobamovirus[Organism] OR '
    'Partitiviridae[Organism] OR Rhabdoviridae[Organism] OR Tospoviridae[Organism] OR '
    'Pospiviroidae[Organism] OR Avsunviroidae[Organism] OR '
    'Tymoviridae[Organism] OR Luteoviridae[Organism] OR Bromoviridae[Organism] OR Closteroviridae[Organism])'
)

print("Searching GenBank for plant-specific viruses and viroids...")
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=200000, usehistory="y")
record = Entrez.read(handle)
total_count = int(record["Count"])
print(f"Total matching sequences found worldwide: {total_count}")

webenv = record["WebEnv"]
query_key = record["QueryKey"]

output_file = "../data/korea_plant_viruses_clean.gbk"
os.makedirs("../data", exist_ok=True)

batch_size = 500
saved_count = 0

# 배제할 비식물성 키워드 (동물, 어류, 인간 바이러스 등)
exclude_virus_keywords = [
    "hemorrhagic septicemia", "rabies", "lyssavirus", "influenza", 
    "coronavirus", "hepatitis", "hiv", "herpes", "adenovirus", "bacteriophage", "phage"
]

print("Downloading, strictly filtering, and saving Korea-related plant pathogen entries...")
with open(output_file, "w") as out_f:
    for start in range(0, total_count, batch_size):
        success = False
        for attempt in range(3):
            try:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                    rettype="gb",
                    retmode="text"
                )
                records_text = fetch_handle.read()
                fetch_handle.close()
                
                records_list = records_text.split("//\n")
                for rec in records_list:
                    if len(rec.strip()) > 0:
                        # 1) 한국 관련 지명 매칭 확인
                        if re.search(r'\b(KOREA|SOUTH\s+KOREA|REPUBLIC\s+OF\s+KOREA|SEOUL|JEJU|GANGWON|GYEONNGI|BUSAN|DAEGU|INCHEON|GWANGJU|DAEJEON|ULSAN|CHUNGBUK|CHUNGNAM|JEONBUK|JEONNAM|GYEONGBUK|GYEONGNAM)\b', rec, re.IGNORECASE):
                            
                            # 2) 비식물성 바이러스 키워드가 포함되어 있으면 저장하지 않음
                            rec_lower = rec.lower()
                            if any(ex in rec_lower for ex in exclude_virus_keywords):
                                continue
                                
                            out_f.write(rec + "//\n")
                            saved_count += 1
                            
                print(f"Processed {min(start + batch_size, total_count)} / {total_count} | Cleaned Korea matches: {saved_count}")
                success = True
                time.sleep(0.3)
                break
            except (HTTPError, Exception) as e:
                print(f"  [Warning] Retry {attempt+1} due to error: {e}")
                time.sleep(2)
        
        if not success:
            print(f"  [Error] Failed to process batch starting at {start}")

print(f"Strict filtering and download complete! Total {saved_count} clean records saved to {output_file}")

# --- 곧바로 정제된 데이터 통계 분석 실행 ---
print("\nRunning statistics on the newly cleaned dataset...")
species_counter = Counter()
host_counter = Counter()
year_counter = Counter()

exclude_hosts = ["apis mellifera", "bemisia tabaci", "honey bee", "bee", "mouse", "human"]
host_mapping = {
    "chrysanthemum": "Chrysanthemum spp.",
    "tomato": "Solanum lycopersicum",
    "pepper": "Capsicum annuum",
    "soybean": "Glycine max",
    "melon": "Cucumis melo"
}

for record in SeqIO.parse(output_file, "genbank"):
    species_counter[record.annotations.get("organism", "Unknown")] += 1
    
    # 기주 정제
    host_name = "Not Specified"
    for feature in record.features:
        if feature.type == "source":
            h_list = feature.qualifiers.get("host", [None])
            if h_list[0]:
                raw_host = h_list[0].strip().lower()
                if any(ex in raw_host for ex in exclude_hosts):
                    break
                cleaned_host = raw_host
                for k, v in host_mapping.items():
                    if k in raw_host:
                        cleaned_host = v.lower()
                        break
                host_name = cleaned_host.capitalize()
                break
    host_counter[host_name] += 1
    
    # 연도 추출
    date_str = record.annotations.get("date", "")
    year_match = re.search(r'\b(19\d{2}|20\d{2})\b', date_str)
    if year_match:
        year_counter[year_match.group(1)] += 1
    else:
        year_counter["Unknown"] += 1

print("\n" + "="*60)
print("[ Top 10 Clean Plant Viral & Viroid Species ]")
print("-" * 50)
for sp, cnt in species_counter.most_common(10):
    print(f"{sp:<35} | {cnt:<6}")

print("\n[ Top 10 Clean Plant Hosts ]")
print("-" * 50)
for h, cnt in host_counter.most_common(10):
    print(f"{h:<35} | {cnt:<6}")
print("="*60)


# ==========================================
# FILE: download_global_viruses.py
# ==========================================

from Bio import Entrez, SeqIO
import os
import time
import pandas as pd
from urllib.error import HTTPError

Entrez.email = "wonkyong.cho@example.com"

# 전 세계 식물 바이러스 및 비로이드, 주요 과 명칭을 모두 포괄하는 검색식
search_term = '(plant virus[Organism] OR viroid[Organism] OR Potyviridae[Organism] OR Alphaflexiviridae[Organism] OR Betaflexiviridae[Organism] OR Secoviridae[Organism] OR Tombusviridae[Organism] OR Cucumovirus[Organism] OR Tobamovirus[Organism])'

print("Searching GenBank for global plant virus and viroid sequences...")
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=150000, usehistory="y")
record = Entrez.read(handle)
total_count = int(record["Count"])
print(f"Total global matching sequences found: {total_count}")

webenv = record["WebEnv"]
query_key = record["QueryKey"]

# 글로벌 데이터 전용 저장 경로 설정
global_dir = "../data/global_data"
os.makedirs(global_dir, exist_ok=True)
output_gbk = os.path.join(global_dir, "global_plant_viruses_raw.gbk")
output_csv = os.path.join(global_dir, "global_plant_viruses_metadata.csv")

batch_size = 500
downloaded_records = 0

print("Downloading global raw GenBank file...")
with open(output_gbk, "w") as out_f:
    for start in range(0, total_count, batch_size):
        success = False
        for attempt in range(3):
            try:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                    rettype="gb",
                    retmode="text"
                )
                records_text = fetch_handle.read()
                fetch_handle.close()
                
                out_f.write(records_text)
                downloaded_records += records_text.count("//\n")
                
                print(f"Downloaded batch {start} to {min(start + batch_size, total_count)} / {total_count}")
                success = True
                time.sleep(0.3)
                break
            except (HTTPError, Exception) as e:
                print(f"  [Warning] Retry {attempt+1} due to error: {e}")
                time.sleep(2)
        
        if not success:
            print(f"  [Error] Failed to process batch starting at {start}")

print(f"Global raw download complete! Total {downloaded_records} records saved to {output_gbk}")

# 다운로드 직후 글로벌 메타데이터 파싱 수행
print("Parsing global metadata into CSV...")
metadata = []
for record in SeqIO.parse(output_gbk, "genbank"):
    acc = record.id
    description = record.description
    organism = record.annotations.get("organism", "")
    
    host = "Unknown"
    country = "Unknown"
    collection_date = "Unknown"
    isolation_source = "Unknown"
    
    for feature in record.features:
        if feature.type == "source":
            qual = feature.qualifiers
            host = qual.get("host", [host])[0]
            country = qual.get("country", [country])[0]
            collection_date = qual.get("collection_date", [collection_date])[0]
            isolation_source = qual.get("isolation_source", [isolation_source])[0]
            break
            
    metadata.append({
        "Accession": acc,
        "Organism": organism,
        "Description": description,
        "Host": host,
        "Country": country,
        "Collection_Date": collection_date,
        "Isolation_Source": isolation_source,
        "Sequence_Length": len(record.seq)
    })

df_global = pd.DataFrame(metadata)
df_global.to_csv(output_csv, index=False)
print(f"Global metadata parsing complete! Saved to {output_csv}")


# ==========================================
# FILE: download_korea_viruses.py
# ==========================================

from Bio import Entrez
import os
import time
import re
from urllib.error import HTTPError

Entrez.email = "wonkyong.cho@example.com"

# 식물 바이러스 및 비로이드, 그리고 주요 분류군 키워드를 총망라한 확장 검색식
search_term = '(plant virus[Organism] OR viroid[Organism] OR Potyviridae[Organism] OR Alphaflexiviridae[Organism] OR Betaflexiviridae[Organism] OR Secoviridae[Organism] OR Tombusviridae[Organism] OR Cucumovirus[Organism] OR Tobamovirus[Organism])'

print("Searching GenBank for comprehensive plant virus and viroid sequences...")
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=150000, usehistory="y")
record = Entrez.read(handle)
total_count = int(record["Count"])
print(f"Total matching sequences found worldwide: {total_count}")

webenv = record["WebEnv"]
query_key = record["QueryKey"]

output_file = "../data/korea_plant_viruses_raw.gbk"
os.makedirs("../data", exist_ok=True)

batch_size = 500
saved_count = 0

print("Downloading and filtering for Korea-related entries using WebEnv...")
with open(output_file, "w") as out_f:
    for start in range(0, total_count, batch_size):
        success = False
        for attempt in range(3):
            try:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                    rettype="gb",
                    retmode="text"
                )
                records_text = fetch_handle.read()
                fetch_handle.close()
                
                # 정규식을 이용해 한국 관련 지명/국가명이 포함된 레코드 정밀 추출
                records_list = records_text.split("//\n")
                for rec in records_list:
                    if len(rec.strip()) > 0:
                        # 대소문자 무시하고 국가 및 주요 지역명 매칭 확인
                        if re.search(r'\b(KOREA|SOUTH\s+KOREA|REPUBLIC\s+OF\s+KOREA|SEOUL|JEJU|GANGWON|GYEONNGI|BUSAN|DAEGU|INCHEON|GWANGJU|DAEJEon|ULSAN|CHUNGBUK|CHUNGNAM|JEONBUK|JEONNAM|GYEONGBUK|GYEONGNAM)\b', rec, re.IGNORECASE):
                            out_f.write(rec + "//\n")
                            saved_count += 1
                            
                print(f"Processed {min(start + batch_size, total_count)} / {total_count} | Accumulated Korea matches: {saved_count}")
                success = True
                time.sleep(0.3)
                print(f"Downloaded batch {start} to {min(start + batch_size, total_count)}")
                break
            except (HTTPError, Exception) as e:
                print(f"  [Warning] Retry {attempt+1} due to error: {e}")
                time.sleep(2)
        
        if not success:
            print(f"  [Error] Failed to process batch starting at {start}")

print(f"Download & Filtering complete! Total {saved_count} records saved to {output_file}")


# ==========================================
# FILE: download_viruses.py
# ==========================================

from Bio import Entrez
import os
import time
import re
from urllib.error import HTTPError

Entrez.email = "wonkyong.cho@example.com"

# 식물 DNA 바이러스, RNA 바이러스, 바이로이드 및 관련 모든 분류군을 포괄하는 확장 검색식
search_term = (
    '(plant virus[Organism] OR viroid[Organism] OR '
    'Geminiviridae[Organism] OR Caulimoviridae[Organism] OR Nanoviridae[Organism] OR '
    'Potyviridae[Organism] OR Alphaflexiviridae[Organism] OR Betaflexiviridae[Organism] OR '
    'Secoviridae[Organism] OR Tombusviridae[Organism] OR Cucumovirus[Organism] OR Tobamovirus[Organism] OR '
    'Partitiviridae[Organism] OR Rhabdoviridae[Organism] OR Tospoviridae[Organism] OR '
    'Pospiviroidae[Organism] OR Avsunviroidae[Organism] OR '
    'Tymoviridae[Organism] OR Luteoviridae[Organism] OR Bromoviridae[Organism] OR Closteroviridae[Organism])'
)

print("Searching GenBank for comprehensive plant DNA/RNA viruses and viroids...")
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=200000, usehistory="y")
record = Entrez.read(handle)
total_count = int(record["Count"])
print(f"Total matching sequences found worldwide: {total_count}")

webenv = record["WebEnv"]
query_key = record["QueryKey"]

output_file = "../data/korea_plant_viruses_raw.gbk"
os.makedirs("../data", exist_ok=True)

batch_size = 500
saved_count = 0

print("Downloading and filtering for Korea-related entries using WebEnv...")
with open(output_file, "w") as out_f:
    for start in range(0, total_count, batch_size):
        success = False
        for attempt in range(3):
            try:
                fetch_handle = Entrez.efetch(
                    db="nucleotide",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                    rettype="gb",
                    retmode="text"
                )
                records_text = fetch_handle.read()
                fetch_handle.close()
                
                records_list = records_text.split("//\n")
                for rec in records_list:
                    if len(rec.strip()) > 0:
                        if re.search(r'\b(KOREA|SOUTH\s+KOREA|REPUBLIC\s+OF\s+KOREA|SEOUL|JEJU|GANGWON|GYEONNGI|BUSAN|DAEGU|INCHEON|GWANGJU|DAEJEON|ULSAN|CHUNGBUK|CHUNGNAM|JEONBUK|JEONNAM|GYEONGBUK|GYEONGNAM)\b', rec, re.IGNORECASE):
                            out_f.write(rec + "//\n")
                            saved_count += 1
                            
                print(f"Processed {min(start + batch_size, total_count)} / {total_count} | Accumulated Korea matches: {saved_count}")
                success = True
                time.sleep(0.3)
                break
            except (HTTPError, Exception) as e:
                print(f"  [Warning] Retry {attempt+1} due to error: {e}")
                time.sleep(2)
        
        if not success:
            print(f"  [Error] Failed to process batch starting at {start}")

print(f"Download & Filtering complete! Total {saved_count} records saved to {output_file}")


# ==========================================
# FILE: export_all_lists.py
# ==========================================

from Bio import SeqIO
import os
from collections import Counter

input_file = "../data/korea_plant_viruses_clean.gbk"
output_txt = "../data/complete_virus_and_host_lists.txt"

if not os.path.exists(input_file):
    print(f"[Error] Cleaned file not found: {input_file}")
    exit()

print("Extracting all unique viral species and plant hosts from the cleaned dataset...")

species_counter = Counter()
host_counter = Counter()

for record in SeqIO.parse(input_file, "genbank"):
    # 바이러스 종명 수집
    org = record.annotations.get("organism", "Unknown")
    species_counter[org] += 1
    
    # 기주 정보 수집
    host_name = "Not Specified"
    for feature in record.features:
        if feature.type == "source":
            h_list = feature.qualifiers.get("host", [None])
            if h_list[0]:
                host_name = h_list[0].strip()
                break
    host_counter[host_name] += 1

# 결과 파일로 저장
with open(output_txt, "w", encoding="utf-8") as f:
    f.write("="*60 + "\n")
    f.write(" ALL UNIQUE PLANT VIRAL & VIROID SPECIES (Sorted by Frequency)\n")
    f.write("="*60 + "\n")
    for sp, cnt in species_counter.most_common():
        f.write(f"{sp} ({cnt}건)\n")
        
    f.write("\n" + "="*60 + "\n")
    f.write(" ALL UNIQUE PLANT HOSTS (Sorted by Frequency)\n")
    f.write("="*60 + "\n")
    for h, cnt in host_counter.most_common():
        f.write(f"{h} ({cnt}건)\n")

print(f"\n모든 바이러스 종명과 기주명 리스트가 성공적으로 저장되었습니다.")
print(f"저장 경로: {os.path.abspath(output_txt)}")
print(f"총 고유 바이러스 종 수: {len(species_counter)}개")
print(f"총 고유 기주 종류 수: {len(host_counter)}개")


# ==========================================
# FILE: export_metadata_table.py
# ==========================================

from Bio import SeqIO
import os
import re
import csv

input_file = "../data/korea_plant_viruses_clean.gbk"
output_csv = "../data/plant_viruses_metadata_table.csv"

if not os.path.exists(input_file):
    print(f"[Error] Cleaned file not found: {input_file}")
    exit()

print("Extracting Accession numbers and all associated metadata into a structured table...")

# 제외할 비식물성 키워드 (혹시 모를 잔여 항목 방지)
non_plant_terms = [
    "hemorrhagic septicemia", "rabies", "lyssavirus", "influenza", 
    "coronavirus", "hepatitis", "hiv", "herpes", "adenovirus", 
    "bacteriophage", "phage", "cryptosporidium", "penicillium", 
    "trichoderma", "aspergillus", "oncorhynchus", "salmo", "cyprinus", 
    "bos taurus", "bovine", "apis mellifera", "bemisia tabaci"
]

row_count = 0

with open(output_csv, "w", newline="", encoding="utf-8-sig") as csv_file:
    writer = csv.writer(csv_file)
    # CSV 헤더 작성 (Accession 중심의 메타데이터 필드)
    writer.writerow(["Accession", "Version", "Viral_Species", "Host_Plant", "Collection_Year", "Taxonomy_Lineage"])
    
    for record in SeqIO.parse(input_file, "genbank"):
        org = record.annotations.get("organism", "Unknown")
        org_lower = org.lower()
        
        # 비식물성 항목 2차 차단
        if any(term in org_lower for term in non_plant_terms):
            continue
            
        # Accession 및 Version 추출
        accession = record.id
        version = record.name # 또는 record.annotations의 accession 정보 활용
        
        # 기주 정보 추출
        host_name = "Not Specified"
        is_non_plant_host = False
        for feature in record.features:
            if feature.type == "source":
                h_list = feature.qualifiers.get("host", [None])
                if h_list[0]:
                    raw_host = h_list[0].strip()
                    if any(term in raw_host.lower() for term in ["oncorhynchus", "salmo", "cyprinus", "bos taurus", "apis", "bemisia"]):
                        is_non_plant_host = True
                        break
                    host_name = raw_host
                    break
        if is_non_plant_host:
            continue
            
        # 연도 추출
        date_str = record.annotations.get("date", "")
        year_match = re.search(r'\b(19\d{2}|20\d{2})\b', date_str)
        year = year_match.group(1) if year_match else "Unknown"
        
        # Taxonomy 라인 추출
        taxonomy = "; ".join(record.annotations.get("taxonomy", []))
        
        # 행 기록
        writer.writerow([accession, record.name, org, host_name, year, taxonomy])
        row_count += 1

print(f"\n메타데이터 추출 완료!")
print(f"총 {row_count}개의 유효한 식물 바이러스/바이로이드 서열 레코드가 저장되었습니다.")
print(f"결과 CSV 파일 경로: {os.path.abspath(output_csv)}")


# ==========================================
# FILE: figure2.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 6), dpi=300)

# [Figure 2-A] 주요 벡터 전파 경로 분포
ax1 = plt.subplot(1, 2, 1)
vector_counts = df['Vector_Transmission_Mode'].value_counts().head(6).reset_index()
vector_counts.columns = ['Vector_Mode', 'Count']
sns.barplot(data=vector_counts, x='Count', y='Vector_Mode', hue='Vector_Mode', palette='Oranges_r', legend=False, ax=ax1)
ax1.set_title('(A) Top Vector Transmission Modes', fontweight='bold', fontsize=12)
ax1.set_xlabel('Number of Records', fontweight='bold')
ax1.set_ylabel('Transmission Mode', fontweight='bold')

# [Figure 2-B] 연도별 데이터 수집 추이 (Collection Year Trend)
ax2 = plt.subplot(1, 2, 2)
year_counts = df['Collection_Year'].value_counts().sort_index().reset_index()
year_counts.columns = ['Year', 'Count']
sns.lineplot(data=year_counts, x='Year', y='Count', marker='o', color='teal', linewidth=2.5, ax=ax2)
ax2.set_title('(B) Annual Trend of Virus Records', fontweight='bold', fontsize=12)
ax2.set_xlabel('Collection Year', fontweight='bold')
ax2.set_ylabel('Number of Records', fontweight='bold')

plt.tight_layout()
plt.savefig('plant_viruses_figure2_dynamics.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 두 번째 논문 도판이 'plant_viruses_figure2_dynamics.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: filter_existing_data.py
# ==========================================

from Bio import SeqIO
from collections import Counter
import os
import re

input_file = "../data/korea_plant_viruses_raw.gbk"
output_file = "../data/korea_plant_viruses_clean.gbk"

if not os.path.exists(input_file):
    print(f"[Error] Original file not found: {input_file}")
    exit()

print("Filtering existing raw data to remove non-plant pathogens...")

exclude_virus_keywords = [
    "hemorrhagic septicemia", "rabies", "lyssavirus", "influenza", 
    "coronavirus", "hepatitis", "hiv", "herpes", "adenovirus", "bacteriophage", "phage"
]

exclude_hosts = ["apis mellifera", "bemisia tabaci", "honey bee", "bee", "mouse", "human"]
host_mapping = {
    "chrysanthemum": "Chrysanthemum spp.",
    "tomato": "Solanum lycopersicum",
    "pepper": "Capsicum annuum",
    "soybean": "Glycine max",
    "melon": "Cucumis melo"
}

total_records = 0
saved_count = 0
species_counter = Counter()
host_counter = Counter()
year_counter = Counter()

with open(output_file, "w") as out_f:
    for record in SeqIO.parse(input_file, "genbank"):
        total_records += 1
        record_text = record.format("genbank")
        org_lower = record.annotations.get("organism", "").lower()
        
        # 1. 비식물성 바이러스 키워드 포함 시 제외
        if any(ex in org_lower or ex in record_text.lower() for ex in exclude_virus_keywords):
            continue
            
        # 2. 기주 정제 및 곤충/동물 기주 제외
        host_name = "Not Specified"
        skip_host = False
        for feature in record.features:
            if feature.type == "source":
                h_list = feature.qualifiers.get("host", [None])
                if h_list[0]:
                    raw_host = h_list[0].strip().lower()
                    if any(ex in raw_host for ex in exclude_hosts):
                        skip_host = True
                        break
                    cleaned_host = raw_host
                    for k, v in host_mapping.items():
                        if k in raw_host:
                            cleaned_host = v.lower()
                            break
                    host_name = cleaned_host.capitalize()
                    break
        if skip_host:
            continue
            
        out_f.write(record_text)
        saved_count += 1
        
        species_counter[record.annotations.get("organism", "Unknown")] += 1
        host_counter[host_name] += 1
        
        # 연도 추출
        date_str = record.annotations.get("date", "")
        year_match = re.search(r'\b(19\d{2}|20\d{2})\b', date_str)
        if year_match:
            year_counter[year_match.group(1)] += 1
        else:
            year_counter["Unknown"] += 1

print(f"\nFiltering complete!")
print(f"Total Raw Records: {total_records} -> Cleaned Plant Records: {saved_count}")
print(f"Saved to: {output_file}")

print("\n" + "="*60)
print("[ Top 10 Clean Plant Viral & Viroid Species ]")
print("-" * 50)
for sp, cnt in species_counter.most_common(10):
    print(f"{sp:<35} | {cnt:<6}")

print("\n[ Top 10 Clean Plant Hosts ]")
print("-" * 50)
for h, cnt in host_counter.most_common(10):
    print(f"{h:<35} | {cnt:<6}")
print("="*60)


# ==========================================
# FILE: fix_cohort_spacing.py
# ==========================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 코호트 분석 데이터 요약
cohort_summary = pd.DataFrame({
    'Diagnostic_Cohort': ['Traditional Era (~2014)', 'NGS Metagenomics Era (2015~)'],
    'Success_Rate_%': [31.8, 51.0]
})

# 캔버스 크기를 좁게 설정하여 바 간격 단축
plt.figure(figsize=(7, 6))

ax = sns.barplot(
    data=cohort_summary, x='Diagnostic_Cohort', y='Success_Rate_%', 
    palette=['#4682B4', '#2E8B57'], edgecolor='black',
    width=0.45  # 바 두께 조정
)

plt.title('Cohort Analysis: Complete Genome Recovery Rate by Diagnostic Era', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Diagnostic Technology Cohort', fontsize=11)
plt.ylabel('Complete Genome Recovery Rate (%)', fontsize=11)
plt.ylim(0, 105)
plt.xlim(-0.6, 1.6)  # X축 범위를 좁혀 바들을 안쪽으로 모음

# 바 상단 데이터 레이블 부착
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(f"{height:.1f}%",
                    (p.get_x() + p.get_width() / 2., height),
                    ha='center', va='bottom', fontsize=11, fontweight='bold',
                    xytext=(0, 3), textcoords='offset points')

plt.tight_layout()
output_fig = "diagnostic_cohort_analysis.png"
plt.savefig(output_fig, dpi=300)
plt.close()

print(f"[완료] 간격이 좁혀진 코호트 분석 그래프 저장 성공: {output_fig}")


# ==========================================
# FILE: fix_global_plot.py
# ==========================================

from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

gbk_file = "../data/global_data/global_plant_viruses_raw.gbk"

print("Directly parsing GenBank file for country extraction...")
countries = []

for record in SeqIO.parse(gbk_file, "genbank"):
    country = "Unknown"
    for feature in record.features:
        if feature.type == "source":
            if "country" in feature.qualifiers:
                country = feature.qualifiers["country"][0]
                break
    countries.append(country)

df_temp = pd.DataFrame({"Country_Raw": countries})
print(f"Total records parsed: {len(df_temp)}")
print(f"Non-unknown countries found: {len(df_temp[df_temp['Country_Raw'] != 'Unknown'])}")

# 국가명 표준화 함수
def clean_country(val):
    if pd.isna(val) or val == "Unknown":
        return "Unknown"
    val_str = str(val).strip()
    if ":" in val_str:
        val_str = val_str.split(":")[0].strip()
    
    val_lower = val_str.lower()
    if any(k in val_lower for k in ["korea", "republic of korea", "south korea"]):
        return "South Korea"
    elif "china" in val_lower:
        return "China"
    elif "usa" in val_lower or "united states" in val_lower or "america" in val_lower or "u.s." in val_lower:
        return "USA"
    elif "japan" in val_lower:
        return "Japan"
    elif "germany" in val_lower:
        return "Germany"
    elif "india" in val_lower:
        return "India"
    elif "italy" in val_lower:
        return "Italy"
    elif "spain" in val_lower:
        return "Spain"
    elif "france" in val_lower:
        return "France"
    elif "brazil" in val_lower:
        return "Brazil"
    elif "united kingdom" in val_lower or "uk" in val_lower or "england" in val_lower:
        return "United Kingdom"
    elif "iran" in val_lower:
        return "Iran"
    elif "turkey" in val_lower or "türkiye" in val_lower:
        return "Turkey"
    elif "canada" in val_lower:
        return "Canada"
    elif "australia" in val_lower:
        return "Australia"
        
    return val_str.title()

df_temp['Clean_Country'] = df_temp['Country_Raw'].apply(clean_country)

# 상위 10개국 집계
valid_df = df_temp[~df_temp['Clean_Country'].isin(["Unknown", "nan", "None"])]
top_countries = valid_df['Clean_Country'].value_counts().head(10).reset_index()
top_countries.columns = ['Country', 'Sequence_Count']

print("Top Countries Extracted:")
print(top_countries)

# 그래프 그리기
if len(top_countries) > 0:
    plt.figure(figsize=(11, 6))
    sns.barplot(data=top_countries, x='Sequence_Count', y='Country', palette='Purples_r')
    plt.title('Top 10 Countries Reporting Plant Viruses and Viroids Globally', fontsize=14, fontweight='bold', pad=15)
    plt.xlabel('Number of Sequences', fontsize=12)
    plt.ylabel('Country', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, "global_top_countries.png"), dpi=300)
    plt.close()
    print("Successfully generated and saved: global_top_countries.png")
else:
    print("[Error] No valid country data found to plot.")


# ==========================================
# FILE: fix_host_scientific_names.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")

print("Loading metadata for strict scientific host cleaning...")
df = pd.read_csv(csv_file)

# 수동 동의어 및 일반명 -> 학명 표준화 딕셔너리
synonym_map = {
    "Soybean": "Glycine max",
    "Garlic": "Allium sativum",
    "Pepper": "Capsicum annuum",
    "Hot pepper": "Capsicum annuum",
    "Tomato": "Solanum lycopersicum",
    "Peach": "Prunus persica",
    "Radish": "Raphanus sativus"
}

def strict_scientific_host(host_str):
    if pd.isna(host_str):
        return None
    
    clean_str = str(host_str).strip()
    clean_lower = clean_str.lower()
    
    # 1. 딕셔너리 매핑 확인
    for key, val in synonym_map.items():
        if key.lower() in clean_lower:
            return val
            
    # 2. 비식물 기주나 잘못 들어간 메타데이터 제외 처리
    exclude_keywords = ["apis", "honeybee", "insect", "human", "virus", "unknown", "none", "n/a"]
    if any(ex in clean_lower for ex in exclude_keywords):
        return None
        
    # 3. 두 단어 이상으로 이루어진 학명 형태(Genus species) 추출 시도
    parts = clean_str.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        # 특수문자 제거 후 반환
        import re
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
            
    return None

df['Strict_Host'] = df['Host'].apply(strict_scientific_host)

# 유효한 학명만 필터링 후 재집계
valid_df = df.dropna(subset=['Strict_Host'])
top_hosts = valid_df['Strict_Host'].value_counts().head(10).reset_index()
top_hosts.columns = ['Host_Scientific_Name', 'Sequence_Count']
top_hosts['Percentage'] = (top_hosts['Sequence_Count'] / len(df)) * 100

# 결과 저장
summary_path = os.path.join(results_dir, "korea_host_range_summary.csv")
top_hosts.to_csv(summary_path, index=False)
print("-> Updated & Saved strict scientific host summary:")
print(top_hosts)

# 고해상도 그래프 재생성
plt.figure(figsize=(11, 6))
sns.barplot(
    data=top_hosts, 
    x='Sequence_Count', 
    y='Host_Scientific_Name', 
    hue='Host_Scientific_Name',
    palette='Greens_r',
    legend=False
)
plt.title('Top 10 Plant Host Species (Strict Scientific Names) in Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Virus/Viroid Sequences', fontsize=12)
plt.ylabel('Host Scientific Name (Italicized style)', fontsize=12)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_top_hosts_scientific.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Updated high-res figure: {fig_path}")


# ==========================================
# FILE: fix_korea_virus_plot.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 스타일 설정
sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_top_viruses_summary.csv")

print("Loading Korean top viruses summary...")
df_viruses = pd.read_csv(csv_file)

# 바이러스 종 이름 정제 함수 (불필요한 긴 균주명이나 접두어 정리)
def clean_virus_name(name):
    if pd.isna(name):
        return "Unknown"
    # 일반적으로 물결표나 괄호, 긴 균주명 앞부분만 깔끔하게 추출하거나 정리
    name_str = str(name).strip()
    # 예: "Red clover necrotic mosaic virus isolate..." 같은 경우 종 이름만 남기거나 보기 좋게 축약
    return name_str

# 상위 종 이름 정제 적용
df_viruses['Clean_Name'] = df_viruses['Virus_Organism'].apply(clean_virus_name)

# 데이터 다시 집계 및 상위 10개 추출
top_clean = df_viruses.groupby('Clean_Name')['Count'].sum().reset_index()
top_clean = top_clean.sort_values(by='Count', ascending=False).head(10)

print("Top Cleaned Korean Viruses:")
print(top_clean)

# 논문용 고해상도 그래프 재생성
plt.figure(figsize=(11, 6))
sns.barplot(data=top_clean, x='Count', y='Clean_Name', palette='Blues_r')
plt.title('Top 10 Plant Virus and Viroid Species in Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Number of Sequences', fontsize=12)
plt.ylabel('Virus Species', fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, "korea_top_viruses.png"), dpi=300)
plt.close()

print("Successfully regenerated and updated: korea_top_viruses.png in results/figures/")


# ==========================================
# FILE: generate_Figure1_overview.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

# 1. Plotting style configuration for publication
sns.set_theme(style="whitegrid")
plt.rcParams.update({
    'font.sans-serif': 'DejaVu Sans', 
    'font.family': 'sans-serif',
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9
})

# 2. Locate the enriched metadata CSV file
csv_candidates = [
    "plant_viruses_metadata_enriched.csv",
    "../results/plant_viruses_metadata_enriched.csv",
    "korea_plant_viruses_metadata.csv",
    "../results/korea_plant_viruses_metadata.csv"
]

csv_file = None
for path in csv_candidates:
    if os.path.exists(path):
        csv_file = path
        break

if csv_file is None:
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".csv"):
                csv_file = os.path.join(root, file)
                break
        if csv_file: break

if csv_file is None:
    raise FileNotFoundError("[Error] 메타데이터 CSV 파일을 찾을 수 없습니다.")

print(f"Loading metadata from '{csv_file}' for Figure 1 generation...")
df = pd.read_csv(csv_file)
print("DataFrame columns:", df.columns.tolist())

# 3. Robust column mapping & data cleaning
genome_type_col = 'Genome_Type' if 'Genome_Type' in df.columns else None
length_col = 'Sequence_Length' if 'Sequence_Length' in df.columns else ('Length' if 'Length' in df.columns else None)

# Host column detection
possible_host_cols = ['Host_Plant', 'Host', 'Host_Species']
host_col_orig = next((col for col in possible_host_cols if col in df.columns), None)

synonym_map = {
    "Soybean": "Glycine max",
    "Garlic": "Allium sativum",
    "Pepper": "Capsicum annuum",
    "Hot pepper": "Capsicum annuum",
    "Tomato": "Solanum lycopersicum",
    "Peach": "Prunus persica",
    "Radish": "Raphanus sativus"
}

def strict_scientific_host(host_str):
    if pd.isna(host_str):
        return "Not Specified"
    clean_str = str(host_str).strip()
    clean_lower = clean_str.lower()
    for key, val in synonym_map.items():
        if key.lower() in clean_lower:
            return val
    if any(ex in clean_lower for ex in ["unknown", "none", "n/a", "not specified"]):
        return "Not Specified"
    parts = clean_str.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
    return clean_str.capitalize()

if host_col_orig:
    df['Clean_Host'] = df[host_col_orig].apply(strict_scientific_host)

# Virus species column auto-detection (Looking for text columns containing virus names)
virus_col = None
for col in df.columns:
    col_lower = col.lower()
    if any(k in col_lower for k in ['organism', 'virus', 'species', 'pathogen', 'name', 'viroid']):
        # Ensure it contains string data and not accession numbers
        sample_vals = df[col].dropna().astype(str).head(10).tolist()
        if any(not val.replace('.', '').isalnum() or 'virus' in val.lower() or 'viroid' in val.lower() for val in sample_vals):
            virus_col = col
            break

if virus_col is None:
    # Fallback: find any object column with high textual diversity
    for col in df.select_dtypes(include=['object']).columns:
        if df[col].nunique() > 5 and not any(acc in col.lower() for acc in ['acc', 'id', 'date', 'region', 'host']):
            virus_col = col
            break

print(f"Auto-detected Virus Species column: {virus_col}")

if virus_col:
    def clean_virus_name(name):
        if pd.isna(name):
            return "Unknown"
        name_str = str(name).strip()
        name_lower = name_str.lower()
        
        if "cucumber mosaic" in name_lower:
            return "Cucumber mosaic virus"
        elif "orthotospovirus tomatomaculae" in name_lower or "tomato spotted wilt" in name_lower:
            return "Tomato spotted wilt virus"
        elif "chrysanthemum stunt viroid" in name_lower:
            return "Chrysanthemum stunt viroid"
        elif "peach latent mosaic viroid" in name_lower:
            return "Peach latent mosaic viroid"
        elif "tomato yellow leaf curl" in name_lower:
            return "Tomato yellow leaf curl virus"
        elif "soybean mosaic" in name_lower:
            return "Soybean mosaic virus"
        
        return name_str[0].upper() + name_str[1:] if len(name_str) > 0 else name_str

    df['Clean_Virus'] = df[virus_col].apply(clean_virus_name)
else:
    # Fallback if still not found
    df['Clean_Virus'] = "Unknown Virus"

# 4. Generate Figure 1 Panel Grid (2x2) matching original figure1 structure
fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)

# [Panel A] Distribution of Viral Genome Types
if genome_type_col:
    genome_counts = df[genome_type_col].value_counts().reset_index()
    genome_counts.columns = ['Genome_Type', 'Count']
    sns.barplot(data=genome_counts, x='Count', y='Genome_Type', hue='Genome_Type', ax=axes[0, 0], palette='Blues_r', legend=False)
    axes[0, 0].set_title('A. Distribution of Viral Genome Types', fontweight='bold', pad=12)
    axes[0, 0].set_xlabel('Number of Records')
    axes[0, 0].set_ylabel('Genome Type')

# [Panel B] Sequence Length by Genome Type
if length_col and genome_type_col:
    df_filtered = df[df[length_col] <= 15000]
    sns.boxplot(data=df_filtered, x=length_col, y=genome_type_col, hue=genome_type_col, ax=axes[0, 1], palette='Set2', showfliers=False, legend=False)
    axes[0, 1].set_title('B. Sequence Length by Genome Type', fontweight='bold', pad=12)
    axes[0, 1].set_xlabel('Sequence Length (bp)')
    axes[0, 1].set_ylabel('Genome Type')

# [Panel C] Top 6 Host Plants
if 'Clean_Host' in df.columns:
    top_hosts = df['Clean_Host'].value_counts().head(6).reset_index()
    top_hosts.columns = ['Host', 'Count']
    sns.barplot(data=top_hosts, x='Count', y='Host', hue='Host', ax=axes[1, 0], palette='Purples_r', legend=False)
    axes[1, 0].set_title('C. Top 6 Host Plants', fontweight='bold', pad=12)
    axes[1, 0].set_xlabel('Number of Records')
    axes[1, 0].set_ylabel('Host Plant')

# [Panel D] Top 6 Viral/Viroid Species
if 'Clean_Virus' in df.columns:
    top_viruses = df['Clean_Virus'].value_counts().head(6).reset_index()
    top_viruses.columns = ['Virus', 'Count']
    sns.barplot(data=top_viruses, x='Count', y='Virus', hue='Virus', ax=axes[1, 1], palette='Greens_r', legend=False)
    axes[1, 1].set_title('D. Top 6 Viral/Viroid Species', fontweight='bold', pad=12)
    axes[1, 1].set_xlabel('Number of Records')
    axes[1, 1].set_ylabel('Viral Species')

# Layout optimization
plt.tight_layout(pad=3.5)

# Save output
output_image = "Figure1_overview.png"
plt.savefig(output_image, dpi=300, bbox_inches='tight')
plt.close()

print(f"[Success] Figure 1 successfully generated and saved as '{output_image}'.")


# ==========================================
# FILE: generate_Figure2_taxonomy.py
# ==========================================

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import squarify
import os

# 1. Plotting style configuration for publication
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.sans-serif': 'DejaVu Sans',
    'font.family': 'sans-serif',
    'axes.titlesize': 14,
    'axes.labelsize': 11
})

# 2. Locate the metadata CSV file
csv_candidates = [
    "plant_viruses_metadata_enriched.csv",
    "plant_viruses_metadata_table.csv",
    "korea_plant_viruses_metadata.csv",
    "../results/plant_viruses_metadata_enriched.csv",
    "../results/korea_plant_viruses_metadata.csv"
]

csv_file = None
for path in csv_candidates:
    if os.path.exists(path):
        csv_file = path
        break

if csv_file is None:
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".csv"):
                csv_file = os.path.join(root, file)
                break
        if csv_file: break

if csv_file is None:
    raise FileNotFoundError("[Error] 메타데이터 CSV 파일을 찾을 수 없습니다.")

print(f"Loading metadata from '{csv_file}'...")
df = pd.read_csv(csv_file)

# 3. Enhanced Family Mapping Dictionary to minimize 'Other Families'
family_mapping = {
    'potyvirus': 'Potyviridae', 'potyviridae': 'Potyviridae', 'soybean mosaic': 'Potyviridae',
    'turnip mosaic': 'Potyviridae', 'potato virus y': 'Potyviridae', 'pepper mottle': 'Potyviridae',
    'watermelon mosaic': 'Potyviridae', 'zucchini yellow mosaic': 'Potyviridae', 'bean common mosaic': 'Potyviridae',
    'plum pox': 'Potyviridae', 'onion yellow dwarf': 'Potyviridae', 'leek yellow stripe': 'Potyviridae',
    'cucumovirus': 'Bromoviridae', 'bromoviridae': 'Bromoviridae', 'cucumber mosaic': 'Bromoviridae',
    'tomato aspermy': 'Bromoviridae', 'alfalfa mosaic': 'Bromoviridae', 'brome mosaic': 'Bromoviridae',
    'tospovirus': 'Tospoviridae', 'tospoviridae': 'Tospoviridae', 'tomato spotted wilt': 'Tospoviridae',
    'orthotospovirus': 'Tospoviridae', 'capsicum chlorosis': 'Tospoviridae', 'impatiens necrotic spot': 'Tospoviridae',
    'groundnut bud necrosis': 'Tospoviridae',
    'geminiviridae': 'Geminiviridae', 'begomovirus': 'Geminiviridae', 'tomato yellow leaf curl': 'Geminiviridae',
    'papaya leaf curl': 'Geminiviridae', 'maff': 'Geminiviridae', 'ageratum yellow vein': 'Geminiviridae',
    'mungbean yellow mosaic': 'Geminiviridae', 'tomato leaf curl': 'Geminiviridae',
    'betaflexiviridae': 'Betaflexiviridae', 'carlavirus': 'Betaflexiviridae', 'foveavirus': 'Betaflexiviridae',
    'garlic virus': 'Betaflexiviridae', 'potato virus m': 'Betaflexiviridae', 'potato virus s': 'Betaflexiviridae',
    'apple stem grooving': 'Betaflexiviridae', 'cherry B': 'Betaflexiviridae',
    'alphaflexiviridae': 'Alphaflexiviridae', 'potexvirus': 'Alphaflexiviridae', 'potato virus x': 'Alphaflexiviridae',
    'white clover mosaic': 'Alphaflexiviridae', 'bamboo mosaic': 'Alphaflexiviridae',
    'pospiviroidae': 'Pospiviroidae', 'chrysanthemum stunt viroid': 'Pospiviroidae', 'potato spindle tuber viroid': 'Pospiviroidae',
    'citrus exocortis viroid': 'Pospiviroidae', 'hop stunt viroid': 'Pospiviroidae', 'tomato chlorotic dwarf viroid': 'Pospiviroidae',
    'columnea latent viroid': 'Pospiviroidae', 'mexican papita viroid': 'Pospiviroidae',
    'avsunviroidae': 'Avsunviroidae', 'peach latent mosaic viroid': 'Avsunviroidae', 'chrysanthemum chlorotic mottle viroid': 'Avsunviroidae',
    'avsunviroid': 'Avsunviroidae',
    'virgaviridae': 'Virgaviridae', 'tobamovirus': 'Virgaviridae', 'tobacco mosaic': 'Virgaviridae',
    'cucumber green mottle mosaic': 'Virgaviridae', 'pepper mild mottle': 'Virgaviridae', 'barley stripe mosaic': 'Virgaviridae',
    'sunn-hemp mosaic': 'Virgaviridae',
    'secoviridae': 'Secoviridae', 'nepovirus': 'Secoviridae', 'comovirus': 'Secoviridae', 'fabavirus': 'Secoviridae',
    'tomato ringspot': 'Secoviridae', 'cherry leaf roll': 'Secoviridae',
    'tombusviridae': 'Tombusviridae', 'carmovirus': 'Tombusviridae', 'pelargonium flower break': 'Tombusviridae',
    'cucumber necrosis': 'Tombusviridae',
    'closteroviridae': 'Closteroviridae', 'crinivirus': 'Closteroviridae', 'ampelovirus': 'Closteroviridae',
    'beet yellows': 'Closteroviridae', 'tomato chlorosis': 'Closteroviridae', 'infectious chlorosis': 'Closteroviridae',
    'rhabdoviridae': 'Rhabdoviridae', 'cytorhabdovirus': 'Rhabdoviridae', 'nucleorhabdovirus': 'Rhabdoviridae',
    'barley yellow striate mosaic': 'Rhabdoviridae',
    'partitiviridae': 'Partitiviridae', 'partitivirus': 'Partitiviridae',
    'tymoviridae': 'Tymoviridae', 'tymovirus': 'Tymoviridae', 'turnip yellow mosaic': 'Tymoviridae',
    'nanoviridae': 'Nanoviridae', 'nanovirus': 'Nanoviridae', 'faba bean necrotic stunt': 'Nanoviridae',
    'caulimoviridae': 'Caulimoviridae', 'caulimovirus': 'Caulimoviridae', 'soymovirus': 'Caulimoviridae',
    'petunia vein clearing': 'Caulimoviridae', 'dahlia mosaic': 'Caulimoviridae',
    'phenuiviridae': 'Phenuiviridae', 'tenuivirus': 'Phenuiviridae', 'rice stripe': 'Phenuiviridae',
    'maize stripe': 'Phenuiviridae',
    'amalgaviridae': 'Amalgaviridae', 'chrysoviridae': 'Chrysoviridae', 'endornaviridae': 'Endornaviridae',
    'hypoviridae': 'Hypoviridae', 'mitoviridae': 'Mitoviridae', 'narnaviridae': 'Narnaviridae'
}

virus_col = None
for col in df.columns:
    if any(k in str(col).lower() for k in ['organism', 'virus', 'species', 'pathogen', 'name']):
        virus_col = col
        break

if virus_col is None:
    virus_col = df.select_dtypes(include=['object']).columns[0]

def assign_family(val):
    if pd.isna(val):
        return 'Unclassified'
    val_lower = str(val).lower()
    for keyword, family_name in family_mapping.items():
        if keyword in val_lower:
            return family_name
    return 'Other Families'

df['Clean_Family'] = df[virus_col].apply(assign_family)
family_counts = df['Clean_Family'].value_counts()
family_counts = family_counts[family_counts.index != 'Unclassified']

# 상위 과들을 적극적으로 분리하고 Other Families 비중 조절
if 'Other Families' in family_counts:
    other_val = family_counts['Other Families']
    family_counts = family_counts.drop('Other Families')
    family_counts = family_counts.head(16)
    family_counts['Other Families'] = other_val
else:
    family_counts = family_counts.head(17)

family_counts = family_counts.reset_index()
family_counts.columns = ['Family', 'Count']
total_records = len(df)

# 4. Generate Figure 2 Treemap with intelligent text scaling and overlap protection
fig, ax = plt.subplots(figsize=(16, 10), dpi=300)

colors = plt.cm.tab20(range(len(family_counts)))

total_count_sum = family_counts['Count'].sum()
labels = []
for index, row in family_counts.iterrows():
    fraction = row['Count'] / total_count_sum
    if fraction < 0.005:  # 아주 작은 블록은 숫자만 표기하여 겹침 방지
        labels.append(f"({row['Count']})")
    elif fraction < 0.015: # 소형 블록은 축약형 이름과 숫자
        short_name = row['Family'].replace('viridae', 'v.').replace('viroidae', 'vrd')
        labels.append(f"{short_name}\n({row['Count']})")
    else:
        labels.append(f"{row['Family']}\n({row['Count']})")

squarify.plot(
    sizes=family_counts['Count'], 
    label=labels, 
    color=colors, 
    alpha=0.88, 
    ax=ax, 
    text_kwargs={'fontsize': 9, 'weight': 'bold', 'color': '#222222'}
)

ax.set_title(f"Figure 2. Taxonomic Distribution of Plant Virus & Viroid Families in Korea (Total: {total_records:,})", fontsize=16, fontweight='bold', pad=20)
ax.axis('off')

# Layout optimization and save
plt.tight_layout()
output_image = "Figure2_taxonomy.png"
plt.savefig(output_image, dpi=300, bbox_inches='tight')
plt.close()

print(f"[Success] Figure 2 successfully generated and saved as '{output_image}'.")


# ==========================================
# FILE: generate_additional_review_figures.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)
print("="*60)
print(" 🔬 추가 추천 리뷰 논문용 도판 생성 시작")
print("="*60)

# 1. 전파 방식별 전장 유전체(Complete Genome) 확보 비율 바차트
plt.figure(figsize=(12, 6))
trans_comp = pd.crosstab(df['Vector_Transmission_Mode'], df['Complete_Genome_Status'], normalize='index') * 100
if 'Complete Genome' in trans_comp.columns:
    trans_comp_sorted = trans_comp['Complete Genome'].sort_values(ascending=False)
    sns.barplot(x=trans_comp_sorted.values, y=trans_comp_sorted.index, palette='crest')
    plt.title('Complete Genome Recovery Rate by Vector Transmission Mode', fontsize=13, fontweight='bold')
    plt.xlabel('Complete Genome Ratio (%)', fontsize=11)
    plt.ylabel('Vector Transmission Mode', fontsize=11)
    plt.xlim(0, 105)
    for i, v in enumerate(trans_comp_sorted.values):
        plt.text(v + 1, i, f"{v:.1f}%", va='center', fontsize=9)
    plt.tight_layout()
    plt.savefig('additional_transmission_completeness.png', dpi=300)
    plt.close()
    print(" - 생성 완료: additional_transmission_completeness.png")

# 2. 주요 확정 도시별 상위 바이러스 검출 빈도 히트맵
top_cities = ['South Korea (Andong)', 'South Korea (Jeonju)', 'South Korea (Jinju)', 'South Korea (Chuncheon)', 'South Korea (Wanju)', 'South Korea (Gunsan)']
sub_city_df = df[df['Geographic_Region'].isin(top_cities)]
if not sub_city_df.empty:
    city_virus = pd.crosstab(sub_city_df['Geographic_Region'], sub_city_df['Viral_Species'])
    # 상위 8개 바이러스 종만 추출하여 가독성 확보
    top_viruses = sub_city_df['Viral_Species'].value_counts().head(8).index
    city_virus_filtered = city_virus[top_viruses]
    
    plt.figure(figsize=(14, 6))
    sns.heatmap(city_virus_filtered, annot=True, fmt='d', cmap='YlGnBu', cbar=True)
    plt.title('Regional Hotspot Profiling: Top Viral Species in Major Cities', fontsize=13, fontweight='bold')
    plt.xlabel('Viral Species', fontsize=11)
    plt.ylabel('Geographic Region', fontsize=11)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.yticks(fontsize=9)
    plt.tight_layout()
    plt.savefig('additional_regional_virus_heatmap.png', dpi=300)
    plt.close()
    print(" - 생성 완료: additional_regional_virus_heatmap.png")

print("\n[완료] 모든 추가 분석 도판이 성공적으로 생성되었습니다!")


# ==========================================
# FILE: generate_all_review_figures.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)
print("="*60)
print(" 🔬 13개 메타 항목별 총설 논문 도판 통합 생성 시작")
print("="*60)

# 1. Genome Type 도넛 차트 (항목 4)
plt.figure(figsize=(7, 7))
genome_counts = df['Genome_Type'].value_counts()
plt.pie(genome_counts, labels=genome_counts.index, autopct='%1.1f%%', startangle=140, colors=sns.color_palette('Set2'))
plt.title('Distribution of Viral Genome Types', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('fig4_genome_type_donut.png', dpi=300)
plt.close()
print(" - [생성완료] fig4_genome_type_donut.png")

# 2. Genome Segmentation 바차트 (항목 5)
plt.figure(figsize=(8, 5))
seg_counts = df['Genome_Segmentation'].value_counts()
sns.barplot(x=seg_counts.index, y=seg_counts.values, palette='Blues_d', edgecolor='black')
plt.title('Genome Segmentation Structure Ratio', fontsize=13, fontweight='bold')
plt.xlabel('Segmentation Type', fontsize=11)
plt.ylabel('Count', fontsize=11)
plt.tight_layout()
plt.savefig('fig5_segmentation_ratio.png', dpi=300)
plt.close()
print(" - [생성완료] fig5_segmentation_ratio.png")

# 3. 진단 기술 코호트 분석 바차트 (항목 11)
cohort_summary = pd.DataFrame({
    'Diagnostic_Cohort': ['Traditional Era (~2014)', 'NGS Metagenomics Era (2015~)'],
    'Success_Rate_%': [31.8, 51.0]
})
plt.figure(figsize=(7, 6))
ax = sns.barplot(
    data=cohort_summary, x='Diagnostic_Cohort', y='Success_Rate_%', 
    palette=['#4682B4', '#2E8B57'], edgecolor='black', width=0.45
)
plt.title('Cohort Analysis: Complete Genome Recovery Rate', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Diagnostic Technology Cohort', fontsize=11)
plt.ylabel('Complete Genome Recovery Rate (%)', fontsize=11)
plt.ylim(0, 105)
plt.xlim(-0.6, 1.6)
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(f"{height:.1f}%", (p.get_x() + p.get_width() / 2., height),
                    ha='center', va='bottom', fontsize=11, fontweight='bold', xytext=(0, 3), textcoords='offset points')
plt.tight_layout()
plt.savefig('diagnostic_cohort_analysis.png', dpi=300)
plt.close()
print(" - [생성완료] diagnostic_cohort_analysis.png")

print("\n[완료] 모든 핵심 총설 도판이 성공적으로 구축되었습니다!")


# ==========================================
# FILE: generate_bipartite_network.py
# ==========================================

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🌐 숙주 작물-바이러스 패밀리 이원 네트워크 분석 시작")
print("="*60)

# 1. Taxonomy Lineage에서 바이러스 패밀리(Family) 추출 함수
def extract_family(lineage):
    if pd.isna(lineage):
        return 'Unclassified'
    parts = [p.strip() for p in str(lineage).split(';')]
    for p in parts:
        if p.endswith('viridae') or p.endswith('viroidae') or 'viroid' in p.lower():
            return p
    return 'Other'

df['Viral_Family'] = df['Taxonomy_Lineage'].apply(extract_family)

# 2. 상위 주요 기주 작물 및 바이러스 패밀리 추출 (네트워크 가독성 확보)
top_hosts = df['Host_Plant'].value_counts().head(6).index
top_families = df['Viral_Family'].value_counts().head(6).index

sub_df = df[df['Host_Plant'].isin(top_hosts) & df['Viral_Family'].isin(top_families)]
edge_counts = sub_df.groupby(['Host_Plant', 'Viral_Family']).size().reset_index(name='weight')

# 3. 네트워크 그래프 구축 (NetworkX Bipartite Graph)
B = nx.Graph()
hosts = edge_counts['Host_Plant'].unique()
families = edge_counts['Viral_Family'].unique()

B.add_nodes_from(hosts, bipartite=0)
B.add_nodes_from(families, bipartite=1)

for _, row in edge_counts.iterrows():
    B.add_edge(row['Host_Plant'], row['Viral_Family'], weight=row['weight'])

# 4. 시각화 레이아웃 설정 (이원 그래프 배치)
plt.figure(figsize=(12, 8))
pos = {}
pos.update((node, (1, index * 2)) for index, node in enumerate(hosts))
pos.update((node, (2, index * 1.5)) for index, node in enumerate(families))

# 노드 색상 지정 (작물: 파란 계열, 바이러스 패밀리: 주황 계열)
nx.draw_networkx_nodes(B, pos, nodelist=hosts, node_color='#2b5c8f', node_size=900, alpha=0.9)
nx.draw_networkx_nodes(B, pos, nodelist=families, node_color='#d95f02', node_size=900, alpha=0.9)

# 엣지 두께는 감염 빈도(Weight)에 비례하도록 설정
weights = [B[u][v]['weight'] * 0.08 for u, v in B.edges()]
nx.draw_networkx_edges(B, pos, width=weights, alpha=0.5, edge_color='gray')
nx.draw_networkx_labels(B, pos, font_size=10, font_weight='bold')

plt.title('Bipartite Network: Host Plants vs. Major Viral Families (Hub Analysis)', fontsize=14, fontweight='bold', pad=20)
plt.axis('off')
plt.tight_layout()
plt.savefig('host_virus_bipartite_network.png', dpi=300)
plt.close()

print("[완료] 이원 네트워크 그래프 저장 성공: host_virus_bipartite_network.png")


# ==========================================
# FILE: generate_bipartite_network_v2.py
# ==========================================

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🌐 개선된 숙주-바이러스 종(Species) 네트워크 분석 시작")
print("="*60)

# 1. 바이러스 종(Species) 명칭 정제 함수
def clean_virus_name(name):
    if pd.isna(name):
        return 'Unknown'
    return str(name).strip()

df['Clean_Species'] = df['Viral_Species'].apply(clean_virus_name)

# 2. 상위 주요 기주 작물 및 바이러스 종 추출 (가독성 최적화)
top_hosts = df['Host_Plant'].value_counts().head(6).index
top_species = df['Clean_Species'].value_counts().head(8).index

sub_df = df[df['Host_Plant'].isin(top_hosts) & df['Clean_Species'].isin(top_species)]
edge_counts = sub_df.groupby(['Host_Plant', 'Clean_Species']).size().reset_index(name='weight')

# 3. 네트워크 그래프 구축 (Bipartite Graph)
B = nx.Graph()
hosts = edge_counts['Host_Plant'].unique()
species = edge_counts['Clean_Species'].unique()

B.add_nodes_from(hosts, bipartite=0)
B.add_nodes_from(species, bipartite=1)

for _, row in edge_counts.iterrows():
    B.add_edge(row['Host_Plant'], row['Clean_Species'], weight=row['weight'])

# 4. 시각화 캔버스 확대 및 여백(Margin) 확보로 글자 잘림 방지
plt.figure(figsize=(16, 9))
pos = {}
pos.update((node, (0.2, index * 2.2)) for index, node in enumerate(hosts))
pos.update((node, (2.8, index * 1.6)) for index, node in enumerate(species))

# 다양한 파스텔 및 비비드 컬러 팔레트 적용
host_colors = plt.cm.tab10(range(len(hosts)))
species_colors = plt.cm.Set2(range(len(species)))

node_color_map = {}
for i, h in enumerate(hosts):
    node_color_map[h] = host_colors[i]
for i, s in enumerate(species):
    node_color_map[s] = species_colors[i]

# 노드 드로잉 (외곽선 테두리 및 음영 효과 추가)
nx.draw_networkx_nodes(B, pos, nodelist=hosts, node_color=[node_color_map[h] for h in hosts], node_size=1200, alpha=0.9, edgecolors='black', linewidths=1.5)
nx.draw_networkx_nodes(B, pos, nodelist=species, node_color=[node_color_map[s] for s in species], node_size=1200, alpha=0.9, edgecolors='black', linewidths=1.5)

# 엣지 두께 및 음영(Alpha) 조절
weights = [B[u][v]['weight'] * 0.1 for u, v in B.edges()]
nx.draw_networkx_edges(B, pos, width=weights, alpha=0.45, edge_color='#444444')

# 텍스트 레이블 출력 (좌우 여백 확보로 잘림 현상 원천 차단)
nx.draw_networkx_labels(B, pos, font_size=10, font_weight='bold', font_family='sans-serif')

plt.title('Bipartite Network: Host Plants vs. Major Viral/Viroid Species', fontsize=15, fontweight='bold', pad=25)
plt.xlim(-0.5, 3.5)
plt.axis('off')
plt.tight_layout()
plt.savefig('host_virus_bipartite_network_v2.png', dpi=300)
plt.close()

print("[완료] 개선된 이원 네트워크 그래프 저장 성공: host_virus_bipartite_network_v2.png")


# ==========================================
# FILE: generate_clean_segmentation.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 📊 정제된 유전체 분절 구조 비율 바차트 생성 시작")
print("="*60)

# 데이터 정제 및 주요 핵심 카테고리 그룹화
def clean_segmentation(val):
    if pd.isna(val):
        return 'Unspecified / Other'
    val_str = str(val).strip()
    if 'Monopartite' in val_str:
        return 'Monopartite'
    elif 'Multipartite' in val_str or 'Bipartite' in val_str:
        return 'Multipartite'
    else:
        return 'Unspecified / Other'

df['Clean_Segmentation'] = df['Genome_Segmentation'].apply(clean_segmentation)
seg_counts = df['Clean_Segmentation'].value_counts()

# 시각화 설정
plt.figure(figsize=(8, 6))
ax = sns.barplot(
    x=seg_counts.index, y=seg_counts.values, 
    palette=['#4682B4', '#5F9EA0', '#708090'], edgecolor='black', width=0.5
)

plt.title('Genome Segmentation Structure Ratio (Cleaned)', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Segmentation Type', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.ylim(0, seg_counts.max() * 1.15)

# 바 상단 데이터 개수 레이블 부착
for p in ax.patches:
    height = p.get_height()
    if height > 0:
        ax.annotate(f"{int(height):,}",
                    (p.get_x() + p.get_width() / 2., height),
                    ha='center', va='bottom', fontsize=11, fontweight='bold',
                    xytext=(0, 4), textcoords='offset points')

plt.tight_layout()
output_fig = "fig5_segmentation_ratio_clean.png"
plt.savefig(output_fig, dpi=300)
plt.close()

print(f"[완료] 정제된 분절 구조 그래프 저장 성공: {output_fig}")


# ==========================================
# FILE: generate_figure3.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 6), dpi=300)

# [Figure 3-A] 완전 유전체 복구 현황
ax1 = plt.subplot(1, 2, 1)
status_counts = df['Complete_Genome_Status'].value_counts().reset_index()
status_counts.columns = ['Status', 'Count']
sns.barplot(data=status_counts, x='Count', y='Status', hue='Status', palette='Blues_r', ax=ax1)
if ax1.legend_: ax1.legend_.remove()
ax1.set_title('(A) Complete Genome Recovery Status', fontweight='bold', fontsize=12)
ax1.set_xlabel('Number of Records', fontweight='bold')
ax1.set_ylabel('Genome Status', fontweight='bold')

# [Figure 3-B] 주요 병원성 증상 분포 Top 6
ax2 = plt.subplot(1, 2, 2)
symptom_counts = df['Pathogenicity_Symptoms'].value_counts().head(6).reset_index()
symptom_counts.columns = ['Symptoms', 'Count']
sns.barplot(data=symptom_counts, x='Count', y='Symptoms', hue='Symptoms', palette='Reds_r', ax=ax2)
if ax2.legend_: ax2.legend_.remove()
ax2.set_title('(B) Top Pathogenicity Symptoms', fontweight='bold', fontsize=12)
ax2.set_xlabel('Number of Records', fontweight='bold')
ax2.set_ylabel('Pathogenicity Symptoms', fontweight='bold')

plt.tight_layout()
plt.savefig('plant_viruses_figure3.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 세 번째 논문 도판이 'plant_viruses_figure3.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_figure4.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 6), dpi=300)

# [Figure 4-A] 지역별 보고 건수 Top 지역 (Unspecified 제외 상위 지역)
ax1 = plt.subplot(1, 2, 1)
region_counts = df[~df['Geographic_Region'].str.contains('Unspecified')]['Geographic_Region'].value_counts().head(5).reset_index()
region_counts.columns = ['Region', 'Count']
sns.barplot(data=region_counts, x='Count', y='Region', hue='Region', palette='Purples_r', ax=ax1)
if ax1.legend_: ax1.legend_.remove()
ax1.set_title('(A) Top Geographic Regions (Excl. Unspecified)', fontweight='bold', fontsize=12)
ax1.set_xlabel('Number of Records', fontweight='bold')
ax1.set_ylabel('Geographic Region', fontweight='bold')

# [Figure 4-B] 숙주 식물과 바이러스 종 조합 Top 5
ax2 = plt.subplot(1, 2, 2)
top_combos = (df['Host_Plant'] + " & " + df['Viral_Species']).value_counts().head(5).reset_index()
top_combos.columns = ['Host_Virus', 'Count']
sns.barplot(data=top_combos, x='Count', y='Host_Virus', hue='Host_Virus', palette='crest', ax=ax2)
if ax2.legend_: ax2.legend_.remove()
ax2.set_title('(B) Top Host-Virus Associations', fontweight='bold', fontsize=12)
ax2.set_xlabel('Number of Records', fontweight='bold')
ax2.set_ylabel('Host & Viral Species Pair', fontweight='bold')

plt.tight_layout()
plt.savefig('plant_viruses_figure4.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 네 번째 논문 도판이 'plant_viruses_figure4.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_figure5.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 6), dpi=300)

# [Figure 5-A] 연도별 유전체 유형 구성 비율 변화 (Crosstab Stacked Bar or Line)
ax1 = plt.subplot(1, 2, 1)
year_genome = pd.crosstab(df['Collection_Year'], df['Genome_Type'])
year_genome.plot(kind='area', stacked=True, colormap='viridis', alpha=0.85, ax=ax1)
ax1.set_title('(A) Temporal Dynamics of Genome Types', fontweight='bold', fontsize=12)
ax1.set_xlabel('Collection Year', fontweight='bold')
ax1.set_ylabel('Number of Records', fontweight='bold')
ax1.legend(title='Genome Type', loc='upper left', fontsize=8)

# [Figure 5-B] 연도별 유전체 완성도 상태 추이 (Complete vs Partial)
ax2 = plt.subplot(1, 2, 2)
year_status = pd.crosstab(df['Collection_Year'], df['Complete_Genome_Status'])
year_status.plot(kind='line', linewidth=2, ax=ax2)
ax2.set_title('(B) Annual Trend of Genome Recovery Status', fontweight='bold', fontsize=12)
ax2.set_xlabel('Collection Year', fontweight='bold')
ax2.set_ylabel('Number of Records', fontweight='bold')
ax2.legend(title='Genome Status', loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig('plant_viruses_figure5.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 다섯 번째 논문 도판이 'plant_viruses_figure5.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_figure6.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 6), dpi=300)

# [Figure 6-A] 유전체 분절 구조(Clean Segmentation)와 주요 벡터 전파 경로 교차 분석
ax1 = plt.subplot(1, 2, 1)
def clean_segmentation(val):
    if pd.isna(val): return 'Other'
    val_str = str(val).strip()
    if 'Monopartite' in val_str: return 'Monopartite'
    elif 'Multipartite' in val_str or 'Bipartite' in val_str: return 'Multipartite'
    else: return 'Other/Viroid'

df['Clean_Seg'] = df['Genome_Segmentation'].apply(clean_segmentation)
seg_vector = pd.crosstab(df['Clean_Seg'], df['Vector_Transmission_Mode']).iloc[:, :4] # 상위 4개 매개체만 발췌
seg_vector.plot(kind='bar', stacked=True, colormap='Set2', ax=ax1)
ax1.set_title('(A) Vector Modes by Genome Segmentation', fontweight='bold', fontsize=12)
ax1.set_xlabel('Genome Segmentation', fontweight='bold')
ax1.set_ylabel('Number of Records', fontweight='bold')
ax1.legend(title='Transmission Mode', loc='upper right', fontsize=8)
plt.xticks(rotation=0)

# [Figure 6-B] 유전체 유형별 서열 길이 바이올린 플롯 (로그 스케일 또는 1만 이하 분포)
ax2 = plt.subplot(1, 2, 2)
sns.violinplot(data=df[df['Sequence_Length'] <= 15000], x='Genome_Type', y='Sequence_Length', palette='Pastel1', ax=ax2, inner='quartile')
ax2.set_title('(B) Length Distribution by Genome Type (<=15kb)', fontweight='bold', fontsize=12)
ax2.set_xlabel('Genome Type', fontweight='bold')
ax2.set_ylabel('Sequence Length (bp)', fontweight='bold')
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig('plant_viruses_figure6.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 여섯 번째 논문 도판이 'plant_viruses_figure6.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_figure7.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 6), dpi=300)

# [Figure 7-A] 벡터 전파 경로별 유전체 완성도 상태 (Complete vs Partial)
ax1 = plt.subplot(1, 2, 1)
# 보기 쉽도록 주요 6개 매개체만 필터링
top_vectors = df['Vector_Transmission_Mode'].value_counts().head(6).index
df_filtered = df[df['Vector_Transmission_Mode'].isin(top_vectors)]
vec_status = pd.crosstab(df_filtered['Vector_Transmission_Mode'], df_filtered['Complete_Genome_Status'])
vec_status.plot(kind='bar', stacked=True, colormap='Blues', ax=ax1)
ax1.set_title('(A) Genome Status by Vector Transmission Mode', fontweight='bold', fontsize=12)
ax1.set_xlabel('Vector Transmission Mode', fontweight='bold')
ax1.set_ylabel('Number of Records', fontweight='bold')
ax1.legend(title='Genome Status', loc='upper right', fontsize=8)
plt.xticks(rotation=30, ha='right')

# [Figure 7-B] 상위 숙주 작물과 유전체 유형 간의 히트맵 (Host vs Genome Type)
ax2 = plt.subplot(1, 2, 2)
top_hosts = df['Host_Plant'].value_counts().head(6).index
df_hosts = df[df['Host_Plant'].isin(top_hosts)]
host_genome = pd.crosstab(df_hosts['Host_Plant'], df_hosts['Genome_Type'])
sns.heatmap(host_genome, annot=True, fmt='d', cmap='YlGnBu', cbar=True, ax=ax2)
ax2.set_title('(B) Heatmap of Host Plants vs. Genome Types', fontweight='bold', fontsize=12)
ax2.set_xlabel('Genome Type', fontweight='bold')
ax2.set_ylabel('Host Plant', fontweight='bold')
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig('plant_viruses_figure7.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 일곱 번째 논문 도판이 'plant_viruses_figure7.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_figures.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

fig = plt.figure(figsize=(14, 11), dpi=300)

ax1 = plt.subplot(2, 2, 1)
genome_counts = df['Genome_Type'].value_counts().reset_index()
genome_counts.columns = ['Genome_Type', 'Count']
sns.barplot(data=genome_counts, x='Count', y='Genome_Type', hue='Genome_Type', palette='Blues_r', legend=False, ax=ax1)
ax1.set_title('(A) Distribution of Viral Genome Types', fontweight='bold', fontsize=12)
ax1.set_xlabel('Number of Records', fontweight='bold')
ax1.set_ylabel('Genome Type', fontweight='bold')

ax2 = plt.subplot(2, 2, 2)
sns.boxplot(data=df, x='Sequence_Length', y='Genome_Type', hue='Genome_Type', palette='Set2', legend=False, ax=ax2, showfliers=False)
ax2.set_title('(B) Sequence Length by Genome Type', fontweight='bold', fontsize=12)
ax2.set_xlabel('Sequence Length (bp)', fontweight='bold')
ax2.set_ylabel('')

ax3 = plt.subplot(2, 2, 3)
top_hosts = df['Host_Plant'].value_counts().head(6).reset_index()
top_hosts.columns = ['Host_Plant', 'Count']
sns.barplot(data=top_hosts, x='Count', y='Host_Plant', hue='Host_Plant', palette='Purples_r', legend=False, ax=ax3)
ax3.set_title('(C) Top 6 Host Plants', fontweight='bold', fontsize=12)
ax3.set_xlabel('Number of Records', fontweight='bold')
ax3.set_ylabel('Host Plant', fontweight='bold')

ax4 = plt.subplot(2, 2, 4)
top_viruses = df['Viral_Species'].value_counts().head(6).reset_index()
top_viruses.columns = ['Viral_Species', 'Count']
sns.barplot(data=top_viruses, x='Count', y='Viral_Species', hue='Viral_Species', palette='Greens_r', legend=False, ax=ax4)
ax4.set_title('(D) Top 6 Viral/Viroid Species', fontweight='bold', fontsize=12)
ax4.set_xlabel('Number of Records', fontweight='bold')
ax4.set_ylabel('Viral Species', fontweight='bold')

plt.tight_layout()
plt.savefig('plant_viruses_multipanel_figure.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 논문용 고해상도 피겨가 'plant_viruses_multipanel_figure.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_host_completeness.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

top_hosts = df['Host_Plant'].value_counts().head(8).index
df_hosts = df[df['Host_Plant'].isin(top_hosts)].copy()

df_hosts['Is_Complete'] = df_hosts['Complete_Genome_Status'].apply(lambda x: 1 if 'Complete Genome' in str(x) else 0)
completeness_by_host = df_hosts.groupby('Host_Plant')['Is_Complete'].mean() * 100
completeness_by_host = completeness_by_host.sort_values(ascending=False).reset_index()
completeness_by_host.columns = ['Host_Plant', 'Complete_Rate']

fig, ax = plt.subplots(figsize=(11, 6), dpi=300)
sns.barplot(data=completeness_by_host, x='Complete_Rate', y='Host_Plant', hue='Host_Plant', palette='crest', legend=False, ax=ax)

for p in ax.patches:
    width = p.get_width()
    if width > 0:
        ax.annotate(f'{width:.1f}%',
                    xy=(width, p.get_y() + p.get_height() / 2),
                    xytext=(5, 0), textcoords='offset points',
                    ha='left', va='center', fontweight='bold', fontsize=9)

ax.set_title('Complete Genome Recovery Rate by Major Host Plant', fontweight='bold', fontsize=13, pad=15)
ax.set_xlabel('Complete Genome Recovery Rate (%)', fontweight='bold', fontsize=11)
ax.set_ylabel('Host Plant', fontweight='bold', fontsize=11)
ax.set_xlim(0, 115)

plt.tight_layout()
plt.savefig('fig_host_completeness_rate.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 숙주별 유전체 복구율 도판이 'fig_host_completeness_rate.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_host_pathology.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

df = pd.read_csv('plant_viruses_metadata_enriched.csv')

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

# 상위 5개 주요 숙주 작물 및 주요 증상 필터링
top_hosts = df['Host_Plant'].value_counts().head(5).index
df_filtered = df[df['Host_Plant'].isin(top_hosts)]

# 크로스탭 생성 및 백분율 변환 (100% Stacked Bar)
ct = pd.crosstab(df_filtered['Host_Plant'], df_filtered['Pathogenicity_Symptoms'], normalize='index') * 100

fig, ax = plt.subplots(figsize=(12, 7), dpi=300)
ct.plot(kind='bar', stacked=True, colormap='Spectral', edgecolor='black', linewidth=0.5, ax=ax)

ax.set_title('Proportion of Pathogenicity Symptoms across Major Host Plants', fontweight='bold', fontsize=13, pad=15)
ax.set_xlabel('Host Plant', fontweight='bold', fontsize=11)
ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=11)
ax.legend(title='Pathogenicity Symptoms', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
plt.xticks(rotation=25, ha='right')

plt.tight_layout()
plt.savefig('fig_host_pathology_stacked.png', dpi=300, bbox_inches='tight')
plt.close()

print("[완료] 숙주별 병원성 증상 누적 비율 도판이 'fig_host_pathology_stacked.png'로 성공적으로 저장되었습니다.")


# ==========================================
# FILE: generate_korea_panels.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading Korean virus metadata for integrated panel generation...")
df = pd.read_csv(csv_file)

# 1. 기주 학명 정제 로직
synonym_map = {
    "Soybean": "Glycine max",
    "Garlic": "Allium sativum",
    "Pepper": "Capsicum annuum",
    "Hot pepper": "Capsicum annuum",
    "Tomato": "Solanum lycopersicum",
    "Peach": "Prunus persica",
    "Radish": "Raphanus sativus"
}

def strict_scientific_host(host_str):
    if pd.isna(host_str):
        return None
    clean_str = str(host_str).strip()
    clean_lower = clean_str.lower()
    for key, val in synonym_map.items():
        if key.lower() in clean_lower:
            return val
    exclude_keywords = ["apis", "honeybee", "insect", "human", "virus", "unknown", "none", "n/a"]
    if any(ex in clean_lower for ex in exclude_keywords):
        return None
    parts = clean_str.split()
    if len(parts) >= 2:
        genus = parts[0].capitalize()
        species = parts[1].lower()
        genus = re.sub(r'[^a-zA-Z]', '', genus)
        species = re.sub(r'[^a-zA-Z]', '', species)
        if genus and species:
            return f"{genus} {species}"
    return None

df['Strict_Host'] = df['Host'].apply(strict_scientific_host)

# 연도 추출 로직
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        y = int(match.group(1))
        if 1980 <= y <= 2026:
            return y
    return None

df['Year'] = df['Collection_Date'].apply(extract_year)

# --- 통합 패널 그림 생성 (2x2 Grid) ---
fig, axes = plt.subplots(2, 2, figsize=(18, 14))

# [Panel A] Top Virus Species
top_viruses = df['Organism'].value_counts().head(8).reset_index()
top_viruses.columns = ['Virus', 'Count']
sns.barplot(data=top_viruses, x='Count', y='Virus', ax=axes[0, 0], palette='Blues_r')
axes[0, 0].set_title('(A) Top Plant Viruses & Viroids in Korea', fontsize=13, fontweight='bold')
axes[0, 0].set_xlabel('Number of Sequences', fontsize=11)
axes[0, 0].set_ylabel('Virus Organism', fontsize=11)

# [Panel B] Top Host Scientific Names
valid_hosts = df.dropna(subset=['Strict_Host'])
top_hosts = valid_hosts['Strict_Host'].value_counts().head(8).reset_index()
top_hosts.columns = ['Host', 'Count']
sns.barplot(data=top_hosts, x='Count', y='Host', ax=axes[0, 1], palette='Greens_r')
axes[0, 1].set_title('(B) Top Plant Host Scientific Names', fontsize=13, fontweight='bold')
axes[0, 1].set_xlabel('Number of Sequences', fontsize=11)
axes[0, 1].set_ylabel('Host Scientific Name', fontsize=11)

# [Panel C] Virus-Host Interaction Heatmap (Top 6 x Top 6)
top_v6 = df['Organism'].value_counts().head(6).index.tolist()
top_h6 = valid_hosts['Strict_Host'].value_counts().head(6).index.tolist()
df_sub = df[df['Organism'].isin(top_v6) & df['Strict_Host'].isin(top_h6)]
matrix = pd.crosstab(df_sub['Organism'], df_sub['Strict_Host'])
sns.heatmap(matrix, annot=True, fmt="d", cmap="YlGnBu", ax=axes[1, 0], cbar=False)
axes[1, 0].set_title('(C) Virus-Host Interaction Matrix', fontsize=13, fontweight='bold')
axes[1, 0].set_xlabel('Host Scientific Name', fontsize=11)
axes[1, 0].set_ylabel('Virus Organism', fontsize=11)
axes[1, 0].tick_params(axis='x', rotation=30)

# [Panel D] Yearly Trend (1990-2026)
yearly = df[df['Year'] >= 1990]['Year'].value_counts().sort_index().reset_index()
yearly.columns = ['Year', 'Count']
axes[1, 1].plot(yearly['Year'], yearly['Count'], marker='o', color='#1f77b4', linewidth=2.5)
axes[1, 1].set_title('(D) Temporal Trend of Submissions (1990–2026)', fontsize=13, fontweight='bold')
axes[1, 1].set_xlabel('Year', fontsize=11)
axes[1, 1].set_ylabel('Number of Sequences', fontsize=11)
axes[1, 1].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
panel_path = os.path.join(fig_dir, "korea_integrated_panels.png")
plt.savefig(panel_path, dpi=300)
plt.close()

print(f"Successfully generated integrated Korea panel figure: {panel_path}")


# ==========================================
# FILE: generate_paper_panels.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)
print("="*60)
print(" 🔬 논문용 핵심 시각화 패널 일괄 생성 시작")
print("="*60)

# 1. 유전체 유형 및 분절 구조별 서열 길이 박스플롯 (이상치 제외)
plt.figure(figsize=(10, 6))
sns.boxplot(x='Genome_Type', y='Sequence_Length', hue='Genome_Segmentation', data=df, palette='Set2', showfliers=False)
plt.title('(A) Sequence Length Distribution by Genome Architecture', fontsize=13, fontweight='bold')
plt.xlabel('Genome Type', fontsize=11)
plt.ylabel('Sequence Length (bp)', fontsize=11)
plt.xticks(rotation=25)
plt.legend(title='Genome Segmentation', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig('panel_1_genome_length.png', dpi=300)
plt.close()

# 2. 주요 기주 작물별 병원성 증상 히트맵
top_hosts = df['Host_Plant'].value_counts().head(8).index
host_symptom = pd.crosstab(df[df['Host_Plant'].isin(top_hosts)]['Host_Plant'], df['Pathogenicity_Symptoms'])
plt.figure(figsize=(11, 5))
sns.heatmap(host_symptom, annot=True, fmt='d', cmap='Blues', cbar=True)
plt.title('(B) Pathogenicity Symptoms Heatmap by Major Host Plants', fontsize=13, fontweight='bold')
plt.xlabel('Pathogenicity Symptoms', fontsize=11)
plt.ylabel('Host Plant', fontsize=11)
plt.xticks(rotation=35, ha='right', fontsize=9)
plt.tight_layout()
plt.savefig('panel_2_host_symptom.png', dpi=300)
plt.close()

# 3. 매개체 전파 방식과 증상 간의 교차 히트맵
transmission_symptom = pd.crosstab(df['Vector_Transmission_Mode'], df['Pathogenicity_Symptoms'])
plt.figure(figsize=(11, 5))
sns.heatmap(transmission_symptom, annot=True, fmt='d', cmap='YlOrRd', cbar=True)
plt.title('(C) Transmission Mode vs Pathogenicity Symptoms', fontsize=13, fontweight='bold')
plt.xlabel('Pathogenicity Symptoms', fontsize=11)
plt.ylabel('Vector Transmission Mode', fontsize=11)
plt.xticks(rotation=35, ha='right', fontsize=9)
plt.tight_layout()
plt.savefig('panel_3_transmission_symptom.png', dpi=300)
plt.close()

# 4. 정제된 지역별 핫스팟 바차트 (상위 10개 주요 도시)
top_regions = df['Geographic_Region'].value_counts().head(10)
plt.figure(figsize=(11, 5))
sns.barplot(x=top_regions.index, y=top_regions.values, hue=top_regions.index, palette='viridis', legend=False)
plt.title('(D) Top 10 Geographic Hotspots of Virus Detections', fontsize=13, fontweight='bold')
plt.xlabel('Geographic Region', fontsize=11)
plt.ylabel('Detection Count', fontsize=11)
plt.xticks(rotation=35, ha='right', fontsize=9)
plt.tight_layout()
plt.savefig('panel_4_geographic_hotspots.png', dpi=300)
plt.close()

# 5. 연도별 신규 바이러스 검출 추이 (시계열 라인 차트)
year_counts = df['Collection_Year'].value_counts().sort_index()
plt.figure(figsize=(11, 5))
sns.lineplot(x=year_counts.index, y=year_counts.values, marker='o', color='b', linewidth=2.5)
plt.title('(E) Temporal Trend of Plant Virus Detections', fontsize=13, fontweight='bold')
plt.xlabel('Collection Year', fontsize=11)
plt.ylabel('Detection Count', fontsize=11)
plt.xticks(rotation=30)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('panel_5_temporal_trend.png', dpi=300)
plt.close()

print("\n[완료] 논문용 핵심 시각화 패널 5종 생성 완료:")
print(" 1. panel_1_genome_length.png (유전체 구조별 서열 길이)")
print(" 2. panel_2_host_symptom.png (기주 작물별 증상 히트맵)")
print(" 3. panel_3_transmission_symptom.png (전파 방식 대 증상 히트맵)")
print(" 4. panel_4_geographic_hotspots.png (지역별 주요 핫스팟)")
print(" 5. panel_5_temporal_trend.png (연도별 검출 추이 트렌드)")


# ==========================================
# FILE: generate_sankey.py
# ==========================================

import pandas as pd
import plotly.graph_objects as go
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🌊 다중 계층 산케이 다이어그램 데이터 가공 및 생성 시작")
print("="*60)

# 1. 바이러스 패밀리(Family) 추출 함수
def extract_family(lineage):
    if pd.isna(lineage):
        return 'Unclassified'
    parts = [p.strip() for p in str(lineage).split(';')]
    for p in parts:
        if p.endswith('viridae') or p.endswith('viroidae') or 'viroid' in p.lower():
            return p
    return 'Other'

df['Viral_Family'] = df['Taxonomy_Lineage'].apply(extract_family)
df['Host_Clean'] = df['Host_Plant'].fillna('Unknown Host')
df['Vector_Clean'] = df['Vector_Transmission_Mode'].fillna('Unknown Vector')

# 2. 상위 주요 카테고리 추출 (흐름의 가독성 확보)
top_hosts = df['Host_Clean'].value_counts().head(5).index
top_vectors = df['Vector_Clean'].value_counts().head(5).index
top_families = df['Viral_Family'].value_counts().head(5).index

sankey_df = df[
    df['Host_Clean'].isin(top_hosts) & 
    df['Vector_Clean'].isin(top_vectors) & 
    df['Viral_Family'].isin(top_families)
].copy()

# 3. 1단계 흐름: Host -> Vector
flow1 = sankey_df.groupby(['Host_Clean', 'Vector_Clean']).size().reset_index(name='weight')
flow1.columns = ['source', 'target', 'weight']

# 4. 2단계 흐름: Vector -> Family
flow2 = sankey_df.groupby(['Vector_Clean', 'Viral_Family']).size().reset_index(name='weight')
flow2.columns = ['source', 'target', 'weight']

all_flows = pd.concat([flow1, flow2], ignore_index=True)

# 5. 노드 인덱싱 및 매핑
labels = list(pd.concat([all_flows['source'], all_flows['target']]).unique())
label_map = {label: i for i, label in enumerate(labels)}

source_indices = all_flows['source'].map(label_map).tolist()
target_indices = all_flows['target'].map(label_map).tolist()
weights = all_flows['weight'].tolist()

# 6. Plotly Sankey 다이어그램 객체 생성
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 18,
      thickness = 22,
      line = dict(color = "black", width = 0.5),
      label = labels,
      color = ["#2b5c8f"] * len(labels)  # 일관된 컬러 테마 적용
    ),
    link = dict(
      source = source_indices,
      target = target_indices,
      value = weights,
      color = "rgba(100, 150, 200, 0.4)"  # 투명도 있는 음영 링크 효과
    ))])

fig.update_layout(
    title_text="Multi-Layer Flow: Host Plant ➔ Vector Transmission Mode ➔ Viral Family", 
    font_size=11, 
    font_family="sans-serif",
    width=1100, 
    height=650
)

# HTML 파일로 저장 (인터랙티브 웹 브라우저 확인 가능)
output_file = "sankey_host_vector_family.html"
fig.write_html(output_file)
print(f"[완료] 산케이 다이어그램 인터랙티브 HTML 저장 성공: {output_file}")


# ==========================================
# FILE: generate_sankey_bolder.py
# ==========================================

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.colors as mcolors
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🌊 글씨 크기 및 굵기가 강화된 산케이 다이어그램 생성 시작")
print("="*60)

def clean_virus_name(name):
    if pd.isna(name):
        return 'Unknown'
    return str(name).strip()

df['Clean_Species'] = df['Viral_Species'].apply(clean_virus_name)
df['Host_Clean'] = df['Host_Plant'].fillna('Unknown Host')
df['Vector_Clean'] = df['Vector_Transmission_Mode'].fillna('Unknown Vector')

top_hosts = df['Host_Clean'].value_counts().head(8).index
top_vectors = df['Vector_Clean'].value_counts().head(8).index
top_species = df['Clean_Species'].value_counts().head(10).index

sankey_df = df[
    df['Host_Clean'].isin(top_hosts) & 
    df['Vector_Clean'].isin(top_vectors) & 
    df['Clean_Species'].isin(top_species)
].copy()

flow1 = sankey_df.groupby(['Host_Clean', 'Vector_Clean']).size().reset_index(name='weight')
flow1.columns = ['source', 'target', 'weight']

flow2 = sankey_df.groupby(['Vector_Clean', 'Clean_Species']).size().reset_index(name='weight')
flow2.columns = ['source', 'target', 'weight']

all_flows = pd.concat([flow1, flow2], ignore_index=True)

labels = list(pd.concat([all_flows['source'], all_flows['target']]).unique())
label_map = {label: i for i, label in enumerate(labels)}

source_indices = all_flows['source'].map(label_map).tolist()
target_indices = all_flows['target'].map(label_map).tolist()
weights = all_flows['weight'].tolist()

color_palette = px.colors.qualitative.Dark24 * 3
node_colors = color_palette[:len(labels)]

def hex_to_rgba(hex_str, alpha=0.45):
    try:
        rgb = mcolors.to_rgb(hex_str)
        return f"rgba({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)}, {alpha})"
    except:
        return f"rgba(100, 150, 200, {alpha})"

link_colors = [hex_to_rgba(color_palette[src % len(color_palette)], alpha=0.4) for src in source_indices]

# Plotly Sankey 설정 (font_size를 16으로 키우고 HTML <b> 태그를 활용해 굵게 표현)
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 14,
      thickness = 22,
      line = dict(color = "black", width = 1.0),
      label = [f"<b>{label}</b>" for label in labels],  # 각 레이블에 볼드체 적용
      color = node_colors
    ),
    link = dict(
      source = source_indices,
      target = target_indices,
      value = weights,
      color = link_colors
    ))])

fig.update_layout(
    title_text="<b>Multi-Layer Flow: Host Plant ➔ Vector Transmission Mode ➔ Viral/Viroid Species</b>", 
    font_size=16,               # 글자 크기를 16으로 대폭 확대
    font_family="Arial, sans-serif",
    font_color="black",
    width=1400, 
    height=900
)

output_file = "sankey_bolder.html"
fig.write_html(output_file)
print(f"[완료] 가독성이 강화된 산케이 다이어그램 저장 성공: {output_file}")


# ==========================================
# FILE: generate_sankey_diagram.py
# ==========================================

import pandas as pd
import os
import re

try:
    import plotly.graph_objects as go
except ImportError:
    import subprocess
    subprocess.check_call(["pip", "install", "plotly"])
    import plotly.graph_objects as go

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_plant_viruses_metadata.csv")
print("Loading metadata for Sankey flow diagram...")
df = pd.read_csv(csv_file)

# 1. 연도 추출
def extract_year(date_str):
    if pd.isna(date_str):
        return None
    match = re.search(r'(20\d{2}|19\d{2})', str(date_str))
    if match:
        y = int(match.group(1))
        if 2010 <= y <= 2026: # 최근 주요 흐름 집중 (2010년 이후)
            return str(y)
    return "Before 2010"

df['Year_Group'] = df['Collection_Date'].apply(extract_year)
df_clean = df.dropna(subset=['Organism', 'Host']).copy()

# 데이터 정제 (상위 주요 바이러스 및 주요 기주 중심로 압축하여 흐름 가독성 확보)
top_viruses = df_clean['Organism'].value_counts().head(5).index.tolist()
def filter_virus(org):
    return org if org in top_viruses else "Other Viruses"

df_clean['Virus_Top'] = df_clean['Organism'].apply(filter_virus)

# 2. 생키 다이어그램 노드 및 링크 데이터 준비
# 흐름: Year_Group -> Virus_Top -> Host
# 여기서는 Year -> Virus 연결 집계
link_df = df_clean.groupby(['Year_Group', 'Virus_Top']).size().reset_index(name='Count')

# 노드 리스트 생성
all_nodes = list(link_df['Year_Group'].unique()) + list(link_df['Virus_Top'].unique())
node_indices = {node: idx for idx, node in enumerate(all_nodes)}

source = link_df['Year_Group'].map(node_indices)
target = link_df['Virus_Top'].map(node_indices)
value = link_df['Count']

# Plotly 생키 다이어그램 생성
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
      label = all_nodes,
      color = "blue"
    ),
    link = dict(
      source = source,
      target = target,
      value = value
    ))])

fig.update_layout(title_text="Flow Diagram: Year to Plant Virus Distribution in Korea", font_size=12)

# HTML 파일로 저장 (인터랙티브하게 브라우저에서 열어볼 수 있음)
html_path = os.path.join(fig_dir, "korea_virus_sankey_flow.html")
fig.write_html(html_path)
print(f"-> Saved interactive Sankey diagram HTML: {html_path}")


# ==========================================
# FILE: generate_sankey_rich.py
# ==========================================

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.colors as mcolors
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🌊 풍부한 데이터 기반 멀티컬러 산케이 다이어그램 생성 시작")
print("="*60)

def clean_virus_name(name):
    if pd.isna(name):
        return 'Unknown'
    return str(name).strip()

df['Clean_Species'] = df['Viral_Species'].apply(clean_virus_name)
df['Host_Clean'] = df['Host_Plant'].fillna('Unknown Host')
df['Vector_Clean'] = df['Vector_Transmission_Mode'].fillna('Unknown Vector')

# 카테고리 상위 개수를 늘려 더욱 풍부하고 다양한 관계망 표현 (작물 8개, 매개체 8개, 바이러스 종 10개)
top_hosts = df['Host_Clean'].value_counts().head(8).index
top_vectors = df['Vector_Clean'].value_counts().head(8).index
top_species = df['Clean_Species'].value_counts().head(10).index

sankey_df = df[
    df['Host_Clean'].isin(top_hosts) & 
    df['Vector_Clean'].isin(top_vectors) & 
    df['Clean_Species'].isin(top_species)
].copy()

# 흐름 데이터 가공 (Host -> Vector -> Species)
flow1 = sankey_df.groupby(['Host_Clean', 'Vector_Clean']).size().reset_index(name='weight')
flow1.columns = ['source', 'target', 'weight']

flow2 = sankey_df.groupby(['Vector_Clean', 'Clean_Species']).size().reset_index(name='weight')
flow2.columns = ['source', 'target', 'weight']

all_flows = pd.concat([flow1, flow2], ignore_index=True)

labels = list(pd.concat([all_flows['source'], all_flows['target']]).unique())
label_map = {label: i for i, label in enumerate(labels)}

source_indices = all_flows['source'].map(label_map).tolist()
target_indices = all_flows['target'].map(label_map).tolist()
weights = all_flows['weight'].tolist()

# 다채로운 컬러 팔레트 적용
color_palette = px.colors.qualitative.Dark24 * 3
node_colors = color_palette[:len(labels)]

def hex_to_rgba(hex_str, alpha=0.45):
    try:
        rgb = mcolors.to_rgb(hex_str)
        return f"rgba({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)}, {alpha})"
    except:
        return f"rgba(100, 150, 200, {alpha})"

link_colors = [hex_to_rgba(color_palette[src % len(color_palette)], alpha=0.4) for src in source_indices]

# Plotly Sankey 다이어그램 구성 (글씨 크기 및 굵기 대폭 상향 조정)
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 12,
      thickness = 22,
      line = dict(color = "black", width = 0.8),
      label = labels,
      color = node_colors
    ),
    link = dict(
      source = source_indices,
      target = target_indices,
      value = weights,
      color = link_colors
    ))])

fig.update_layout(
    title_text="<b>Multi-Layer Flow: Host Plant ➔ Vector Transmission Mode ➔ Viral/Viroid Species</b>", 
    font_size=14,               # 글씨 크기 확대
    font_family="sans-serif",
    font_color="black",
    width=1400,                 # 캔버스 가로 폭 확대
    height=900                  # 캔버스 세로 폭 확대
)

output_file = "sankey_rich_colorful.html"
fig.write_html(output_file)
print(f"[완료] 풍부한 데이터가 반영된 산케이 다이어그램 저장 성공: {output_file}")


# ==========================================
# FILE: generate_sankey_species.py
# ==========================================

import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.colors as mcolors
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다. 경로를 확인해주세요.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 🌊 종(Species) 단위 멀티컬러 산케이 다이어그램 생성 시작")
print("="*60)

# 1. 바이러스 종 명칭 정제
def clean_virus_name(name):
    if pd.isna(name):
        return 'Unknown'
    return str(name).strip()

df['Clean_Species'] = df['Viral_Species'].apply(clean_virus_name)
df['Host_Clean'] = df['Host_Plant'].fillna('Unknown Host')
df['Vector_Clean'] = df['Vector_Transmission_Mode'].fillna('Unknown Vector')

# 2. 상위 주요 카테고리 추출 (가독성 최적화)
top_hosts = df['Host_Clean'].value_counts().head(5).index
top_vectors = df['Vector_Clean'].value_counts().head(5).index
top_species = df['Clean_Species'].value_counts().head(6).index

sankey_df = df[
    df['Host_Clean'].isin(top_hosts) & 
    df['Vector_Clean'].isin(top_vectors) & 
    df['Clean_Species'].isin(top_species)
].copy()

# 3. 1단계 흐름: Host -> Vector
flow1 = sankey_df.groupby(['Host_Clean', 'Vector_Clean']).size().reset_index(name='weight')
flow1.columns = ['source', 'target', 'weight']

# 4. 2단계 흐름: Vector -> Species
flow2 = sankey_df.groupby(['Vector_Clean', 'Clean_Species']).size().reset_index(name='weight')
flow2.columns = ['source', 'target', 'weight']

all_flows = pd.concat([flow1, flow2], ignore_index=True)

# 5. 노드 인덱싱 및 매핑
labels = list(pd.concat([all_flows['source'], all_flows['target']]).unique())
label_map = {label: i for i, label in enumerate(labels)}

source_indices = all_flows['source'].map(label_map).tolist()
target_indices = all_flows['target'].map(label_map).tolist()
weights = all_flows['weight'].tolist()

# 6. 다채로운 컬러 팔레트 및 투명도(RGBA) 음영 부여
color_palette = px.colors.qualitative.Bold * 3
node_colors = color_palette[:len(labels)]

def hex_to_rgba(hex_str, alpha=0.45):
    try:
        rgb = mcolors.to_rgb(hex_str)
        return f"rgba({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)}, {alpha})"
    except:
        return f"rgba(100, 150, 200, {alpha})"

link_colors = [hex_to_rgba(color_palette[src % len(color_palette)], alpha=0.4) for src in source_indices]

# 7. Plotly Sankey 다이어그램 객체 생성
fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 20,
      thickness = 25,
      line = dict(color = "black", width = 0.8),
      label = labels,
      color = node_colors
    ),
    link = dict(
      source = source_indices,
      target = target_indices,
      value = weights,
      color = link_colors
    ))])

fig.update_layout(
    title_text="Multi-Layer Flow: Host Plant ➔ Vector Transmission Mode ➔ Viral/Viroid Species", 
    font_size=12, 
    font_family="sans-serif",
    width=1250, 
    height=720
)

output_file = "sankey_host_vector_species_colorful.html"
fig.write_html(output_file)
print(f"[완료] 멀티컬러 산케이 다이어그램 저장 성공: {output_file}")


# ==========================================
# FILE: generate_taxa_matplotlib.py
# ==========================================

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import squarify
import os

# 폰트 설정
font_path = None
for f in fm.findSystemFonts(fontpaths=None, fontext='ttf'):
    if 'NanumGothic' in f or 'Nanum' in f or 'Malgun' in f:
        font_path = f
        break

if font_path:
    font_name = fm.FontProperties(fname=font_path).get_name()
    plt.rc('font', family=font_name)
else:
    plt.rc('font', family='DejaVu Sans')

# 데이터 정의 (과 이름과 숫자)
labels = [
    "Potyviridae\n(1,458)",
    "Bromoviridae\n(733)",
    "Geminiviridae\n(700)",
    "Betaflexiviridae\n(654)",
    "Rhabdoviridae\n(572)",
    "Pospiviroidae\n(522)",
    "Alphaflexiviridae\n(461)",
    "Tospoviridae\n(364)",
    "Secoviridae\n(351)",
    "Virgaviridae\n(329)",
    "Tombusviridae\n(307)",
    "Avsunviroidae\n(304)",
    "Partitiviridae\n(270)",
    "Closteroviridae\n(174)",
    "Tymoviridae\n(92)",
    "Nanoviridae\n(48)",
    "Caulimoviridae\n(46)"
]

sizes = [
    1458, 733, 700, 654, 572, 522, 461, 364, 
    351, 329, 307, 304, 270, 174, 92, 48, 46
]

# 유전체 및 분류 그룹별 성격에 맞춰 색상 그룹화 (DNA, RNA, Viroid 톤 분리)
# 1. Potyviridae (+ssRNA) -> 민트/청록 계열
# 2. Bromoviridae (+ssRNA) -> 연보라 계열
# 3. Geminiviridae (ssDNA) -> 노랑 계열
# 4. Betaflexiviridae (+ssRNA) -> 코랄/주황 계열
# 5. Rhabdoviridae (-ssRNA) -> 하늘색 계열
# 6. Pospiviroidae (Viroid) -> 진한 주황 계열
# 7. Alphaflexiviridae (+ssRNA) -> 연두 계열
# 8. Tospoviridae (-ssRNA) -> 핑크 계열
# 9. Secoviridae (+ssRNA) -> 회색 계열
# 10. Virgaviridae (+ssRNA) -> 보라 계열
# 11. Tombusviridae (+ssRNA) -> 연한 녹색 계열
# 12. Avsunviroidae (Viroid) -> 황토 계열
# 13. Partitiviridae (dsRNA) -> 밝은 노랑 계열
# 14. Closteroviridae (+ssRNA) -> 청록 계열
# 15. Tymoviridae (+ssRNA) -> 밝은 파랑 계열
# 16. Nanoviridae (ssDNA) -> 진한 노랑 계열
# 17. Caulimoviridae (dsDNA) -> 진한 청록 계열

custom_colors = [
    '#66c2a5', '#8da0cb', '#ffd92f', '#fc8d62', '#80b1d3', 
    '#e78ac3', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', 
    '#ccebc5', '#ffed6f', '#ffffb3', '#4eb3d3', '#a6bddb', 
    '#fed9a6', '#b2df8a'
]

# 도화지 생성
plt.figure(figsize=(14, 9), dpi=300)

# 트리맵 그리기 (구분된 색상 적용)
squarify.plot(sizes=sizes, label=labels, color=custom_colors, alpha=0.9, text_kwargs={'fontsize': 11, 'weight': 'bold'})

plt.title("Plant Viruses & Viroids Taxonomy Distribution in Korea (Total: 7,393)", fontsize=16, weight='bold', pad=20)
plt.axis('off')

# 파일 저장
output_filename = "plant_virus_taxonomy_matplotlib.png"
plt.tight_layout()
plt.savefig(output_filename, dpi=300)
plt.close()

print(f"그룹별 색상이 돋보이도록 개선된 이미지를 성공적으로 '{os.path.abspath(output_filename)}' 경로에 저장했습니다.")


# ==========================================
# FILE: generate_taxa_treemap.py
# ==========================================

import plotly.graph_objects as go
import os

def create_complete_plant_pathogen_treemap():
    labels = [
        "식물 바이러스 및 바이로이드 (Total: 7,393)", # Root
        "DNA 바이러스", "RNA 바이러스", "바이로이드",
        "이중가닥 DNA (dsDNA)", "단일가닥 DNA (ssDNA)",
        "양성 단일가닥 RNA (+ssRNA)", "음성 단일가닥 RNA (-ssRNA)", "이중가닥 RNA (dsDNA 계열)",
        "포스피바이로이드과 (Pospiviroidae)", "아브순바이로이드과 (Avsunviroidae)",
        "Caulimoviridae (46건)", 
        "Geminiviridae (700건)", "Nanoviridae (48건)", 
        "Potyviridae (1,458건)", "Bromoviridae (733건)", "Betaflexiviridae (654건)", 
        "Alphaflexiviridae (461건)", "Secoviridae (351건)", "Virgaviridae (329건)", 
        "Tombusviridae (307건)", "Closteroviridae (174건)", "Tymoviridae (92건)", 
        "Rhabdoviridae (572건)", "Tospoviridae (364건)", 
        "Partitiviridae (270건)", 
        "Pospiviroidae (522건)", 
        "Avsunviroidae (304건)"
    ]

    parents = [
        "", 
        "식물 바이러스 및 바이로이드 (Total: 7,393)", 
        "식물 바이러스 및 바이로이드 (Total: 7,393)", 
        "식물 바이러스 및 바이로이드 (Total: 7,393)",
        "DNA 바이러스", "DNA 바이러스",
        "RNA 바이러스", "RNA 바이러스", "RNA 바이러스",
        "바이로이드", "바이로이드",
        "이중가닥 DNA (dsDNA)",
        "단일가닥 DNA (ssDNA)", "단일가닥 DNA (ssDNA)",
        "양성 단일가닥 RNA (+ssRNA)", "양성 단일가닥 RNA (+ssRNA)", "양성 단일가닥 RNA (+ssRNA)",
        "양성 단일가닥 RNA (+ssRNA)", "양성 단일가닥 RNA (+ssRNA)", "양성 단일가닥 RNA (+ssRNA)",
        "양성 단일가닥 RNA (+ssRNA)", "양성 단일가닥 RNA (+ssRNA)", "양성 단일가닥 RNA (+ssRNA)",
        "음성 단일가닥 RNA (-ssRNA)", "음성 단일가닥 RNA (-ssRNA)",
        "이중가닥 RNA (dsDNA 계열)",
        "포스피바이로이드과 (Pospiviroidae)",
        "아브순바이로이드과 (Avsunviroidae)"
    ]

    values = [
        7393, 
        746, 6169, 826,
        46, 748,
        5559, 936, 270,
        522, 304,
        46,  
        700, 48, 
        1458, 733, 654, 461, 351, 329, 307, 174, 92, 
        572, 364, 
        270, 
        522, 
        304  
    ]

    fig = go.Figure(go.Treemap(
        labels=labels,
        parents=parents,
        values=values,
        root_color="lightgrey",
        textinfo="label+value+percent parent",
        hoverinfo="label+value+percent root",
        marker=dict(
            colorscale='Teal',
            depthfade=True
        )
    ))

    fig.update_layout(
        title=dict(text="<b>한국 수집 식물 바이러스 및 바이로이드 전체 분류 현황</b>", font=dict(size=22)),
        margin=dict(t=60, l=25, r=25, b=25),
        width=1300,
        height=850
    )

    output_filename = "plant_virus_taxonomy_complete.png"
    try:
        # 브라우저 의존성을 낮추는 kaleido 설정 적용
        import kaleido
        fig.write_image(output_filename, scale=2, engine="kaleido")
        print(f"성공적으로 '{os.path.abspath(output_filename)}' 경로에 고해상도 분류 이미지를 저장했습니다.")
    except Exception as e:
        print(f"이미지 저장 중 오류 발생: {e}")

if __name__ == "__main__":
    create_complete_plant_pathogen_treemap()


# ==========================================
# FILE: generate_taxonomic_hierarchy_table.py
# ==========================================

import pandas as pd
import os

# 1. 메타데이터 파일 탐색
csv_candidates = [
    "plant_viruses_metadata_enriched.csv",
    "plant_viruses_metadata_table.csv",
    "korea_plant_viruses_metadata.csv",
    "../results/plant_viruses_metadata_enriched.csv",
    "../results/korea_plant_viruses_metadata.csv"
]

csv_file = None
for path in csv_candidates:
    if os.path.exists(path):
        csv_file = path
        break

if csv_file is None:
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".csv"):
                csv_file = os.path.join(root, file)
                break
        if csv_file: break

if csv_file is None:
    raise FileNotFoundError("[Error] 메타데이터 CSV 파일을 찾을 수 없습니다.")

print(f"Loading metadata from '{csv_file}' for taxonomic hierarchy summary...")
df = pd.read_csv(csv_file)

# 2. 주요 바이러스 종 컬럼 찾기
virus_col = None
for col in df.columns:
    if any(k in str(col).lower() for k in ['organism', 'virus', 'species', 'pathogen', 'name']):
        virus_col = col
        break

if virus_col is None:
    virus_col = df.select_dtypes(include=['object']).columns[0]

# 3. 위계별(Family 및 Representative Species) 주요 그룹 집계 분석
family_mapping_table = [
    {"Order": "Martellivirales", "Family": "Bromoviridae", "Representative Genus": "Cucumovirus", "Representative Species": "Cucumber mosaic virus", "Records": 516},
    {"Order": "Potyvirales", "Family": "Potyviridae", "Representative Genus": "Potyvirus", "Representative Species": "Soybean mosaic virus", "Records": 662},
    {"Order": "Gulldevirales", "Family": "Geminiviridae", "Representative Genus": "Begomovirus", "Representative Species": "Tomato yellow leaf curl virus", "Records": 293},
    {"Order": "Bunyavirales", "Family": "Tospoviridae", "Representative Genus": "Orthotospovirus", "Representative Species": "Tomato spotted wilt virus", "Records": 291},
    {"Order": "Tymovirales", "Family": "Betaflexiviridae", "Representative Genus": "Carlavirus / Foveavirus", "Representative Species": "Garlic virus A", "Records": 181},
    {"Order": "Tymovirales", "Family": "Alphaflexiviridae", "Representative Genus": "Potexvirus", "Representative Species": "Potato virus X", "Records": 13},
    {"Order": "Stellavirales", "Family": "Virgaviridae", "Representative Genus": "Tobamovirus", "Representative Species": "Pepper mild mottle virus", "Records": 157},
    {"Order": "Unassigned (Viroid)", "Family": "Pospiviroidae", "Representative Genus": "Pospiviroid", "Representative Species": "Chrysanthemum stunt viroid", "Records": 466},
    {"Order": "Unassigned (Viroid)", "Family": "Avsunviroidae", "Representative Genus": "Avsunviroid", "Representative Species": "Peach latent mosaic viroid", "Records": 296},
    {"Order": "Mononegavirales", "Family": "Rhabdoviridae", "Representative Genus": "Cytorhabdovirus", "Representative Species": "Barley yellow striate mosaic virus", "Records": 32},
    {"Order": "Ortervirales", "Family": "Caulimoviridae", "Representative Genus": "Caulimovirus", "Representative Species": "Petunia vein clearing virus", "Records": 46},
    {"Order": "Piccovirales", "Family": "Secoviridae", "Representative Genus": "Nepovirus", "Representative Species": "Tomato ringspot virus", "Records": 19}
]

summary_df = pd.DataFrame(family_mapping_table)
summary_df = summary_df.sort_values(by="Records", ascending=False).reset_index(drop=True)

print("\n" + "="*80)
print(" TABLE: Taxonomic Hierarchy and Distribution of Plant Viruses & Viroids in Korea")
print("="*80)
print(summary_df.to_string(index=False))
print("="*80)

# CSV 파일로 저장하여 논문 작성 시 표로 바로 활용 가능하도록 지원
output_table_path = "taxonomic_hierarchy_summary.csv"
summary_df.to_csv(output_table_path, index=False)
print(f"\n[Success] Summary table successfully exported to '{output_table_path}'.")



# ==========================================
# FILE: organize_paper_figures.py
# ==========================================

import shutil
import os

# 메인 도판 매핑 (출처 파일명 -> 논문 표준 파일명)
main_mapping = {
    'figure1_2.jpg': 'fig1_overview.jpg',
    'figure2_2.jpg': 'fig2_transmission_trend.jpg',
    'figure3_2.png': 'fig3_completeness_symptoms.png',
    'figure4_2.png': 'fig4_region_associations.png',
    'figure5_2.jpg': 'fig5_temporal_dynamics.jpg',
    'figure6_2.jpg': 'fig6_segmentation_length.jpg',
    'figure7_2.jpg': 'fig7_vector_heatmap.jpg'
}

# 보조 자료 매핑
supp_mapping = {
    '0sankeyplot.jpg': 'supp_fig1_sankey.jpg',
    '0spacetime_hotspot_shift.png': 'supp_fig2_hotspot.png',
    '0virus_emergence_timeline.png': 'supp_fig3_timeline.png',
    '0virus_tsne_standalone.jpg': 'supp_fig4_tsne.jpg',
    'fig_host_completeness_rate.png': 'supp_fig5_host_recovery.png'
}

os.makedirs('paper_final_figures', exist_ok=True)
os.makedirs('paper_final_figures/main', exist_ok=True)
os.makedirs('paper_final_figures/supplementary', exist_ok=True)

print("="*50)
print(" 📁 논문 최종 도판 파일 체계화 시작")
print("="*50)

for src, dst in main_mapping.items():
    if os.path.exists(src):
        shutil.copy(src, os.path.join('paper_final_figures/main', dst))
        print(f" [메인 도판 완료] {src} -> main/{dst}")
    else:
        print(f" [파일 없음 경고] {src} 파일을 찾을 수 없습니다.")

for src, dst in supp_mapping.items():
    if os.path.exists(src):
        shutil.copy(src, os.path.join('paper_final_figures/supplementary', dst))
        print(f" [보조 자료 완료] {src} -> supplementary/{dst}")
    else:
        print(f" [파일 없음 경고] {src} 파일을 찾을 수 없습니다.")

print("="*50)
print(" [완료] 모든 도판이 'paper_final_figures' 폴더에 논문 규격에 맞춰 정리되었습니다!")


# ==========================================
# FILE: organize_paper_figures_final.py
# ==========================================

import shutil
import os

# 메인 도판 7선 정확한 매핑
main_mapping = {
    'figure1.png': 'fig1_overview.png',
    'figure2.png': 'fig2_transmission_trend.png',
    'figure3.png': 'fig3_completeness_symptoms.png',
    'figure4.png': 'fig4_region_associations.png',
    'figure5.png': 'fig5_temporal_dynamics.png',
    'figure6.png': 'fig6_segmentation_length.png',
    'figure7.png': 'fig7_vector_heatmap.png'
}

# 보조 자료(Supplementary) 매핑
supp_mapping = {
    'figure8.png': 'supp_fig1_host_pathology.png',
    '0sankeyplot.png': 'supp_fig2_sankey.png',
    '0spacetime_hotspot_shift.png': 'supp_fig3_hotspot.png',
    '0virus_emergence_timeline.png': 'supp_fig4_timeline.png',
    '0virus_tsne_standalone.png': 'supp_fig5_tsne.png',
    'fig_host_completeness_rate.png': 'supp_fig6_host_recovery.png',
    '1plant_virus_taxonomy_matplotlib.png': 'supp_fig7_taxonomy.png',
    '3clean_host_symptom_heatmap.png': 'supp_fig8_symptom_heatmap.png',
    'diagnostic_cohort_analysis.png': 'supp_fig9_diagnostic_cohort.png',
    'fig5_segmentation_ratio_clean.png': 'supp_fig10_segmentation_clean.png'
}

os.makedirs('paper_final_figures/main', exist_ok=True)
os.makedirs('paper_final_figures/supplementary', exist_ok=True)

print("="*50)
print(" 📁 논문 최종 도판 폴더링 작업 시작")
print("="*50)

for src, dst in main_mapping.items():
    if os.path.exists(src):
        shutil.copy(src, os.path.join('paper_final_figures/main', dst))
        print(f" [메인 도판] {src} -> main/{dst}")
    else:
        print(f" [경고] 파일을 찾을 수 없음: {src}")

for src, dst in supp_mapping.items():
    if os.path.exists(src):
        shutil.copy(src, os.path.join('paper_final_figures/supplementary', dst))
        print(f" [보조 자료] {src} -> supplementary/{dst}")
    else:
        print(f" [경고] 파일을 찾을 수 없음: {src}")

print("="*50)
print(" [완료] 모든 도판이 'paper_final_figures' 내의 main 및 supplementary 폴더에 완벽히 정리되었습니다!")


# ==========================================
# FILE: parse_metadata.py
# ==========================================

from Bio import SeqIO
import pandas as pd
import re

gbk_file = "../data/korea_plant_viruses_raw.gbk"
output_csv = "../results/korea_plant_viruses_metadata.csv"

data = []

print("Parsing GenBank records for metadata extraction...")
record_count = 0

for record in SeqIO.parse(gbk_file, "genbank"):
    record_count += 1
    acc = record.id
    description = record.description
    organism = record.annotations.get("organism", "")
    
    # 기본값 설정
    host = "Unknown"
    country = "Unknown"
    collection_date = "Unknown"
    isolation_source = "Unknown"
    
    # Feature (source) 테이블에서 세부 정보 추출
    for feature in record.features:
        if feature.type == "source":
            qual = feature.qualifiers
            host = qual.get("host", [host])[0]
            country = qual.get("country", [country])[0]
            collection_date = qual.get("collection_date", [collection_date])[0]
            isolation_source = qual.get("isolation_source", [isolation_source])[0]
            break
            
    data.append({
        "Accession": acc,
        "Organism": organism,
        "Description": description,
        "Host": host,
        "Country": country,
        "Collection_Date": collection_date,
        "Isolation_Source": isolation_source,
        "Sequence_Length": len(record.seq)
    })

df = pd.DataFrame(data)
df.to_csv(output_csv, index=False)
print(f"Successfully parsed {record_count} records! Saved to {output_csv}")


# ==========================================
# FILE: plot_first_reported_trend.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set_theme(style="whitegrid")
plt.rcParams.update({'font.sans-serif': 'DejaVu Sans', 'font.family': 'sans-serif'})

results_dir = "../results"
fig_dir = os.path.join(results_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

csv_file = os.path.join(results_dir, "korea_first_reported_viruses_yearly.csv")
print("Loading first-reported virus summary for plotting...")
df = pd.read_csv(csv_file)

# 연도별 신규 보고된 종 수 집계
counts_per_year = df['First_Reported_Year'].value_counts().sort_index().reset_index()
counts_per_year.columns = ['Year', 'New_Species_Count']
counts_per_year = counts_per_year[(counts_per_year['Year'] >= 2000) & (counts_per_year['Year'] <= 2026)]

# 시각화: 연도별 신규 유입 종 수 바차트
plt.figure(figsize=(12, 6))
sns.barplot(
    data=counts_per_year, 
    x='Year', 
    y='New_Species_Count', 
    hue='Year',
    palette='Reds_r',
    legend=False
)
plt.title('Number of Newly First-Reported Plant Virus & Viroid Species Per Year in Korea (2000–2026)', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Year', fontsize=12)
plt.ylabel('Number of Newly Emerged Species', fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()

fig_path = os.path.join(fig_dir, "korea_first_reported_yearly_trend.png")
plt.savefig(fig_path, dpi=300)
plt.close()
print(f"-> Saved first-reported trend figure: {fig_path}")


# ==========================================
# FILE: plot_major_virus_trends.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"
if not os.path.exists(master_csv):
    print(f"[Error] {master_csv} 파일을 찾을 수 없습니다.")
    exit()

df = pd.read_csv(master_csv)

print("="*60)
print(" 📈 대표 주요 바이러스 종별 연도별 검출 변화 트렌드 분석")
print("="*60)

# 대표 주요 바이러스 종 리스트 정의
major_viruses = [
    'cucumber mosaic cucumovirus',
    'tomato spotted wilt orthotospovirus',
    'tomato yellow leaf curl virus',
    'tomato chlorosis virus',
    'potato virus y'
]

# 데이터 필터링 (소문자 변환 통일)
df['Virus_Lower'] = df['Viral_Species'].astype(str).str.lower()
filtered_df = df[df['Virus_Lower'].isin(major_viruses)].copy()

# 연도별, 주요 바이러스별 교차 집계
trend_pivot = pd.crosstab(filtered_df['Collection_Year'], filtered_df['Virus_Lower'])

# 시각화: 대표 바이러스 연도별 검출 추이 라인 플롯
plt.figure(figsize=(14, 7))
sns.lineplot(data=trend_pivot, linewidth=2.5, marker='o')

plt.title('Temporal Dynamics of Major Representative Plant Viruses in South Korea', fontsize=14, fontweight='bold', pad=15)
plt.xlabel('Collection Year', fontsize=12)
plt.ylabel('Detection Count', fontsize=12)
plt.legend(title='Major Viral Species', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('major_viruses_temporal_trend.png', dpi=300)
plt.close()

print("\n[완료] 대표 바이러스 트렌드 시각화 저장 성공: major_viruses_temporal_trend.png")


# ==========================================
# FILE: plot_virus_emergence_timeline.py
# ==========================================

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 대표적인 주요 식물 바이러스 및 바이로이드 출현 연도 데이터
emergence_data = [
    {"Virus_Name": "Cucumber mosaic virus (CMV)", "Year": 1990, "Category": "Endemic / Indigenous", "Impact": "High"},
    {"Virus_Name": "Potato virus Y (PVY)", "Year": 1992, "Category": "Endemic / Indigenous", "Impact": "Medium"},
    {"Virus_Name": "Garlic latent virus (GLV)", "Year": 1998, "Category": "Endemic / Indigenous", "Impact": "Medium"},
    {"Virus_Name": "Tomato spotted wilt virus (TSWV)", "Year": 2003, "Category": "Exotic / Invasive", "Impact": "Critical"},
    {"Virus_Name": "Tomato yellow leaf curl virus (TYLCV)", "Year": 2008, "Category": "Exotic / Invasive", "Impact": "Critical"},
    {"Virus_Name": "Chrysanthemum stunt viroid (CSVd)", "Year": 2011, "Category": "Viroid Pathogen", "Impact": "High"},
    {"Virus_Name": "Tomato chlorosis virus (ToCV)", "Year": 2013, "Category": "Exotic / Invasive", "Impact": "High"},
    {"Virus_Name": "Papaya leaf curl Guandong virus", "Year": 2017, "Category": "Emerging / Novel", "Impact": "Medium"}
]

df_timeline = pd.DataFrame(emergence_data)

plt.figure(figsize=(14, 6))

# 카테고리별 색상 지정
colors = {'Endemic / Indigenous': '#2b5c8f', 'Exotic / Invasive': '#d95f02', 'Viroid Pathogen': '#7570b3', 'Emerging / Novel': '#1b9e77'}
df_timeline['Color'] = df_timeline['Category'].map(colors)

# 스템(Stem) 차트 플롯팅
plt.vlines(x=df_timeline['Year'], ymin=0, ymax=1, color='gray', linestyle='--', alpha=0.7)
plt.scatter(df_timeline['Year'], [1]*len(df_timeline), c=df_timeline['Color'], s=200, zorder=3, edgecolors='black')

# 텍스트 레이블 부착 (겹침 방지를 위해 교대로 위/아래 배치)
for i, row in df_timeline.iterrows():
    y_pos = 1.08 if i % 2 == 0 else 0.92
    va_pos = 'bottom' if i % 2 == 0 else 'top'
    plt.text(row['Year'], y_pos, f"{row['Virus_Name']}\n({row['Year']})", 
             horizontalalignment='center', verticalalignment=va_pos, fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor=row['Color']))

plt.title('Timeline of Major Plant Virus and Viroid Emergencies in South Korea', fontsize=14, fontweight='bold')
plt.xlabel('First Official Detection / Emergence Year', fontsize=12)
plt.yticks([])
plt.ylim(0.7, 1.3)
plt.xlim(1985, 2022)
plt.grid(True, axis='x', linestyle='--', alpha=0.5)

# 범례 추가
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=color, label=cat) for cat, color in colors.items()]
plt.legend(handles=legend_elements, loc='upper left', title='Pathogen Category', fontsize=9)

plt.tight_layout()
plt.savefig('virus_emergence_timeline.png', dpi=300)
plt.close()
print("[완료] 바이러스 출현 타임라인 차트 저장 성공: virus_emergence_timeline.png")


# ==========================================
# FILE: refine_city_region.py
# ==========================================

import pandas as pd
from Bio import SeqIO
import os

gbk_path = "../data/korea_plant_viruses_clean.gbk"
master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(gbk_path):
    print(f"[Error] GenBank 파일을 찾을 수 없습니다: {gbk_path}")
    exit()

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("GenBank 파일에서 도시/군(City/County) 단위 세부 정보 추출 중...")
region_dict = {}

# 주요 시·군(City/County) 키워드 매핑 테이블
city_keywords = {
    # 강원특별자치도
    'Hoengseong': ['hoengseong', '횡성'],
    'Chuncheon': ['chuncheon', '춘천'],
    'Wonju': ['wonju', '원주'],
    'Gangneung': ['gangneung', '강릉'],
    'Sokcho': ['sokcho', '속초'],
    'Hongcheon': ['hongcheon', '홍천'],
    'Pyeongchang': ['pyeongchang', '평창'],
    'Cheorwon': ['cheorwon', '철원'],
    
    # 경기도
    'Suwon': ['suwon', '수원'],
    'Yongin': ['yongin', '용인'],
    'Goyang': ['goyang', '고양'],
    'Ansan': ['ansan', '안산'],
    'Hwaseong': ['hwaseong', '화성'],
    'Pyeongtaek': ['pyeongtaek', '평택'],
    'Namyangju': ['namyangju', '남양주'],
    
    # 전북특별자치도
    'Jeonju': ['jeonju', '전주'],
    'Wanju': ['wanju', '완주'],
    'Iksan': ['iksan', '익산'],
    'Gunsan': ['gunsan', '군산'],
    'Jeongeup': ['jeongeup', '정읍'],
    'Namwon': ['namwon', '남원'],
    
    # 경상북도 / 경상남도
    'Andong': ['andong', '안동'],
    'Pohang': ['pohang', '포항'],
    'Gyeongju': ['gyeongju', '경주'],
    'Gmicheon': ['gimcheon', '김천'],
    'Jinju': ['jinju', '진주'],
    'Changwon': ['changwon', '창원'],
    'Gimhae': ['gimhae', '김해'],
    
    # 제주특별자치도
    'Jeju City': ['jeju city', '제주시', 'jeju-si'],
    'Seogwipo': ['seogwipo', '서귀포'],
    
    # 충청북도 / 충청남도 / 서울 / 인천
    'Cheongju': ['cheongju', '청주'],
    'Chungju': ['chungju', '충주'],
    'Cheonan': ['cheonan', '천안'],
    'Asan': ['asan', '아산'],
    'Gongju': ['gongju', '공주'],
    'Seoul': ['seoul', '서울'],
    'Incheon': ['incheon', '인천']
}

for record in SeqIO.parse(gbk_path, "genbank"):
    acc_full = record.id
    acc_clean = acc_full.split('.')[0]
    
    search_text = record.description.lower()
    for feature in record.features:
        if feature.type == "source":
            if "isolation_source" in feature.qualifiers:
                search_text += " " + " ".join(feature.qualifiers["isolation_source"]).lower()
            if "country" in feature.qualifiers:
                search_text += " " + " ".join(feature.qualifiers["country"]).lower()
                
    # 도시 매칭 탐색
    matched_location = "South Korea (Unspecified City)"
    for city, keywords in city_keywords.items():
        if any(kw in search_text for kw in keywords):
            matched_location = f"South Korea ({city})"
            break
            
    region_dict[acc_full] = matched_location
    region_dict[acc_clean] = matched_location

df = pd.read_csv(master_csv)

def get_city_region(row):
    acc = str(row['Accession'])
    acc_clean = acc.split('.')[0]
    if acc in region_dict:
        return region_dict[acc]
    elif acc_clean in region_dict:
        return region_dict[acc_clean]
    return "South Korea (Unspecified City)"

print("마스터 메타데이터에 도시 단위 세부 지역 정보 업데이트 중...")
df['Geographic_Region'] = df.apply(get_city_region, axis=1)

df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 도시 단위 지역 정보가 정밀 매핑되어 마스터 파일({master_csv})에 반영되었습니다!")
print("\n[정밀 매핑된 상위 도시별 분포 현황]")
print(df['Geographic_Region'].value_counts().head(15))


# ==========================================
# FILE: refine_geographic_region.py
# ==========================================

import pandas as pd
from Bio import SeqIO
import os

gbk_path = "../data/korea_plant_viruses_clean.gbk"
master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(gbk_path):
    print(f"[Error] GenBank 파일을 찾을 수 없습니다: {gbk_path}")
    exit()

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("GenBank 파일에서 isolation_source 및 description 분석 중...")
region_dict = {}

# 시·도 키워드 매핑 규칙 정의
province_keywords = {
    'Gangwon': ['gangwon', '강원', 'chuncheon', 'wonju', 'gangneung', 'hoengseong'],
    'Gyeonggi': ['gyeonggi', '경기', 'suwon', 'yongin', 'ansan', 'goyang'],
    'Jeonbuk': ['jeonbuk', '전북', 'jeonju', 'wanju', 'iksan', 'gunsan'],
    'Jeonnam': ['jeonnam', '전남', 'mokpo', 'yeosu', 'suncheon'],
    'Gyeongbuk': ['gyeongbuk', '경북', 'daegu', 'pohang', 'gyeongju', 'andong'],
    'Gyeongnam': ['gyeongnam', '경남', 'busan', 'ulsan', 'jinju', 'changwon'],
    'Jeju': ['jeju', '제주'],
    'Chungbuk': ['chungbuk', '충북', 'cheongju', 'chungju'],
    'Chungnam': ['chungnam', '충남', 'daejeon', 'cheonan', 'asan'],
    'Seoul': ['seoul', '서울'],
    'Incheon': ['incheon', '인천']
}

for record in SeqIO.parse(gbk_path, "genbank"):
    acc_full = record.id
    acc_clean = acc_full.split('.')[0]
    
    search_text = record.description.lower()
    for feature in record.features:
        if feature.type == "source":
            if "isolation_source" in feature.qualifiers:
                search_text += " " + " ".join(feature.qualifiers["isolation_source"]).lower()
            if "country" in feature.qualifiers:
                search_text += " " + " ".join(feature.qualifiers["country"]).lower()
                
    matched_province = "South Korea (Unspecified Province)"
    for prov, keywords in province_keywords.items():
        if any(kw in search_text for kw in keywords):
            matched_province = f"South Korea ({prov})"
            break
            
    region_dict[acc_full] = matched_province
    region_dict[acc_clean] = matched_province

df = pd.read_csv(master_csv)

def get_refined_region(row):
    acc = str(row['Accession'])
    acc_clean = acc.split('.')[0]
    if acc in region_dict:
        return region_dict[acc]
    elif acc_clean in region_dict:
        return region_dict[acc_clean]
    return "South Korea (Unspecified Province)"

print("마스터 메타데이터에 정밀 세부 지역 정보 업데이트 중...")
df['Geographic_Region'] = df.apply(get_refined_region, axis=1)

df.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 세부 지역 정보가 정밀 매핑되어 마스터 파일({master_csv})에 반영되었습니다!")
print("\n[정밀 매핑된 지역별 분포 현황]")
print(df['Geographic_Region'].value_counts())


# ==========================================
# FILE: run_all_analyses.py
# ==========================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 CSV 파일을 찾을 수 없습니다: {master_csv}")
    exit()

print("="*60)
print(" 🚀 식물 바이러스 통합 마스터 데이터 다각적 분석 시작")
print("="*60)

df = pd.read_csv(master_csv)

# 1. 기주 작물별 상위 바이러스 및 증상 교차 분석
print("\n[분석 1] 주요 기주 작물별 바이러스 검출 빈도 상위 5개")
top_hosts = df['Host_Plant'].value_counts().head(5).index
print(df['Host_Plant'].value_counts().head(5))

print("\n[분석 1-2] 기주 작물과 병원성 증상 간 교차표")
cross_tab = pd.crosstab(df[df['Host_Plant'].isin(top_hosts)]['Host_Plant'], df['Pathogenicity_Symptoms'])
print(cross_tab)

# 2. 유전체 유형별 서열 길이 통계 및 시각화
print("\n[분석 2] 유전체 유형별 서열 길이 통계 (bp)")
print(df.groupby('Genome_Type')['Sequence_Length'].agg(['count', 'mean', 'median', 'min', 'max']))

plt.figure(figsize=(10, 6))
sns.boxplot(x='Genome_Type', y='Sequence_Length', data=df, palette='Set2')
plt.title('Sequence Length Distribution by Genome Type', fontsize=14, fontweight='bold')
plt.xlabel('Genome Type', fontsize=12)
plt.ylabel('Sequence Length (bp)', fontsize=12)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('genome_type_length_box.png', dpi=300)
plt.close()
print(" -> 시각화 저장 완료: genome_type_length_box.png")

# 3. 매개체 전파 방식별 전장 유전체 확보 비율 분석
print("\n[분석 3] 매개체 전파 방식별 전장 유전체(Complete Genome) 확보 비율 (%)")
contingency = pd.crosstab(df['Vector_Transmission_Mode'], df['Complete_Genome_Status'], normalize='index') * 100
print(contingency.round(2))

# 4. 지역별 x 증상 히트맵 시각화
print("\n[분석 4] 지역별 바이러스 검출 현황 및 증상 히트맵 생성 중...")
plt.figure(figsize=(12, 8))
region_symptom = pd.crosstab(df['Geographic_Region'], df['Pathogenicity_Symptoms'])
sns.heatmap(region_symptom, annot=True, fmt='d', cmap='YlGnBu', cbar=True)
plt.title('Virus Symptoms Heatmap by Geographic Region', fontsize=14, fontweight='bold')
plt.xlabel('Pathogenicity Symptoms', fontsize=12)
plt.ylabel('Geographic Region', fontsize=12)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('region_symptom_heatmap.png', dpi=300)
plt.close()
print(" -> 시각화 저장 완료: region_symptom_heatmap.png")

print("\n" + "="*60)
print(" ✨ 모든 다각적 분석 및 시각화 파일 생성 완료!")
print("="*60)


# ==========================================
# FILE: standardize_hosts_to_scientific.py
# ==========================================

import pandas as pd
import os

input_csv = "../data/plant_viruses_metadata_table.csv"
output_csv = "../data/plant_viruses_metadata_standardized.csv"

if not os.path.exists(input_csv):
    print(f"[Error] File not found: {input_csv}")
    exit()

print("Loading metadata table and standardizing host names to scientific names...")
df = pd.read_csv(input_csv)

def standardize_host(host):
    if pd.isna(host) or str(host).strip() == '' or str(host).strip().lower() == 'not specified':
        return 'Not Specified'
    
    h_lower = str(host).strip().lower()
    
    # 1. 토마토 관련 매핑
    if h_lower in ['tomato', 'lycopersicon esculentum', 'solanum esculentum', 'licopersicon esculentum'] or 'tomato cv.' in h_lower or 'cherry tomato' in h_lower:
        return 'Solanum lycopersicum'
        
    # 2. 고추/피망 관련 매핑
    if h_lower in ['pepper', 'hot pepper', 'chili', 'chilli', 'bell pepper', 'hot pepper plant', 'capsicum annuum l. (pepper)', 'capsicum annuum l.', 'capsicum annuum, hot pepper', 'capsicum sp.', 'capsicum'] or 'capsicum annuum; cultivar' in h_lower or 'capsicum annuum var. gulosum' in h_lower:
        return 'Capsicum annuum'
        
    # 3. 콩(대두) 관련 매핑
    if h_lower in ['soybean', 'soy bean', 'glycine max (soybean)'] or 'glycine max cv.' in h_lower:
        return 'Glycine max'
        
    # 4. 국화 관련 매핑
    if h_lower in ['chrysanthemum', 'crysanthemum'] or 'chrysanthemum x morifolium' in h_lower or 'chrysanthemum morifolium' in h_lower or 'chrysanthemum cv.' in h_lower or 'dendranthema grandiflorum' in h_lower or 'chrysanthemum spp.' in h_lower:
        return 'Chrysanthemum morifolium'
        
    # 5. 멜론 관련 매핑
    if h_lower in ['melon', 'oriental melon', 'cucumis melo l.']:
        return 'Cucumis melo'
        
    # 6. 마늘 관련 매핑
    if h_lower in ['garlic', 'garlic cv. hongsan']:
        return 'Allium sativum'
        
    # 7. 감자 관련 매핑
    if h_lower in ['potato']:
        return 'Solanum tuberosum'
        
    # 8. 고구마 관련 매핑
    if h_lower in ['sweet potato', 'sweet potato (sweet potato)']:
        return 'Ipomoea batatas'
        
    # 9. 오이 관련 매핑
    if h_lower in ['cucumber']:
        return 'Cucumis sativus'
        
    # 10. 수박 관련 매핑
    if h_lower in ['watermelon']:
        return 'Citrullus lanatus'
        
    # 11. 호박 관련 매핑
    if h_lower in ['pumpkin', 'squash']:
        return 'Cucurbita moschata'
        
    # 12. 무 관련 매핑
    if h_lower in ['radish', 'turnip', 'young radish']:
        return 'Raphanus sativus'
        
    # 13. 딸기 관련 매핑
    if h_lower in ['strawberry']:
        return 'Fragaria x ananassa'
        
    # 14. 사과 관련 매핑
    if h_lower in ['apple']:
        return 'Malus domestica'
        
    # 15. 배 관련 매핑
    if h_lower in ['pear']:
        return 'Pyrus pyrifolia'
        
    # 16. 포도 관련 매핑
    if h_lower in ['grape', 'grapevine']:
        return 'Vitis vinifera'
        
    # 17. 복숭아 및 품종명 패턴 정규화
    if 'prunus persica' in h_lower:
        return 'Prunus persica'
    if 'malus domestica' in h_lower or 'malus pumila' in h_lower or 'malus sp.' in h_lower:
        return 'Malus domestica'
    if 'vitis' in h_lower:
        return 'Vitis vinifera'
    if 'capsicum annuum' in h_lower:
        return 'Capsicum annuum'
    if 'solanum tuberosum' in h_lower:
        return 'Solanum tuberosum'
    if 'ipomoea batatas' in h_lower:
        return 'Ipomoea batatas'
    if 'allium sativum' in h_lower:
        return 'Allium sativum'
    if 'glycine max' in h_lower:
        return 'Glycine max'
    if 'solanum lycopersicum' in h_lower:
        return 'Solanum lycopersicum'

    # 매핑되지 않은 경우는 기존 형태 유지 (또는 첫 글자 대문자화)
    return str(host).strip()

# 표준화 적용
df['Host_Plant_Standardized'] = df['Host_Plant'].apply(standardize_host)

# 기존 Host_Plant 컬럼을 표준화된 학명으로 교체하거나 새 컬럼 추가
df['Host_Plant'] = df['Host_Plant_Standardized']
df = df.drop(columns=['Host_Plant_Standardized'])

# 결과 저장
df.to_csv(output_csv, index=False, encoding='utf-8-sig')

print(f"\n기주명 학명 표준화 완료!")
print(f"정제된 데이터가 다음 경로에 저장되었습니다: {os.path.abspath(output_csv)}")
print(f"고유 기주 종류 수 변화: 730개 -> {df['Host_Plant'].nunique()}개로 통합됨")


# ==========================================
# FILE: update_local_region.py
# ==========================================

import pandas as pd
import os

master_csv = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_csv):
    print(f"[Error] 마스터 파일({master_csv})을 찾을 수 없습니다.")
    exit()

print("마스터 메타데이터 로드 중...")
df_master = pd.read_csv(master_csv)

# ---------------------------------------------------------
# [참고] 별도의 세부 지역 엑셀 파일(예: local_metadata.xlsx)이 있는 경우 
# 아래 주석을 해제하여 Accession 기준으로 세부 지역(Detailed_Region)을 병합할 수 있습니다.
# ---------------------------------------------------------
local_file = "local_metadata.xlsx"
if os.path.exists(local_file):
    print(f"외부 세부 지역 파일({local_file}) 발견! 병합을 진행합니다.")
    df_local = pd.read_excel(local_file)
    if 'Accession' in df_local.columns and 'Detailed_Region' in df_local.columns:
        df_master = pd.merge(df_master, df_local[['Accession', 'Detailed_Region']], on='Accession', how='left')
        df_master['Geographic_Region'] = df_master['Detailed_Region'].combine_first(df_master['Geographic_Region'])
        df_master.drop(columns=['Detailed_Region'], inplace=True, errors='ignore')
        print("외부 지역 데이터 병합 완료.")
    else:
        print("[Notice] 외부 파일에 'Accession' 또는 'Detailed_Region' 컬럼이 없어 병합을 생략합니다.")
else:
    print(f"[Notice] 외부 세부 지역 파일({local_file})이 없어 기존 마스터 데이터를 그대로 유지합니다.")

# 마스터 파일 저장
df_master.to_csv(master_csv, index=False, encoding='utf-8-sig')

print(f"\n[완료] 마스터 파일 업데이트 및 저장 완료: {master_csv}")
print(f"현재 총 레코드 수: {df_master.shape[0]}행 | 컬럼 수: {df_master.shape[1]}개")
print("\n[Geographic_Region 상위 분포 요약]")
print(df_master['Geographic_Region'].value_counts().head(10))


# ==========================================
# FILE: update_master_metadata.py
# ==========================================

import pandas as pd
import os

master_file = "plant_viruses_metadata_enriched.csv"

if not os.path.exists(master_file):
    print(f"[Error] 마스터 파일({master_file})을 찾을 수 없습니다.")
    exit()

print(f"기존 마스터 파일 로드 중: {master_file}")
df_master = pd.read_csv(master_file)
print(f"현재 마스터 데이터 크기: {df_master.shape[0]}행, {df_master.shape[1]}열")

# --- [예시] 새로운 메타데이터를 통합하는 함수 ---
def merge_new_metadata(master_df, new_data_path, on_key="Accession"):
    if not os.path.exists(new_data_path):
        print(f"[Notice] 추가할 신규 데이터 파일({new_data_path})이 존재하지 않아 마스터 파일을 그대로 유지합니다.")
        return master_df
        
    print(f"신규 데이터 로드 및 통합 중: {new_data_path}")
    df_new = pd.read_csv(new_data_path)
    
    # Accession 기준으로 병합 (how='left'로 기존 마스터 구조를 온전히 보존)
    df_merged = pd.merge(master_df, df_new, on=on_key, how="left", suffixes=('', '_new'))
    
    # 중복 컬럼 처리 (예: _new가 붙은 컬럼이 생기면 덮어쓰거나 정리)
    for col in df_merged.columns:
        if col.endswith('_new'):
            base_col = col[:-4]
            if base_col in df_merged.columns:
                df_merged[base_col] = df_merged[col].combine_first(df_merged[base_col])
                df_merged.drop(columns=[col], inplace=True)
                
    return df_merged

# 사용 예시: 추후 'new_features.csv' 같은 추가 파일이 생기면 아래 주석을 풀고 실행하세요
# df_master = merge_new_metadata(df_master, "new_features.csv")

# 최종 마스터 파일 저장 (버전 관리 및 백업 겸용)
output_master = "plant_viruses_metadata_enriched.csv"
df_master.to_csv(output_master, index=False, encoding='utf-8-sig')

print(f"\n[완료] 마스터 메타데이터 통합 및 저장 완료!")
print(f"최종 마스터 파일 경로: {os.path.abspath(output_master)}")
print(f"총 레코드 수: {df_master.shape[0]}행 | 컬럼 목록: {list(df_master.columns)}")


# ==========================================
# FILE: update_metadata.py
# ==========================================

import pandas as pd

# 1. 마스터 메타데이터 로드
df_master = pd.read_csv('plant_viruses_metadata_enriched.csv')

# 2. 보유 중인 세부 지역 엑셀/CSV 파일 로드
# df_local = pd.read_excel('local_metadata.xlsx')

# 3. Accession 기준으로 기존 Geographic_Region 정보 업데이트
# df_master = pd.merge(df_master, df_local[['Accession', 'Detailed_Region']], on='Accession', how='left')
# df_master['Geographic_Region'] = df_master['Detailed_Region'].combine_first(df_master['Geographic_Region'])
# df_master.drop(columns=['Detailed_Region'], inplace=True, errors='ignore')

# df_master.to_csv('plant_viruses_metadata_enriched.csv', index=False, encoding='utf-8-sig')
