import pandas as pd

# 1. 병합된 고급 메타데이터 파일 로드
df = pd.read_csv('korea_plant_viruses_advanced_merged.csv')

print("=" * 70)
print(f" [Geographic Metadata Analysis] Total Records: {len(df)}")
print("=" * 70)

# 2. Location_Info 정제 및 세부 지역 분리 함수
def clean_location(loc_str):
    if pd.isna(loc_str) or not str(loc_str).strip():
        return 'Unspecified'
    loc = str(loc_str).strip()
    if loc.lower() in ['korea', 'republic of korea']:
        return 'South Korea'
    return loc

df['Clean_Location'] = df['Location_Info'].apply(clean_location)

# 3. 전체 지역 분포 상위 30개 출력
print("\n[Top 30 Geographic Locations in Raw GenBank Metadata]:")
print(df['Clean_Location'].value_counts().head(30))

# 4. 'South Korea' 내부의 세부 도시/도(Province) 별 분포 집계
kr_subset = df[df['Clean_Location'].str.startswith('South Korea', na=False)]
print("\n[South Korea Internal Regional Breakdown]:")
print(kr_subset['Clean_Location'].value_counts().head(20))

# 세부 정리된 CSV 저장
output_cleaned_csv = "korea_plant_viruses_geography_cleaned.csv"
df.to_csv(output_cleaned_csv, index=False, encoding='utf-8-sig')
print(f"\nSuccessfully saved cleaned geography data to: {output_cleaned_csv}")
