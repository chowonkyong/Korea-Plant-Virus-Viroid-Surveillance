import pandas as pd
import os

# 1. 파일 존재 여부 확인 및 로드
fully_merged_path = 'korea_plant_viruses_fully_merged_metadata.csv'
advanced_merged_path = 'korea_plant_viruses_advanced_merged.csv'
clean_inst_path = 'korea_plant_viruses_clean_institutions.csv'

if not os.path.exists(fully_merged_path):
    print(f"Error: '{fully_merged_path}' not found.")
    exit(1)

df_fully = pd.read_csv(fully_merged_path)
print(f"Loaded {fully_merged_path}: {len(df_fully)} rows")

# 2. 고급 메타데이터(도시/지역 정보) 병합용 데이터 로드
if os.path.exists(advanced_merged_path):
    df_adv = pd.read_csv(advanced_merged_path)
    if 'Base_Accession' in df_adv.columns and 'Location_Info' in df_adv.columns:
        df_fully['Base_Accession'] = df_fully['Accession'].astype(str).str.split('.').str[0]
        # Location_Info 업데이트
        loc_map = dict(zip(df_adv['Base_Accession'], df_adv['Location_Info']))
        df_fully['Refined_Location_Info'] = df_fully['Base_Accession'].map(loc_map).fillna(df_fully['Extracted_Country'])
        print("Successfully integrated refined location/city metadata.")

# 3. 정제된 연구기관 정보 병합용 데이터 로드
if os.path.exists(clean_inst_path):
    df_inst = pd.read_csv(clean_inst_path)
    if len(df_inst) == len(df_fully):
        df_fully['Refined_Institution_Source'] = df_inst['Full_Institution']
        print("Successfully integrated refined institution source metadata.")
    else:
        # Accession 기준으로 매핑 시도
        print("Row count mismatch with institution file; appending or mapping by index.")
        df_fully['Refined_Institution_Source'] = df_inst['Full_Institution'].values[:len(df_fully)]

# 4. 최종 업데이트된 메타데이터 저장
output_filename = 'korea_plant_viruses_fully_merged_metadata_updated.csv'
df_fully.to_csv(output_filename, index=False, encoding='utf-8-sig')

print("=" * 70)
print(f"Successfully created updated metadata file: {output_filename}")
print(f"Total Columns: {list(df_fully.columns)}")
print("=" * 70)
