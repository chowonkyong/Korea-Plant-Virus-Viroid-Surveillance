import pandas as pd
import re
import os

input_file = 'korea_plant_viruses_final_clean.csv'
if not os.path.exists(input_file):
    print(f"Error: '{input_file}' not found.")
    exit(1)

df = pd.read_csv(input_file)

# 1. 지리적 위치 정보(Location) 표준화 함수
def clean_location(loc):
    if pd.isna(loc) or not str(loc).strip():
        return 'Unspecified'
    loc_str = str(loc).strip()
    if loc_str.lower() in ['korea', 'republic of korea']:
        return 'South Korea'
    if 'choong' in loc_str.lower():
        return 'South Korea: Choong-Nam'
    if 'suwon' in loc_str.lower():
        return 'South Korea: Suwon'
    if 'jincheon' in loc_str.lower():
        return 'South Korea: Jincheon'
    if 'andong' in loc_str.lower():
        return 'South Korea: Andong'
    if 'gwang' in loc_str.lower() and 'ju' in loc_str.lower():
        return 'South Korea: Gwangju'
    return loc_str

# 2. 연구기관 정보(Institution) 표준화 함수
def clean_institution(inst):
    if pd.isna(inst) or not str(inst).strip():
        return 'Unspecified'
    inst_str = str(inst).strip()
    
    if any(k in inst_str for k in ['Chonnam', 'Yongbong']):
        return 'Chonnam National University'
    elif any(k in inst_str for k in ['Seoul National', 'Nowon-Gu, Seoul 01797', 'Heukseok-ro', 'Pyeongchang']):
        return 'Seoul National University'
    elif any(k in inst_str for k in ['Kyungpook']):
        return 'Kyungpook National University'
    elif any(k in inst_str for k in ['Andong', 'Gyeongdong']):
        return 'Andong National University'
    elif any(k in inst_str for k in ['Horticultural and Herbal Science', 'NIHHS']):
        return 'NIHHS, RDA'
    elif any(k in inst_str for k in ['Agricultural Science', 'NAAS', 'Nongsaengmyeong']):
        return 'NAAS, RDA'
    elif any(k in inst_str for k in ['Chung-Ang']):
        return 'Chung-Ang University'
    elif any(k in inst_str for k in ['Crop Science']):
        return 'National Institute of Crop Science, RDA'
    elif any(k in inst_str for k in ['Jeju']):
        return 'Jeju National University'
    elif any(k in inst_str for k in ['Kangwon']):
        return 'Kangwon National University'
    elif any(k in inst_str for k in ['Chungbuk']):
        return 'Chungbuk National University'
    elif 'JOURNAL' in inst_str or 'Sequencing Technology' in inst_str:
        if 'University' in inst_str:
            match = re.search(r'([A-Za-z\s]+University)', inst_str)
            if match:
                return match.group(1).strip()
        return 'General GenBank Submission'
        
    return inst_str

# 컬럼 적용
if 'Refined_Location_Info' in df.columns:
    df['Location_Clean'] = df['Refined_Location_Info'].apply(clean_location)
if 'Refined_Institution_Source' in df.columns:
    df['Institution_Clean'] = df['Refined_Institution_Source'].apply(clean_institution)

# 최종 정제된 메타데이터 저장
df.to_csv(input_file, index=False, encoding='utf-8-sig')

print("=" * 70)
print(f"Successfully standardized location & institution in: {input_file}")
print("Top Locations:")
print(df['Location_Clean'].value_counts().head(5))
print("\nTop Institutions:")
print(df['Institution_Clean'].value_counts().head(5))
print("=" * 70)
