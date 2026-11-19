import pandas as pd
import re

# 1. 코어 메타데이터 로드
df = pd.read_csv('korea_plant_viruses_core_perfect.csv')

print("=" * 70)
print(f" [Advanced Institution Extraction] Total Records: {len(df)}")
print("=" * 70)

# 2. 대학 및 연구소 관련 키워드 및 주요 대학 영문 약어 정의
inst_patterns = [
    'university', 'college', 'institute', 'hospital', 'center', 'lab',
    'CNU', 'KNU', 'SNU', 'GNU', 'JNU', 'PNU', 'KRIBB', 'NIAS', 'RDA'
]

# 3. Description 필드에서 기관명 및 연구소 관련 키워드 추출 함수
def extract_detailed_institution(desc):
    if pd.isna(desc):
        return 'Unspecified'
    found = []
    desc_str = str(desc)
    for pat in inst_patterns:
        # 대소문자 구분 없이 단어 또는 약어 매칭
        if re.search(r'\b' + pat + r'\b', desc_str, re.IGNORECASE):
            found.append(pat.upper())
    return ', '.join(sorted(set(found))) if found else 'Unspecified'

df['Extracted_Institution_Tags'] = df['Description'].apply(extract_detailed_institution)

# 결과 요약 출력
print("\n[Extracted Institution / Acronym Distribution]:")
dist = df['Extracted_Institution_Tags'].value_counts()
print(dist.head(15))

# 결과 저장
output_csv = "korea_plant_viruses_institutions_detailed.csv"
df.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"\nSuccessfully extracted and saved detailed institution data to: {output_csv}")
