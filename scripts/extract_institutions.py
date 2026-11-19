import pandas as pd
import re

# 1. 코어 메타데이터 로드
df = pd.read_csv('korea_plant_viruses_core_perfect.csv')

print("=" * 70)
print(f" [Institution Extraction Script] Total Records: {len(df)}")
print("=" * 70)

# 2. 국내 주요 대학, 연구소, 농업기술원 관련 키워드 정의
institution_keywords = [
    'University', 'College', 'Institute', 'Academy', 'Hospital', 'Foundation',
    'National', 'Agricultural', 'Technology', 'RDA', 'Korean', 'Korea',
    '대학교', '연구소', '원예', '농업', '기술원', '진흥청', '국립'
]

# 3. Description 또는 관련 필드에서 기관명 후보 추출 함수
def extract_institution(desc):
    if pd.isna(desc):
        return 'Unspecified'
    found = [kw for kw in institution_keywords if kw.lower() in str(desc).lower()]
    return ', '.join(set(found)) if found else 'Unspecified'

df['Extracted_Institution_Keywords'] = df['Description'].apply(extract_institution)

# 결과 요약 출력
print("\n[Extracted Institution Keywords Distribution]:")
print(df['Extracted_Institution_Keywords'].value_counts().head(15))

# 결과 저장
output_csv = "korea_plant_viruses_institutions_extracted.csv"
df.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"\nSuccessfully extracted and saved institution metadata to: {output_csv}")
