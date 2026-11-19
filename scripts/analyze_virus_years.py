import pandas as pd

# 1. 실제 메타데이터 파일 로드
df = pd.read_csv('korea_plant_viruses_core_perfect.csv')

print("=" * 70)
print(f" [Metadata Year Analysis] 메타데이터 컬럼 확인 및 연도 관련 필드 조사")
print("=" * 70)

# 2. 컬럼명 전체 출력
print("\n[현재 메타데이터 컬럼 목록]:")
print(list(df.columns))

# 3. 연도(Year, Date, Collection_Date 등)와 관련된 컬럼이 있는지 탐색
year_columns = [col for col in df.columns if any(k in col.lower() for k in ['year', 'date', 'time', 'period'])]
print(f"\n[연도/날짜 관련 잠재적 컬럼]: {year_columns}")

if year_columns:
    for col in year_columns:
        print(f"\n- 컬럼 '{col}' 상위 샘플 값:")
        print(df[col].head(10))
else:
    print("\n[알림] 메타데이터 파일 내에 연도나 날짜를 직접 나타내는 명시적인 컬럼이 존재하지 않습니다.")
    print("샘플 ID, 파일명, 혹은 텍스트 필드에 연도 정보가 포함되어 있는지 확인합니다.")
    
    # 예시로 첫 번째 컬럼이나 텍스트 컬럼 스캔
    for col in df.columns:
        sample_vals = df[col].astype(str).str.contains('20|19', regex=True).sum()
        if sample_vals > 0:
            print(f"  > 컬럼 '{col}'에 연도 추정 패턴(19xx/20xx)이 {sample_vals}건 발견되었습니다.")

print("=" * 70)
