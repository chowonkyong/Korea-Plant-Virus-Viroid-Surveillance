import pandas as pd
import re
import os

# 1. 코어 메타데이터 파일 로드
csv_path = 'korea_plant_viruses_core_perfect.csv'
if not os.path.exists(csv_path):
    for root, dirs, files in os.walk('.'):
        if 'korea_plant_viruses_core_perfect.csv' in files:
            csv_path = os.path.join(root, 'korea_plant_viruses_core_perfect.csv')
            break

df_core = pd.read_csv(csv_path)
print("=" * 70)
print(f" [Advanced Metadata Extraction] Core CSV Records: {len(df_core)}")
print("=" * 70)

# 2. GenBank 파일 전체 스캔 및 다양한 위치/출처 태그 동적 추출 함수
def parse_gbk_advanced(gbk_path):
    extracted_records = []
    
    if not os.path.exists(gbk_path):
        print(f"Error: '{gbk_path}' not found in current directory.")
        return None

    print(f"Scanning GenBank file with flexible tag matching: {gbk_path} ...")
    with open(gbk_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        entries = content.split('LOCUS')
        
        for entry in entries:
            if not entry.strip():
                continue
            
            # Accession 번호 추출
            acc_match = re.search(r'ACCESSION\s+([A-Za-z0-9_.]+)', entry, re.IGNORECASE)
            if not acc_match:
                acc_match = re.search(r'Accession\s+([A-Za-z0-9_.]+)', entry, re.IGNORECASE)
            
            full_acc = acc_match.group(1) if acc_match else 'Unknown'
            base_acc = full_acc.split('.')[0]
            
            # 다양한 지리적/출처 관련 GenBank 태그 탐색
            country_match = re.search(r'/country="([^"]+)"', entry, re.IGNORECASE)
            geo_match = re.search(r'/geo_loc_name="([^"]+)"', entry, re.IGNORECASE)
            isolation_match = re.search(r'/isolation_source="([^"]+)"', entry, re.IGNORECASE)
            institution_match = re.search(r'/institution="([^"]+)"', entry, re.IGNORECASE)
            lab_host_match = re.search(r'/lab_host="([^"]+)"', entry, re.IGNORECASE)
            note_match = re.search(r'/note="([^"]+)"', entry, re.IGNORECASE)
            
            # 매칭된 값 중 유효한 것 선택
            country = country_match.group(1) if country_match else (geo_match.group(1) if geo_match else '')
            isolation = isolation_match.group(1) if isolation_match else ''
            institution = institution_match.group(1) if institution_match else ''
            lab_host = lab_host_match.group(1) if lab_host_match else ''
            note = note_match.group(1) if note_match else ''
            
            combined_text = f"Country/Geo: {country} | Isolation: {isolation} | Inst: {institution} | Host: {lab_host} | Note: {note}"
            
            extracted_records.append({
                'Base_Accession': base_acc,
                'GenBank_Accession': full_acc,
                'Location_Info': country if country else isolation,
                'Institution_Info': institution,
                'Lab_Host': lab_host,
                'Note_Info': note,
                'Raw_Blob': combined_text
            })
            
    result_df = pd.DataFrame(extracted_records)
    print(f"Successfully processed {len(result_df)} records.")
    return result_df

# 3. 실행 및 병합
gbk_filename = 'korea_plant_viruses_raw.gbk'
gbk_raw_df = parse_gbk_advanced(gbk_filename)

if gbk_raw_df is not None and not gbk_raw_df.empty:
    df_core['Base_Accession'] = df_core['Accession'].astype(str).str.split('.').str[0]
    
    merged_df = pd.merge(df_core, gbk_raw_df[['Base_Accession', 'Location_Info', 'Institution_Info', 'Lab_Host', 'Note_Info', 'Raw_Blob']], on='Base_Accession', how='left')
    
    output_csv = "korea_plant_viruses_advanced_merged.csv"
    merged_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    
    print("=" * 70)
    print(f"Successfully saved advanced merged metadata to: {output_csv}")
    print("=" * 70)
    
    print("\n[Location Info Distribution (Non-empty)]:")
    loc_counts = merged_df['Location_Info'].replace('', pd.NA).dropna().value_counts()
    print(loc_counts.head(15) if not loc_counts.empty else "All location fields are currently unpopulated in raw gbk tags.")
else:
    print("Extraction failed.")
