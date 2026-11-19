import pandas as pd
import os

input_file = 'korea_plant_viruses_fully_merged_metadata_updated.csv'
if os.path.exists(input_file):
    df = pd.read_csv(input_file)
    
    # 불필요한 중복 및 결측 컬럼 정의
    cols_to_drop = [
        'Extracted_Country', 
        'Extracted_Institution', 
        'Raw_Metadata_Blob', 
        'Base_Accession', 
        'Extracted_Lab_Host', 
        'Extracted_Note'
    ]
    
    # 존재하는 컬럼만 선택하여 드롭
    existing_cols_to_drop = [col for col in cols_to_drop if col in df.columns]
    df_cleaned = df.drop(columns=existing_cols_to_drop)
    
    # 최종 정제된 파일 저장
    output_file = 'korea_plant_viruses_final_clean.csv'
    df_cleaned.to_csv(output_file, index=False, encoding='utf-8-sig')
    
    print("=" * 70)
    print(f"Successfully cleaned metadata! Saved to: {output_file}")
    print(f"Total columns remaining: {len(df_cleaned.columns)}")
    print("Columns:", df_cleaned.columns.tolist())
    print("=" * 70)
else:
    print(f"Error: '{input_file}' not found.")
