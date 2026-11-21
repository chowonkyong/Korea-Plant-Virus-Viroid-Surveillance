import os
import pandas as pd

def enrich_master_dataset():
    core_file = "korea_plant_viruses_core_perfect.csv"
    enriched_file = "plant_viruses_metadata_enriched.csv"
    output_file = "korea_plant_viruses_core_perfect.csv"
    
    if not os.path.exists(core_file) or not os.path.exists(enriched_file):
        print("오류: 필수 파일('korea_plant_viruses_core_perfect.csv' 또는 'plant_viruses_metadata_enriched.csv')을 찾을 수 없습니다.")
        return
        
    df_core = pd.read_csv(core_file)
    df_enriched = pd.read_csv(enriched_file)
    print(f"기존 마스터 레코드 수: {len(df_core):,}개")
    
    # 추가할 가치 있는 컬럼 정의
    cols_to_add = ['Accession', 'Geographic_Region', 'Pathogenicity_Symptoms', 'Genome_Segmentation']
    df_subset = df_enriched[[c for c in cols_to_add if c in df_enriched.columns]].copy()
    
    # 중복 Accession 방지
    df_subset = df_subset.drop_duplicates(subset=['Accession'])
    
    # Accession 기준으로 병합 (Left Join)
    df_merged = pd.merge(df_core, df_subset, on='Accession', how='left')
    
    # 최종 이상적인 컬럼 순서 재배열 (새로 추가된 컬럼 포함)
    desired_order = [
        'Accession',
        'Organism',
        'Standardized_Host',
        'Family',
        'Genus',
        'Order',
        'Taxonomy_Lineage',
        'Genome_Type',
        'Genome_Segmentation',
        'Vector_Mode',
        'Pathogenicity_Symptoms',
        'Geographic_Region',
        'Collection_Year',
        'Length',
        'Completeness_Flag',
        'Description'
    ]
    
    existing_columns = [col for col in desired_order if col in df_merged.columns]
    remaining_columns = [col for col in df_merged.columns if col not in existing_columns]
    final_column_order = existing_columns + remaining_columns
    
    df_merged = df_merged[final_column_order]
    
    # 결측치 방어 (혹시 모를 NaN은 빈 문자열 또는 Unknown으로 처리)
    df_merged['Geographic_Region'] = df_merged['Geographic_Region'].fillna('South Korea (Unspecified City)')
    df_merged['Pathogenicity_Symptoms'] = df_merged['Pathogenicity_Symptoms'].fillna('Unspecified / General Symptoms')
    df_merged['Genome_Segmentation'] = df_merged['Genome_Segmentation'].fillna('Unspecified')
    
    # 최종 마스터 파일 저장
    df_merged.to_csv(output_file, index=False, encoding="utf-8-sig")
    
    print("=" * 60)
    print(f" 데이터 보강 및 최종 마스터 갱신 완료: {output_file}")
    print("=" * 60)
    print(f" - 최종 총 컬럼 수: {len(df_merged.columns)}개")
    print(" - 추가된 주요 컬럼:")
    print("   1. Genome_Segmentation (유전체 분절 구조)")
    print("   2. Pathogenicity_Symptoms (주요 병징 및 증상)")
    print("   3. Geographic_Region (세부 채집 지역)")

if __name__ == "__main__":
    enrich_master_dataset()
