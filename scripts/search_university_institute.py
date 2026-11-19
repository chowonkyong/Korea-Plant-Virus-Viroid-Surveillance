import pandas as pd
import re
import os

# 1. GenBank 원본 파일 경로 설정
gbk_path = 'korea_plant_viruses_raw.gbk'

def search_institutions_in_gbk(path):
    if not os.path.exists(path):
        print(f"Error: '{path}' not found in current directory.")
        return
    
    print("=" * 70)
    print(f" [GenBank Institution & University Search] Scanning: {path}")
    print("=" * 70)
    
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        entries = content.split('LOCUS')
        
        universities = []
        institutes = []
        
        for entry in entries:
            if not entry.strip():
                continue
            
            # /institution 태그 검색
            inst_match = re.search(r'/institution="([^"]+)"', entry, re.IGNORECASE)
            if inst_match:
                val = inst_match.group(1)
                if 'univ' in val.lower():
                    universities.append(val)
                else:
                    institutes.append(val)
            
            # AUTHORS 또는 SOURCE 주석 내 대학/연구소 명칭 검색
            authors_match = re.findall(r'(?:university|institute|college|academy|research center|국립대학교|대학교|연구소)[^,.\n]*', entry, re.IGNORECASE)
            for m in authors_match:
                if 'univ' in m.lower():
                    universities.append(m.strip())
                else:
                    institutes.append(m.strip())
                    
    print(f"\n[Detected Universities / Higher Education]:")
    if universities:
        print(pd.Series(universities).value_counts().head(20))
    else:
        print("No explicit university tags found in raw gbk text.")
        
    print(f"\n[Detected Institutes / Research Centers]:")
    if institutes:
        print(pd.Series(institutes).value_counts().head(20))
    else:
        print("No explicit institute tags found in raw gbk text.")

# 실행
search_institutions_in_gbk(gbk_path)
