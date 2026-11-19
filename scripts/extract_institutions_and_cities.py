import pandas as pd
import re

input_csv = 'korea_plant_viruses_clean_institutions.csv'
df = pd.read_csv(input_csv)

print(f"Loaded {len(df)} rows from {input_csv}")

# Keyword dictionary mapping address markers or fragment terms to standardized institutions
def clean_and_map_institution(row):
    notes = str(row.get('Notes', ''))
    isolation = str(row.get('Isolation_Source', ''))
    
    # Check if there's already an extracted institution or text containing clues
    full_text = f"{notes} {isolation}"
    
    # Mapping specific addresses or keywords to target Universities/Institutes
    if 'yongbong-ro' in full_text.lower() or 'chonnam' in full_text.lower():
        return 'Chonnam National University'
    elif 'hwarangro' in full_text.lower() or 'seoul women' in full_text.lower():
        return 'Seoul Women\'s University'
    elif 'gyeongdong-ro' in full_text.lower() or 'andong' in full_text.lower():
        return 'Andong National University'
    elif 'daehak-ro' in full_text.lower() or 'daehakro' in full_text.lower() or 'kyungpook' in full_text.lower():
        if 'daegu' in full_text.lower() or 'knu' in full_text.lower():
            return 'Kyungpook National University'
    elif 'heukseok-ro' in full_text.lower() or 'chung-ang' in full_text.lower():
        return 'Chung-Ang University'
    elif 'suin-ro' in full_text.lower() or 'nongsaengmyeong' in full_text.lower() or 'national institute of agricultural sciences' in full_text.lower() or 'nihas' in full_text.lower() or 'national institute of horticultural' in full_text.lower():
        return 'National Institute of Agricultural Sciences (NIAS / NIHHS)'
    elif 'cornell' in full_text.lower():
        return 'Cornell University'
    elif 'veterinary research institute' in full_text.lower():
        return 'Animal and Plant Quarantine Agency (APQA / Veterinary Research Institute)'
    
    # Generic regex fallback for universities/institutes
    match_univ = re.search(r'([A-Za-z\s\-]+University|[A-Za-z\s\-]+Univ\.|[가-힣\s\-]+대학교)', full_text)
    if match_univ:
        return match_univ.group(0).strip()
        
    return ''

df['Standardized_Institution'] = df.apply(clean_and_map_institution, axis=1)

print("\n--- Top Standardized Institutions ---")
print(df['Standardized_Institution'].value_counts().head(15))

output_csv = 'korea_plant_viruses_institutions_detailed.csv'
df.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"\nSaved refined institution mapping to {output_csv}.")
