import pandas as pd
import re

input_csv = 'korea_plant_viruses_clean_institutions.csv'
df = pd.read_csv(input_csv)

print(f"Loaded {len(df)} rows from {input_csv}")
print("Available columns:", df.columns.tolist())

def clean_and_map_institution(row):
    # Combine all text columns available in the row
    combined_text = " ".join([str(val) for val in row.values if pd.notna(val)])
    text_lower = combined_text.lower()
    
    # Mapping specific addresses or keywords to target Universities/Institutes
    if 'yongbong-ro' in text_lower or 'chonnam' in text_lower:
        return 'Chonnam National University'
    elif 'hwarangro' in text_lower or 'seoul women' in text_lower:
        return 'Seoul Women\'s University'
    elif 'gyeongdong-ro' in text_lower or 'andong' in text_lower:
        return 'Andong National University'
    elif 'daehak-ro' in text_lower or 'daehakro' in text_lower or 'kyungpook' in text_lower:
        if 'daegu' in text_lower or 'knu' in text_lower or 'buk-gu' in text_lower:
            return 'Kyungpook National University'
        return 'Kyungpook National University'
    elif 'heukseok-ro' in text_lower or 'chung-ang' in text_lower:
        return 'Chung-Ang University'
    elif 'suin-ro' in text_lower or 'nongsaengmyeong' in text_lower or 'national institute of agricultural sciences' in text_lower or 'nihas' in text_lower:
        return 'National Institute of Agricultural Sciences (NIAS)'
    elif 'horticultural' in text_lower or 'nihhs' in text_lower:
        return 'National Institute of Horticultural and Herbal Science (NIHHS)'
    elif 'cornell' in text_lower:
        return 'Cornell University'
    elif 'veterinary research institute' in text_lower or 'animal and plant quarantine agency' in text_lower:
        return 'Animal and Plant Quarantine Agency (APQA)'
    
    # Generic regex fallback for universities/institutes if address mapping didn't trigger
    match_univ = re.search(r'([A-Za-z\s\-]+(?:University|Univ\.?|College|Institute|Institution|Lab(?:oratory)?))', combined_text, re.IGNORECASE)
    if match_univ:
        val = match_univ.group(0).strip()
        if len(val) > 4 and not 'the' in val.lower()[:3]:
            return val
            
    match_kr = re.search(r'([가-힣\s\-]+(?:대학교|대학|연구소|연구원|진흥원|시험장|농업기술원|농업기술센터))', combined_text)
    if match_kr:
        return match_kr.group(0).strip()
        
    return ''

df['Standardized_Institution'] = df.apply(clean_and_map_institution, axis=1)

print("\n--- Top Standardized Institutions ---")
print(df['Standardized_Institution'].value_counts().head(20))

output_csv = 'korea_plant_viruses_institutions_detailed.csv'
df.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"\nSaved refined institution mapping to {output_csv}.")
