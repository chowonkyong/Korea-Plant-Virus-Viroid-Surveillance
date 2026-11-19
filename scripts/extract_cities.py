from Bio import SeqIO
import pandas as pd

gbk_file = 'korea_plant_viruses_raw.gbk'
records = []

print("Reading GenBank file and extracting geo_loc_name, organism, host, etc...")
for record in SeqIO.parse(gbk_file, 'genbank'):
    accession = record.id
    description = record.description
    
    geo_loc = ""
    organism = ""
    host = ""
    isolate = ""
    notes = []
    
    for feature in record.features:
        if feature.type == 'source':
            quals = feature.qualifiers
            geo_loc = quals.get('geo_loc_name', [''])[0]
            organism = quals.get('organism', [''])[0]
            host = quals.get('host', [''])[0]
            isolate = quals.get('isolate', [''])[0]
            notes.extend(quals.get('note', []))
            
    records.append({
        'Accession': accession,
        'Organism': organism,
        'Geo_Loc_Name': geo_loc,
        'Host': host,
        'Isolate': isolate,
        'Notes': " ; ".join(notes)
    })

df = pd.DataFrame(records)
print(f"\nTotal records parsed: {len(df)}")
print("\nTop 20 values in 'Geo_Loc_Name':")
print(df['Geo_Loc_Name'].value_counts().head(20))

output_csv = 'extracted_cities_from_gbk.csv'
df.to_csv(output_csv, index=False, encoding='utf-8-sig')
print(f"\nUpdated data saved to {output_csv}.")
