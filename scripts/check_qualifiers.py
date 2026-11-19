from Bio import SeqIO

gbk_file = 'korea_plant_viruses_raw.gbk'

print("Inspecting qualifiers in the first 3 records...")
for i, record in enumerate(SeqIO.parse(gbk_file, 'genbank')):
    if i >= 3:
        break
    print(f"\n--- Record {i+1}: {record.id} ---")
    for feature in record.features:
        if feature.type == 'source':
            print("Source feature qualifiers:")
            for key, val in feature.qualifiers.items():
                print(f"  {key}: {val}")
        else:
            # Let's check if there are other interesting features like 'CDS' or 'misc_feature' containing notes
            notes = feature.qualifiers.get('note', [])
            for n in notes:
                if any(kw in n.lower() for kw in ['south korea', 'korea', 'univ', 'institute', 'lab', 'city', 'Seoul', 'Suwon', 'Jeonju']):
                    print(f"  [{feature.type}] note keyword match: {n}")
