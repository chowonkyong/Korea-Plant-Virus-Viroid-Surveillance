import pandas as pd
import numpy as np
import plotly.graph_objects as go

# 1. 데이터 로드 및 'Unknown' 숙주 제거
try:
    df = pd.read_csv('korea_plant_viruses_core_perfect.csv')
    if 'Standardized_Host' in df.columns:
        df = df[df['Standardized_Host'].str.lower() != 'unknown']
except Exception as e:
    print(f"Dataset load error: {e}")

# 2. 기주 이름 정제 (Chrysanthemum sp. -> Chrysanthemum morifolium)
def clean_host_name(host):
    name_str = str(host)
    if 'Chrysanthemum' in name_str:
        return 'Chrysanthemum morifolium'
    return name_str

# 바이러스 과(Family) 매핑 함수
def get_virus_family(name):
    name_str = str(name).lower()
    if 'viroid' in name_str:
        if 'latent mosaic' in name_str or 'plmvd' in name_str:
            return 'Avsunviroidae (Viroid)'
        return 'Pospiviroidae (Viroid)'
    elif 'tswv' in name_str or 'tospovirus' in name_str or 'orthotospovirus' in name_str:
        return 'Tospoviridae'
    elif 'cmv' in name_str or 'cucumber mosaic' in name_str or 'cucumovirus' in name_str:
        return 'Bromoviridae'
    elif 'smv' in name_str or 'wmv' in name_str or 'lysv' in name_str or 'mosaic' in name_str or 'potyvirus' in name_str:
        return 'Potyviridae'
    elif 'tylcv' in name_str or 'yellow leaf curl' in name_str or 'gemini' in name_str:
        return 'Geminiviridae'
    elif 'capillo' in name_str or 'betaflexi' in name_str:
        return 'Betaflexiviridae'
    return 'Other Plant Viruses'

df['Clean_Host'] = df['Standardized_Host'].apply(clean_host_name)
df['Clean_Vector'] = df['Vector_Mode']
df['Virus_Family'] = df['Organism'].apply(get_virus_family)

# 상위 기주 및 벡터 필터링
top_n_hosts = df['Clean_Host'].value_counts().head(15).index
top_n_vectors = df['Clean_Vector'].value_counts().head(6).index

df_filtered = df[
    df['Clean_Host'].isin(top_n_hosts) & 
    df['Clean_Vector'].isin(top_n_vectors)
].copy()

df_filtered['Clean_Virus'] = df_filtered['Virus_Family']

# 벡터 모드별 빈도 및 비율(%) 계산 (df_filtered 기준)
vector_total_n = len(df_filtered)
vector_counts = df_filtered['Clean_Vector'].value_counts().to_dict()
vector_pcts = {
    vec: (cnt / vector_total_n) * 100 if vector_total_n > 0 else 0.0
    for vec, cnt in vector_counts.items()
}

print(f"Total filtered records in Sankey (n): {vector_total_n}")
for vec, cnt in vector_counts.items():
    print(f"  {vec}: n = {cnt} ({vector_pcts[vec]:.1f}%)")

# 노드 정의 및 3단 레이어 좌표 설정 (좌: 0.01, 중: 0.50, 우: 0.99)
hosts = sorted(df_filtered['Clean_Host'].unique())
vectors = sorted(df_filtered['Clean_Vector'].unique())
viruses = sorted(df_filtered['Clean_Virus'].unique())

node_raw_labels = list(hosts) + list(vectors) + list(viruses)

node_labels = []
for label in list(hosts):
    node_labels.append(f"<b>{label}</b>")
for label in list(vectors):
    cnt = vector_counts.get(label, 0)
    pct = vector_pcts.get(label, 0.0)
    node_labels.append(f"<b>{label}</b> ({cnt}, {pct:.1f}%)")
for label in list(viruses):
    node_labels.append(f"<b>{label}</b>")

node_dict = {label: i for i, label in enumerate(node_raw_labels)}

node_x = [0.01] * len(hosts) + [0.50] * len(vectors) + [0.99] * len(viruses)
node_y = [None] * len(node_labels)

flow1 = df_filtered.groupby(['Clean_Host', 'Clean_Vector']).size().reset_index(name='Count')
flow2 = df_filtered.groupby(['Clean_Vector', 'Clean_Virus']).size().reset_index(name='Count')

source_indices, target_indices, values, link_colors = [], [], [], []
host_palette = ['#457b9d', '#1d3557', '#2a9d8f', '#264653', '#3a86ff', '#00b4d8', '#0077b6', '#17c3b2', '#e07a5f', '#e63946', '#f4a261', '#2a9d8f', '#264653', '#457b9d', '#1d3557']
vector_palette = ['#e76f51', '#f4a261', '#f77f00', '#fcbf49', '#d62828', '#b5179e']
virus_palette = ['#2b9348', '#55a630', '#80b918', '#aacc00', '#bfd200', '#38b000', '#70e000', '#3a86ff']

node_colors = (
    [host_palette[i % len(host_palette)] for i in range(len(hosts))] + 
    [vector_palette[i % len(vector_palette)] for i in range(len(vectors))] + 
    [virus_palette[i % len(virus_palette)] for i in range(len(viruses))]
)

for _, row in flow1.iterrows():
    src_idx = node_dict[row['Clean_Host']]
    col_hex = host_palette[src_idx % len(host_palette)]
    source_indices.append(src_idx)
    target_indices.append(node_dict[row['Clean_Vector']])
    values.append(row['Count'])
    r, g, b = int(col_hex[1:3], 16), int(col_hex[3:5], 16), int(col_hex[5:7], 16)
    link_colors.append(f"rgba({r}, {g}, {b}, 0.3)")

for _, row in flow2.iterrows():
    src_idx = node_dict[row['Clean_Vector']]
    col_hex = vector_palette[(src_idx - len(hosts)) % len(vector_palette)]
    source_indices.append(src_idx)
    target_indices.append(node_dict[row['Clean_Virus']])
    values.append(row['Count'])
    r, g, b = int(col_hex[1:3], 16), int(col_hex[3:5], 16), int(col_hex[5:7], 16)
    link_colors.append(f"rgba({r}, {g}, {b}, 0.3)")

fig = go.Figure(data=[go.Sankey(
    arrangement = "perpendicular",
    node = dict(
        pad = 15, 
        thickness = 20, 
        line = dict(color = "#ffffff", width = 1.0), 
        label = node_labels, 
        color = node_colors, 
        x = node_x, 
        y = node_y,
        hovertemplate='%{label}<br>Total flow: %{value}<extra></extra>'
    ),
    link = dict(
        source = source_indices, 
        target = target_indices, 
        value = values, 
        color = link_colors,
        hovertemplate='%{source.label} → %{target.label}<br>Records: %{value}<extra></extra>'
    )
)])

# 타이틀 제거 및 논문용 여백 최적화
fig.update_layout(
    font=dict(size=13, color="#111111", family="Arial Black, Arial, sans-serif"),
    plot_bgcolor='#ffffff', 
    paper_bgcolor='#ffffff', 
    width=1650, 
    height=1000,
    margin=dict(t=30, b=30, l=30, r=30)
)

# HTML 저장
fig.write_html("Figure_3_Sankey.html")
print("Successfully generated HTML: Figure_3_Sankey.html")

# 벡터 포맷(SVG, PDF) 및 고해상도 이미지 저장 (kaleido 필요)
try:
    fig.write_image("Figure_3_Sankey.svg")
    fig.write_image("Figure_3_Sankey.pdf")
    print("Successfully generated vector formats: Figure_3_Sankey.svg, Figure_3_Sankey.pdf")
except Exception as e:
    print(f"Vector export note (requires kaleido): {e}")