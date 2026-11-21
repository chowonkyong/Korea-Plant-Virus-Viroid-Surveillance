import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Verified core plant virus timeline data up to 2026 (Polydnavirus citrulli and unverified taxons excluded)
timeline_data = [
    (1990, 'Cucumber mosaic virus'),
    (1994, 'Potato virus Y'),
    (1995, 'Cymbidium mosaic virus'),
    (1996, 'Turnip mosaic virus'),
    (1996, 'Garlic latent virus'),
    (1997, 'Tobacco mosaic virus'),
    (1997, 'Potato virus S'),
    (1999, 'Hop stunt viroid'),
    (1999, 'Cucumber green mottle mosaic virus'),
    (2000, 'Broad bean wilt virus 2'),
    (2000, 'Alfalfa mosaic virus'),
    (2001, 'Chrysanthemum stunt viroid'),
    (2001, 'Tomato spotted wilt orthotospovirus'),
    (2002, 'Soybean mosaic virus'),
    (2002, 'Garlic virus A'),
    (2003, 'Soybean yellow mottle mosaic virus'),
    (2003, 'Barley mild mosaic virus'),
    (2004, 'Bean common mosaic virus'),
    (2004, 'Pear black necrotic leaf spot virus'),
    (2005, 'Peanut stunt virus'),
    (2005, 'Daphne virus S'),
    (2006, 'Garlic common latent virus'),
    (2006, 'Zucchini yellow mosaic virus'),
    (2007, 'Lily mottle virus'),
    (2007, 'Youcai mosaic virus'),
    (2008, 'Tomato yellow leaf curl virus'),
    (2008, 'Capillovirus mali'),
    (2009, 'Watermelon mosaic virus'),
    (2009, 'Trichovirus mali'),
    (2010, 'Plantago asiatica mosaic virus'),
    (2010, 'Citrus tristeza virus'),
    (2011, 'Pepper cryptic virus 2'),
    (2011, 'Leek yellow stripe virus'),
    (2012, 'Peanut mottle virus'),
    (2012, 'Marafivirus mali'),
    (2013, 'Cherry green ring mottle virus'),
    (2013, 'Tomato chlorosis virus'),
    (2014, 'Impatiens necrotic spot virus'),
    (2014, 'Cucurbit chlorotic yellows virus'),
    (2015, 'Shallot latent virus'),
    (2015, 'Papaya leaf curl Guandong virus'),
    (2016, 'Onion yellow dwarf virus'),
    (2016, 'Potato virus X'),
    (2017, 'Citrus leaf blotch virus'),
    (2017, 'Tomato leaf curl New Delhi virus'),
    (2018, 'Citrus viroid I-LSS'),
    (2019, 'Chinese artichoke mosaic virus'),
    (2019, 'Nerine latent virus'),
    (2020, 'Pepino mosaic virus'),
    (2020, 'Narcissus late season yellows virus'),
    (2021, 'Blackberry chlorotic ringspot virus'),
    (2021, 'Arabis mosaic virus'),
    (2022, 'Citrus exocortis viroid'),
    (2023, 'Strawberry necrotic shock virus'),
    (2023, 'Grapevine virus E'),
    (2024, 'Avocado sunblotch viroid'),
    (2025, 'Blueberry shock virus'),
    (2025, 'Tomato brown rugose fruit virus'),
    (2026, 'Cucumber fruit mottle mosaic virus')
]

# 캔버스 가로/세로 크기 확장 및 600 DPI 설정
fig, ax = plt.subplots(figsize=(92, 48), dpi=600)
fig.patch.set_facecolor('#ffffff')

# 우측 끝 긴 학명('2026 | Cucumber fruit mottle mosaic virus') 여백 확보
min_year, max_year = 1989, 2038
ax.set_xlim(min_year, max_year)
ax.set_ylim(-8.5, 8.5)

# 시대별 배경색 적용
ax.axvspan(1989, 1999.5, ymin=0, ymax=1, color='#eff6ff', alpha=0.6, zorder=0)
ax.axvspan(1999.5, 2009.5, ymin=0, ymax=1, color='#fef3c7', alpha=0.6, zorder=0)
ax.axvspan(2009.5, 2019.5, ymin=0, ymax=1, color='#ecfdf5', alpha=0.6, zorder=0)
ax.axvspan(2019.5, max_year, ymin=0, ymax=1, color='#fef2f2', alpha=0.6, zorder=0)

# 5년 단위 세로 그리드 선
for yr in range(1990, 2027, 5):
    ax.axvline(x=yr, color='#cbd5e1', linestyle=':', linewidth=1.2, alpha=0.7, zorder=1)

# 중심 기준선
ax.hlines(y=0, xmin=min_year, xmax=max_year, color='#1e293b', linewidth=6, zorder=2)

# 위/아래 각각 22단계 트랙 풀 설정
top_tiers = np.linspace(0.65, 8.00, 22)
bot_tiers = np.linspace(-0.65, -8.00, 22)

# 34pt 텍스트 가로 길이(연도 단위)를 넉넉하게 추정
def estimate_x_span(yr, virus_name):
    label_str = f"{int(yr)}  |  {virus_name}"
    x_start = yr + 0.15
    x_end = x_start + len(label_str) * 0.235
    return x_start, x_end

top_occupied = {i: [] for i in range(len(top_tiers))}
bot_occupied = {i: [] for i in range(len(bot_tiers))}
level_map = {}

golden_ratio = 0.618033988749895

for idx, (yr, virus) in enumerate(timeline_data):
    x0, x1 = estimate_x_span(yr, virus)
    pad = 0.6

    if idx % 2 == 0:
        n_tiers = len(top_tiers)
        base_idx = int(((idx * golden_ratio) % 1) * n_tiers)
        candidate_indices = [(base_idx + k * 5) % n_tiers for k in range(n_tiers)]
        chosen_idx = None
        for t_i in candidate_indices:
            if all((x1 + pad < a) or (x0 - pad > b) for (a, b) in top_occupied[t_i]):
                chosen_idx = t_i
                break
        if chosen_idx is None:
            best_t = 0
            best_ov = float('inf')
            for t_i in range(n_tiers):
                ov = sum(max(0.0, min(x1 + pad, b) - max(x0 - pad, a)) for (a, b) in top_occupied[t_i])
                if ov < best_ov:
                    best_ov = ov
                    best_t = t_i
            chosen_idx = best_t
        top_occupied[chosen_idx].append((x0, x1))
        level_map[(yr, virus)] = top_tiers[chosen_idx]
    else:
        n_tiers = len(bot_tiers)
        base_idx = int(((idx * golden_ratio) % 1) * n_tiers)
        candidate_indices = [(base_idx + k * 5) % n_tiers for k in range(n_tiers)]
        chosen_idx = None
        for t_i in candidate_indices:
            if all((x1 + pad < a) or (x0 - pad > b) for (a, b) in bot_occupied[t_i]):
                chosen_idx = t_i
                break
        if chosen_idx is None:
            best_t = 0
            best_ov = float('inf')
            for t_i in range(n_tiers):
                ov = sum(max(0.0, min(x1 + pad, b) - max(x0 - pad, a)) for (a, b) in bot_occupied[t_i])
                if ov < best_ov:
                    best_ov = ov
                    best_t = t_i
            chosen_idx = best_t
        bot_occupied[chosen_idx].append((x0, x1))
        level_map[(yr, virus)] = bot_tiers[chosen_idx]

for yr, virus in timeline_data:
    level_y = level_map[(yr, virus)]

    if yr < 2000:
        pt_color = '#2563eb'
    elif yr < 2010:
        pt_color = '#d97706'
    elif yr < 2020:
        pt_color = '#059669'
    else:
        pt_color = '#dc2626'

    ax.scatter(yr, 0, color=pt_color, s=420, zorder=4, edgecolor='white', linewidth=3)
    ax.vlines(x=yr, ymin=0, ymax=level_y, color='#64748b', linestyle='--', linewidth=1.6, zorder=3)

    ax.text(
        yr + 0.15,
        level_y,
        f"{int(yr)}  |  {virus}",
        rotation=0,
        ha='left',
        va='center',
        fontsize=34,
        fontweight='heavy',
        color='#000000'
    )

# 범례(Legend) 배치
legend_elements = [
    Patch(facecolor='#2563eb', edgecolor='white', label='1990s Era'),
    Patch(facecolor='#d97706', edgecolor='white', label='2000s Era'),
    Patch(facecolor='#059669', edgecolor='white', label='2010s Era'),
    Patch(facecolor='#dc2626', edgecolor='white', label='2020s Era')
]
ax.legend(
    handles=legend_elements,
    loc='upper left',
    bbox_to_anchor=(0.01, 0.98),
    fontsize=32,
    frameon=True,
    facecolor='#ffffff',
    edgecolor='#cbd5e1',
    framealpha=0.95,
    borderpad=1.3
)

ax.axis('off')
plt.tight_layout()

# 600 DPI 고해상도로 PNG, SVG, PDF 포맷 동시 저장
base_filename = "Figure_5_Timeline_Perfect_Final"
for ext in ['png', 'svg', 'pdf']:
    output_filename = f"{base_filename}.{ext}"
    plt.savefig(output_filename, dpi=600, bbox_inches='tight', facecolor=fig.get_facecolor(), edgecolor='none')
    print(f"Successfully generated {output_filename} (600 DPI)")

plt.close(fig)