import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#333333'

fig, ax = plt.subplots(figsize=(11, 8.5), dpi=600)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')
fig.patch.set_facecolor('white')

def draw_box(x, y, w, h, title, body_lines, bg_color='#F4F8FB', border_color='#2B5C8F', title_color='#1A3A5C', lw=1.8):
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.6,rounding_size=1.5",
        facecolor=bg_color, edgecolor=border_color, linewidth=lw
    )
    ax.add_patch(rect)

    ax.text(
        x + w / 2, y + h - 2.6, title,
        ha='center', va='center', fontsize=10.5, fontweight='bold', color=title_color
    )

    ax.plot(
        [x + 1.5, x + w - 1.5], [y + h - 4.8, y + h - 4.8],
        color=border_color, lw=0.8, alpha=0.5
    )

    start_y = y + h - 7.5
    line_step = (h - 8.5) / max(len(body_lines), 1)
    for idx, line in enumerate(body_lines):
        is_bold = line.startswith("**")
        clean_line = line.replace("**", "")
        ax.text(
            x + w / 2, start_y - idx * line_step, clean_line,
            ha='center', va='center',
            fontsize=9.2 if not is_bold else 9.8,
            fontweight='bold' if is_bold else 'normal',
            color='#111111'
        )

def draw_arrow(x1, y1, x2, y2, color='#2B5C8F', lw=1.8):
    ax.annotate(
        '', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle='->', color=color, lw=lw, shrinkA=0, shrinkB=0, mutation_scale=15)
    )

# Title header
ax.text(
    50, 96.5,
    "PRISMA-Style Data Curation and Filtering Workflow",
    ha='center', va='center', fontsize=14, fontweight='bold', color='#111111'
)
ax.text(
    50, 93.2,
    "National Surveillance Inventory of Plant Viruses and Viroids in South Korea (1990–2026)",
    ha='center', va='center', fontsize=10, style='italic', color='#444444'
)

# Box 1: Identification
draw_box(
    x=14, y=76, w=48, h=13,
    title="1. Identification (NCBI GenBank Retrieval)",
    body_lines=[
        '("Viruses"[Organism] OR "Viroids"[Organism]) & Korea filters',
        "**Initial Records Retrieved: n = 7,393**"
    ],
    bg_color='#EEF4FA', border_color='#2B5C8F', title_color='#163A5F'
)

# Exclusion Box 1
draw_box(
    x=62, y=55, w=34, h=14,
    title="Excluded: Step 1 QC (n = 855)",
    body_lines=[
        "• Non-plant viruses (animal/insect/fungal)",
        "• Reference strains / imported isolates",
        "• Duplicates & out-of-scope entries"
    ],
    bg_color='#FDF2F2', border_color='#C0392B', title_color='#922B21', lw=1.4
)

# Box 2: Curated National Inventory
draw_box(
    x=14, y=54, w=40, h=13.5,
    title="2. Curated National Inventory",
    body_lines=[
        "ICTV-standardized plant virus & viroid records",
        "**Included National Inventory: n = 6,538**"
    ],
    bg_color='#EBF7F2', border_color='#1E8449', title_color='#145A32'
)

# Exclusion Box 2
draw_box(
    x=62, y=31, w=34, h=14.5,
    title="Excluded: Incomplete Metadata (n = 4,947)",
    body_lines=[
        "Missing or ambiguous entries across ≥1 of 4",
        "mandatory fields: host, genome type,",
        "vector transmission mode, or symptom"
    ],
    bg_color='#FDF2F2', border_color='#C0392B', title_color='#922B21', lw=1.4
)

# Box 3: Core Analytical Subset
draw_box(
    x=14, y=29, w=40, h=14.5,
    title="3. Core Analytical Subset (Sankey & Multi-Var)",
    body_lines=[
        "Complete high-confidence metadata across all",
        "four mandatory dimensions",
        "**Final Core Observations: n = 1,591**"
    ],
    bg_color='#FEF9E7', border_color='#B7950B', title_color='#7D6608'
)

# Bottom Transparency Box
draw_box(
    x=14, y=6, w=82, h=15,
    title="Data Transparency & Reproducibility Notes",
    body_lines=[
        "• Retrieval date: August 15, 2026 | Taxonomy standard: ICTV Master Species List (MSL39)",
        "• National inventory (n = 6,538) supports taxonomic, temporal, and genome/vector analyses.",
        "• Core subset (n = 1,591) provides complete multi-dimensional records for Sankey & association modeling."
    ],
    bg_color='#F8F9F9', border_color='#5D6D7E', title_color='#2C3E50', lw=1.2
)

# Arrows
draw_arrow(34, 75.3, 34, 68.2, color='#2B5C8F')
ax.plot([34, 57], [71.8, 71.8], color='#C0392B', lw=1.4)
ax.plot([57, 57], [71.8, 62], color='#C0392B', lw=1.4)
draw_arrow(57, 62, 61.3, 62, color='#C0392B', lw=1.4)

draw_arrow(34, 53.3, 34, 44.2, color='#1E8449')
ax.plot([34, 57], [48.8, 48.8], color='#C0392B', lw=1.4)
ax.plot([57, 57], [48.8, 38.2], color='#C0392B', lw=1.4)
draw_arrow(57, 38.2, 61.3, 38.2, color='#C0392B', lw=1.4)

draw_arrow(34, 28.3, 34, 21.7, color='#5D6D7E')

plt.tight_layout()
plt.savefig('prisma_filtering_flow.png', dpi=600, bbox_inches='tight')
plt.savefig('prisma_filtering_flow.pdf', bbox_inches='tight')
plt.savefig('prisma_filtering_flow.svg', bbox_inches='tight')
print("Saved prisma_filtering_flow.png, .pdf, and .svg successfully.")
