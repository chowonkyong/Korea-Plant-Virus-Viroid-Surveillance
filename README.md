# National Meta-Analysis of Plant Viruses and Viroids in South Korea (1990–2026)

This repository contains the complete dataset, curation scripts, and multi-panel visualization code supporting our study published in *Physiological and Molecular Plant Pathology*.

## 1. Overview & Data Retrieval Parameters
To ensure full reproducibility and transparency, all surveillance records and associated metadata were harvested from the NCBI GenBank database.
* **Retrieval Date:** August 15, 2026
* **Database Sources:** NCBI Nucleotide and GenBank Flat File (GBFF) databases
* **Entrez Search Query:** 
  `("Viruses"[Organism] OR "Viroids"[Organism]) AND ("South Korea"[Country] OR "Republic of Korea"[Country] OR "Korea"[All Fields])`
* **Taxonomic Standard:** Official ICTV Master Species List (MSL39 / 2024–2025 release)

## 2. Quality Control & Curation Criteria (Exclusion Pipeline)
During the automated and manual metadata curation process, strict inclusion and exclusion criteria were applied:
* **Non-Plant Viruses:** Sequences belonging to non-plant viruses (e.g., mycoviruses, entomopathogenic viruses, and bacteriophages) inadvertently captured by broad host/country keywords were systematically excluded.
* **Reference Strains & Imported Isolates:** International reference strains lacking confirmed domestic field origin and intercepted quarantine isolates not established in South Korean agroecosystems were filtered out.
* **PRISMA-Style Stepwise Filtering:** 
  * **Initial Raw Extraction:** $n = 6,538$ total national records.
  * **Core Analytical Subset:** Entries containing incomplete, ambiguous, or missing metadata across mandatory fields (host species, genome type, vector transmission mode, and symptom phenotypes) ($n = 4,947$) were excluded, yielding a final high-confidence core subset of **$n = 1,591$ observations** for multi-variable network and Sankey analyses.

## 3. Repository Structure
* **`/scripts/`**: Custom Python v3.8+ scripts utilized for automated metadata extraction, data standardization, quality control, statistical analysis, and high-resolution (600 DPI) figure generation (including Figure 8 and Sankey diagrams).
* **`/data/`**: Curated CSV datasets and supplementary Excel files, including the refined core analytical dataset (`korea_plant_viruses_core_perfect.csv`) and the complete master metadata table (`Supplementary Table S1.xlsx`).

## 4. Code Validation & Transparency
While generative AI assisted in drafting boilerplate syntax for initial script development, all codes underwent rigorous human inspection, unit testing, and independent cross-validation against primary GenBank records (achieving a 96.5% metadata audit accuracy). All scripts and raw inputs are provided herein to enable complete independent verification.
