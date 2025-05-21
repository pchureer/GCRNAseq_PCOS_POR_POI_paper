## GCRNAseq_PCOS_POR_POI_paper
## This repository contains scripts to reproduce the analyses from "Granulosa Cell RNA-seq Reveals Overlapping and Unique Differential Expression Signatures in PCOS, POR, and POI Using Pooled vs. Matched Controls"

### Prerequisites

### 1. Clone the repo
git clone https://github.com/yourusername/GCRNAseq_PCOS_POR_POI_paper.git

cd GCRNAseq_PCOS_POR_POI_paper

### 2. Install Conda environment
conda env create -f environment.yml

conda activate gc_rnaseq

### 3. Directory layout

    ├── data/                  # Reference files
    │   ├── adapters.fa
    │   ├── OpenTargets_PCOS.csv
    │   ├── OpenTargets_OvarianDysfunction.csv
    │   └── OpenTargets_POI.csv
    ├── results/               # Outputs of each step
    ├── figures/               # Plots generated
    └── scripts/
        ├── 01_preprocessing/
        ├── 02_dge/
        ├── 03_venn/
        ├── 04_heatmap/
        ├── 05_intersection/
        └── 06_enrichment/
    
### Script templates

### 4.1. Trim adapters & align reads
**Script**: scripts/01_preprocessing/01_trim_and_align.sh

**Purpose**: Remove adapters, filter low-quality bases, and align to GRCh38 with STAR, producing sorted BAM and gene counts.
    
    bash scripts/01_preprocessing/01_trim_and_align.sh \
        SRR10239199 \
        data/fastq/SRR10239199_R1.fq.gz \
        data/fastq/SRR10239199_R2.fq.gz

    Args:
        Sample ID (e.g. SRR10239199)
        R1 FASTQ path
        R2 FASTQ path

    Outputs (in results/SRR10239199/):
        SRR10239199_R1_paired.fq.gz and _R2_paired.fq.gz
        SRR10239199.Aligned.sortedByCoord.out.bam
        ReadsPerGene.out.tab (raw gene counts)

### 4.2. Differential expression with edgeR
    - Script: scripts/02_dge/02_edgeR_DE.R
    - Purpose: Normalize counts, model batch + group, and identify DE genes.

    - Rscript scripts/02_dge/02_edgeR_DE.R \
          results/counts_matrix.csv \
          results/metadata.csv \
          results/DE_results_PCOS_pooled.csv

    - Inputs:
        counts_matrix.csv: genes × samples raw counts (from all ReadsPerGene.out.tab files).
        metadata.csv: sample metadata (sampleID, group, batch).

    - Output:
        DE_results_*.csv: table with gene, log2FC, PValue, FDR for each contrast.

### 4.3. Venn diagrams of DE gene intersections
    - Script: scripts/03_venn/03_plot_venn.R
    - Purpose: Show overlap of significant genes across phenotypes.

    - Rscript scripts/03_venn/03_plot_venn.R \
        results/DE_results_PCOS_pooled.csv \
        results/DE_results_POR_pooled.csv \
        results/DE_results_POI_pooled.csv \
        figures/venn_pooled.png

    - Inputs: three DE result CSVs (pooled-control contrasts).
    - Output: venn_pooled.png illustrating shared/unique genes.

### 4.4. Heatmap of core 17-gene signature
    - Script: scripts/04_heatmap/04_plot_heatmap.R
    - Purpose: Visualize log₂FC of the 17 shared signature genes across all six contrasts.

    - Rscript scripts/04_heatmap/04_plot_heatmap.R \
        results/log2FC_17core.csv \
        figures/core17_heatmap.png

    - Input: log2FC_17core.csv (rows = genes, columns = contrasts).
    - Output: core17_heatmap.png.

### 4.5. Intersection-count & Jaccard similarity heatmaps
    - Script: scripts/05_intersection/05_intersection_heatmap_jaccard.R
    - Purpose: Compare your DE sets (six contrasts) against Open Targets reference lists.

    - Rscript scripts/05_intersection/05_intersection_heatmap_jaccard.R \
        results/DE_results_PCOS_pooled.csv \
        results/DE_results_PCOS_matched.csv \
        results/DE_results_POR_pooled.csv \
        results/DE_results_POR_matched.csv \
        results/DE_results_POI_pooled.csv \
        results/DE_results_POI_matched.csv \
        data/OpenTargets_PCOS.csv \
        data/OpenTargets_OvarianDysfunction.csv \
        data/OpenTargets_POI.csv \
        figures/intersect_count_heatmap.png \
        figures/jaccard_heatmap.png

    - Outputs:
        intersect_count_heatmap.png: absolute overlaps.
        jaccard_heatmap.png: normalized similarities.

### 4.6. Text enrichment (GO terms word-cloud)
    - Script: scripts/06_enrichment/06_text_enrichment.R
    - Purpose: Generate word-clouds of significantly enriched GO Biological Process terms.

    - Rscript scripts/06_enrichment/06_text_enrichment.R \
        results/GO_enrichment_PCOS_pooled.csv \
        figures/PCOS_pooled_wordcloud.png

    - Input: GO enrichment CSV (term, Count, FDR).
    - Output: PCOS_pooled_wordcloud.png.

## Cite us: 

## License
    This project is licensed under the MIT License.
