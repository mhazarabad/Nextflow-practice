nextflow.enable.dsl = 2

params.outdir = params.outdir ?: "results"
params.cache_dir = params.cache_dir ?: "cache"
params.bulk_matrix_url = params.bulk_matrix_url ?: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60424/matrix/GSE60424_series_matrix.txt.gz"
params.bulk_counts_url = params.bulk_counts_url ?: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60424/suppl/GSE60424_GEOSubmit_FC1to11_normalized_counts.txt.gz"
params.celltype = params.celltype ?: "Monocytes"
params.conditions = params.conditions ?: "Healthy Control,MS posttreatment"
params.peaks_url = params.peaks_url ?: "https://www.encodeproject.org/files/ENCFF002CDY/@@download/ENCFF002CDY.bed.gz"
params.max_chip_peaks = params.max_chip_peaks ?: 800
params.gtf_url = params.gtf_url ?: "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
params.scrna_cluster_key = params.scrna_cluster_key ?: "louvain"

process GET_BULK {
    publishDir "${params.outdir}/bulk", mode: 'copy'

    output:
        path "bulk_expression.tsv", emit: expression
        path "bulk_metadata.tsv", emit: metadata
        path "bulk_summary.json", emit: summary

    script:
    """
    python ${projectDir}/bin/parse_geo_matrix.py \\
      --matrix-url '${params.bulk_matrix_url}' \\
      --counts-url '${params.bulk_counts_url}' \\
      --celltype '${params.celltype}' \\
      --conditions '${params.conditions}' \\
      --out-expression bulk_expression.tsv \\
      --out-metadata bulk_metadata.tsv \\
      --summary-json bulk_summary.json \\
      --cache-dir ${params.cache_dir}
    """
}

process GET_GTF {
    publishDir "${params.outdir}/reference", mode: 'copy'

    output:
        path "gencode.v19.annotation.gtf.gz", emit: gtf

    script:
    """
    curl -L '${params.gtf_url}' -o gencode.v19.annotation.gtf.gz
    """
}

process GET_CHIP {
    publishDir "${params.outdir}/chipseq", mode: 'copy'

    output:
        path "stat1_peaks.bed", emit: peaks

    script:
    """
    curl -L '${params.peaks_url}' -o stat1_peaks.bed.gz
    gunzip -c stat1_peaks.bed.gz > stat1_peaks.bed
    rm stat1_peaks.bed.gz
    """
}

process BULK_DEG {
    publishDir "${params.outdir}/bulk", mode: 'copy'

input:
    path expression_file
    path metadata_file
    path gtf_file

    output:
        path "bulk_deg.tsv", emit: deg
        path "bulk_toplist.tsv", emit: toplist
        path "bulk_deg_summary.json", emit: summary

    script:
    """
    python ${projectDir}/bin/bulk_deg.py \\
      --expression ${expression_file} \\
      --metadata ${metadata_file} \\
      --gtf ${gtf_file} \\
      --out-deg bulk_deg.tsv \\
      --out-toplist bulk_toplist.tsv \\
      --summary-json bulk_deg_summary.json
    """
}

process ANNOTATE_CHIP {
    publishDir "${params.outdir}/chipseq", mode: 'copy'

input:
    path peak_file
    path gtf_file

    output:
        path "chip_peak_annotations.tsv", emit: annotations
        path "chip_gene_targets.tsv", emit: gene_summary
        path "chip_annotation_summary.json", emit: summary

    script:
    """
    python ${projectDir}/bin/annotate_chip.py \\
      --peaks ${peak_file} \\
      --gtf ${gtf_file} \\
      --max-peaks ${params.max_chip_peaks} \\
      --out-annotated chip_peak_annotations.tsv \\
      --out-genes chip_gene_targets.tsv \\
      --summary-json chip_annotation_summary.json
    """
}

process GET_SCRNA {
    publishDir "${params.outdir}/scrna", mode: 'copy'

    output:
        path "pbmc3k_processed.h5ad", emit: adata
        path "scrna_download_summary.json", emit: summary

    script:
    """
    python ${projectDir}/bin/download_scrna.py \\
      --out-file pbmc3k_processed.h5ad \\
      --summary-json scrna_download_summary.json
    """
}

process SCRNA_SIGNATURE {
    publishDir "${params.outdir}/scrna", mode: 'copy'

input:
    path adata_file
    path gene_table

    output:
        path "scrna_gene_signature.tsv", emit: signature
        path "scrna_signature_summary.json", emit: summary

    script:
    """
    python ${projectDir}/bin/scrna_signature.py \\
      --adata ${adata_file} \\
      --gene-table ${gene_table} \\
      --cluster-key '${params.scrna_cluster_key}' \\
      --out-signature scrna_gene_signature.tsv \\
      --summary-json scrna_signature_summary.json
    """
}

process BUILD_REPORT {
    publishDir "${params.outdir}/integration", mode: 'copy'

input:
    path deg_table
    path chip_gene_table
    path scrna_signature

    output:
        path "multiomic_summary.tsv", emit: summary_table
        path "multiomic_summary.json", emit: summary_json
    path "multiomic_heatmap.png", emit: heatmap
    path "chip_target_bar.png", emit: chip_bar

    script:
    """
    python ${projectDir}/bin/build_report.py \\
      --deg ${deg_table} \\
      --chip-genes ${chip_gene_table} \\
      --scrna-signature ${scrna_signature} \\
      --out-summary multiomic_summary.tsv \\
      --out-json multiomic_summary.json
    """
}

workflow {
    bulk = GET_BULK()
    gtf_reference = GET_GTF()
    chip_peaks = GET_CHIP()

    deg_results = BULK_DEG(bulk.expression, bulk.metadata, gtf_reference.gtf)
    chip_results = ANNOTATE_CHIP(chip_peaks.peaks, gtf_reference.gtf)
    scrna_results = GET_SCRNA()
    signature_results = SCRNA_SIGNATURE(scrna_results.adata, chip_results.gene_summary)

    BUILD_REPORT(deg_results.deg, chip_results.gene_summary, signature_results.signature)
}

