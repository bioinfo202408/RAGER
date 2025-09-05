# Load configuration
configfile: "config.yaml"

rule all:
    input:
        f"{config['paths']['dir']}/custom_genes_ensemblID.txt",
        f"{config['paths']['dir']}/custom_gene_gsea_GOenrich.csv",
        f"{config['paths']['dir']}/custom_gene_gsea_KEGGenrich.csv",
        f"{config['paths']['dir']}/promoter_seqs.fa",
        f"{config['paths']['dir']}/motif_enrich_result/analysis_complete.txt",
        f"{config['paths']['dir']}/TF_gene_network.txt",
        f"{config['paths']['dir']}/enriched_tfs_list.txt",
        f"{config['paths']['dir']}/TF_expr_heatmap.pdf"

rule convert_customGene_to_ensemlID:
    input:
        gene=f"{config['input_files']['custom_genes']}"
    output:
        ensembl=f"{config['paths']['dir']}/custom_genes_ensemblID.txt"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/convert_gene_to_ensemblID.R \
            {input.gene} \
            {output.ensembl}
        """

rule GSEA_enrichment:
    input:
        gene=f"{config['paths']['dir']}/custom_genes_ensemblID.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        GOenrich=f"{config['paths']['dir']}/custom_gene_gsea_GOenrich.csv",
        KEGGenrich=f"{config['paths']['dir']}/custom_gene_gsea_KEGGenrich.csv"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/GSEA_enrichment.R \
            {input.gene} \
            {input.DEgene} \
            {output.GOenrich} \
            {output.KEGGenrich}
        """

rule extr_promoterSeq:
    input:
        gene=f"{config['paths']['dir']}/custom_genes_ensemblID.txt"
    output:
        fa=f"{config['paths']['dir']}/promoter_seqs.fa"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/extr_promoterSeq.R \
            {input.gene} \
            {output.fa}
        """

rule motif_enrich:
    input:
        fa=f"{config['paths']['dir']}/promoter_seqs.fa",
        motif_dir=config["paths"]["jaspar_dir"]
    output:
        complete=f"{config['paths']['dir']}/motif_enrich_result/analysis_complete.txt",
        flag=touch("flags/motif_enrich_for_customGenes.done")
    params:
        output_dir=f"{config['paths']['dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/tf_motif_enrich_analysis.pl \
            --fastafile {input.fa} \
            --outputdir {params.output_dir} \
            --motifdir {input.motif_dir}
        
        echo "Analysis completed on $(date)" > {output.complete}
        """

rule TF_Gene_network:
    input:
        motif_analysis=f"{config['paths']['dir']}/motif_enrich_result/analysis_complete.txt"
    output:
        network=f"{config['paths']['dir']}/TF_gene_network.txt"
    params:
        motif_dir=f"{config['paths']['dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/build_TF_Gene_motif_network.pl \
            --motifdir {params.motif_dir} \
            --outfile {output.network}
        """

rule extract_enriched_motifs:
    input:
        flag="flags/motif_enrich_for_customGenes.done",
        tfmapfile=config["paths"]["tf_mapping_file"]
    output:
        enriched_tfs=f"{config['paths']['dir']}/enriched_tfs_list.txt"
    params:
        motif_dir=f"{config['paths']['dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/extr_enriched_motifs.pl \
            --tfmapfile {input.tfmapfile} \
            --motifdir {params.motif_dir} \
            --outfile {output.enriched_tfs}
        """

rule TF_heatmap:
    input:
        enriched_tfs=f"{config['paths']['dir']}/enriched_tfs_list.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        heatmap=f"{config['paths']['dir']}/TF_expr_heatmap.pdf"
    params:
        up_log2fc=config["params"]["log2fc_up"],
        up_padj=config["params"]["padj"],
        down_log2fc=config["params"]["log2fc_down"],
        down_padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/TFheatmap.R \
            {input.enriched_tfs} \
            {input.DEgene} \
            {params.up_log2fc} \
            {params.up_padj} \
            {params.down_log2fc} \
            {params.down_padj} \
            {output.heatmap}
        """