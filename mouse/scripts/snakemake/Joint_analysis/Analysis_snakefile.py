# Load configuration
configfile: "config.yaml"

# Extract samples from config
SAMPLES = list(config["samples"].keys())

rule all:
    input:
        f"{config['paths']['promoter_analysis_dir']}/UpVenn.pdf",
        f"{config['paths']['promoter_analysis_dir']}/shared_up_genes_INpromote.txt",
        f"{config['paths']['promoter_analysis_dir']}/DownVenn.pdf",
        f"{config['paths']['promoter_analysis_dir']}/shared_down_genes_INpromote.txt",
        f"{config['paths']['promoter_analysis_dir']}/shared_genes_INpromote_volcano_plot.pdf",
        f"{config['paths']['promoter_analysis_dir']}/shared_up_peak_INpromote.bed",
        f"{config['paths']['promoter_analysis_dir']}/shared_down_peak_INpromote.bed",
        f"{config['paths']['promoter_analysis_dir']}/shared_peak_INpromote.bed",
        f"{config['paths']['promoter_analysis_dir']}/shared_peak_scaled.tab",
        f"{config['paths']['promoter_analysis_dir']}/Expr_ccess_corr_plot.pdf",
        f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_GOenrich.csv",
        f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_GOenrich.rds",
        f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_KEGGenrich.csv",
        f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_KEGGenrich.rds",
        f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_GOenrich.csv",
        f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_GOenrich.rds",
        f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_KEGGenrich.csv",
        f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_KEGGenrich.rds",
        f"{config['paths']['jaspar_dir']}/Alx1.meme",
        f"{config['params']['regulation_analysis']['promoter_up_dir']}/promoter_seqs.fa",
        f"{config['params']['regulation_analysis']['promoter_up_dir']}/motif_enrich_result/analysis_complete.txt",
        f"{config['params']['regulation_analysis']['promoter_up_dir']}/SharedUPgene_in_promoter_TF_gene_network.txt",
        f"{config['params']['regulation_analysis']['promoter_up_dir']}/enriched_tfs_list.txt",
        f"{config['params']['regulation_analysis']['promoter_down_dir']}/promoter_seqs.fa",
        f"{config['params']['regulation_analysis']['promoter_down_dir']}/motif_enrich_result/analysis_complete.txt",
        f"{config['params']['regulation_analysis']['promoter_down_dir']}/SharedDOWNgene_in_promoter_TF_gene_network.txt",
        f"{config['params']['regulation_analysis']['promoter_down_dir']}/enriched_tfs_list.txt",
        f"{config['params']['regulation_analysis']['promoter_up_dir']}/shared_UpGene_InPromote_tfs_heatmap_RNAseqExpr.pdf",
        f"{config['params']['regulation_analysis']['promoter_down_dir']}/shared_DownGene_InPromote_tfs_heatmap_RNAseqExpr.pdf",
        f"{config['paths']['enhancer_analysis_dir']}/UpGene_anno.bed",
        f"{config['paths']['enhancer_analysis_dir']}/DownGene_anno.bed",
        f"{config['paths']['enhancer_analysis_dir']}/enh_UpGene_within_1mb.txt",
        f"{config['paths']['enhancer_analysis_dir']}/enh_DownGene_within_1mb.txt",
        f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_UpGene.bed",
        f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_DownGene.bed",
        f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_UPgene.gz",
        f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_UPgene.tab",
        f"{config['paths']['enhancer_analysis_dir']}/enh_sortedRegions_for_UPgene.bed",
        f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_DOWNgene.gz",
        f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_DOWNgene.tab",
        f"{config['paths']['enhancer_analysis_dir']}/enh_sortedRegions_for_DOWNgene.bed",
        f"{config['paths']['enhancer_analysis_dir']}/SharedUPenhancer.txt",
        f"{config['paths']['enhancer_analysis_dir']}/SharedUPgene_INenhancer.txt",
        f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNenhancer.txt",
        f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNgene_INenhancer.txt",
        f"{config['paths']['enhancer_analysis_dir']}/shared_genes_INenhancer_volcano_plot.pdf",
        f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_GOenrich.csv",
        f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_GOenrich.rds",
        f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_KEGGenrich.csv",
        f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_KEGGenrich.rds",
        f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_GOenrich.csv",
        f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_GOenrich.rds",
        f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_KEGGenrich.csv",
        f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_KEGGenrich.rds",
        f"{config['params']['regulation_analysis']['enhancer_up_dir']}/SharedUP_enhancer_sequences.fa",
        f"{config['params']['regulation_analysis']['enhancer_down_dir']}/SharedDOWN_enhancer_sequences.fa",
        f"{config['params']['regulation_analysis']['enhancer_up_dir']}/motif_enrich_result/analysis_complete.txt",
        f"{config['params']['regulation_analysis']['enhancer_up_dir']}/SharedUPenhancer_TF_gene_network.txt",
        f"{config['params']['regulation_analysis']['enhancer_up_dir']}/enriched_tfs_list.txt",
        f"{config['params']['regulation_analysis']['enhancer_down_dir']}/motif_enrich_result/analysis_complete.txt",
        f"{config['params']['regulation_analysis']['enhancer_down_dir']}/SharedDOWNenhancer_TF_gene_network.txt",
        f"{config['params']['regulation_analysis']['enhancer_down_dir']}/enriched_tfs_list.txt",
        f"{config['params']['regulation_analysis']['enhancer_up_dir']}/shared_UpEnhancer_tfs_heatmap_RNAseqExpr.pdf",
        f"{config['params']['regulation_analysis']['enhancer_down_dir']}/shared_DownEnhancer_tfs_heatmap_RNAseqExpr.pdf"

rule shared_up_gene_promoter:
    input:
        UPpeak=f"{config['input_files']['up_peak_csv']}",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        Venn=f"{config['paths']['promoter_analysis_dir']}/UpVenn.pdf",
        gene=f"{config['paths']['promoter_analysis_dir']}/shared_up_genes_INpromote.txt"
    params:
        Log2FC=config["params"]["log2fc_up"],
        padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/identificate_shared_UpGenes_InPromoter.R {input.UPpeak} {input.DEgene} {params.Log2FC} {params.padj} {output.Venn} {output.gene}
        """

rule shared_down_gene_promoter:
    input:
        DOWNpeak=f"{config['input_files']['down_peak_csv']}",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        Venn=f"{config['paths']['promoter_analysis_dir']}/DownVenn.pdf",
        gene=f"{config['paths']['promoter_analysis_dir']}/shared_down_genes_INpromote.txt"
    params:
        Log2FC=config["params"]["log2fc_down"],
        padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/identificate_shared_DownGenes_InPromoter.R {input.DOWNpeak} {input.DEgene} {params.Log2FC} {params.padj} {output.Venn} {output.gene}
        """

rule shared_genes_INpromote_volcano_plot:
    input:
        UPgene=f"{config['paths']['promoter_analysis_dir']}/shared_up_genes_INpromote.txt",
        DOWNgene=f"{config['paths']['promoter_analysis_dir']}/shared_down_genes_INpromote.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        volcano_plot=f"{config['paths']['promoter_analysis_dir']}/shared_genes_INpromote_volcano_plot.pdf"
    params:
        Log2FC=config["params"]["log2fc_up"],
        padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/shared_gene_volcano_plot.R \
            {input.DEgene} \
            {input.UPgene} \
            {input.DOWNgene} \
            {params.Log2FC} \
            {params.padj} \
            {output.volcano_plot}
        """

rule extract_shared_up_peak_anno:
    input:
        UPpeak=f"{config['input_files']['up_peak_csv']}",
        DEgene=f"{config['input_files']['deg_csv']}",
        Allpeak_anno=f"{config['input_files']['up_peak_xls']}"
    output:
        peakanno=f"{config['paths']['promoter_analysis_dir']}/shared_up_peak_INpromote.bed"
    params:
        Log2FC=config["params"]["log2fc_up"],
        padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/expr_shared_up_peak_anno.R {input.UPpeak} {input.DEgene} {params.Log2FC} {params.padj} {input.Allpeak_anno} {output.peakanno}
        """

rule expr_shared_peak_anno:
    input:
        DOWNpeak=f"{config['input_files']['down_peak_csv']}",
        DEgene=f"{config['input_files']['deg_csv']}",
        Allpeak_anno=f"{config['input_files']['down_peak_xls']}",
        UPpeakanno=f"{config['paths']['promoter_analysis_dir']}/shared_up_peak_INpromote.bed"
    output:
        DOWNpeakanno=f"{config['paths']['promoter_analysis_dir']}/shared_down_peak_INpromote.bed",
        peakanno=f"{config['paths']['promoter_analysis_dir']}/shared_peak_INpromote.bed" 
    params:
        Log2FC=config["params"]["log2fc_down"],
        padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/expr_sharedpeak_anno.R {input.DOWNpeak} {input.DEgene} {params.Log2FC} {params.padj} {input.Allpeak_anno} {output.DOWNpeakanno} {input.UPpeakanno} {output.peakanno}
        """

rule calATACexpr:
    input:
        bw=expand(f"{config['paths']['atacseq_bowtie2_dir']}/{{sample}}.bw", sample=SAMPLES),
        bedfile=f"{config['paths']['promoter_analysis_dir']}/shared_peak_INpromote.bed" 
    output:
        gz=f"{config['paths']['promoter_analysis_dir']}/shared_peak_scaled.gz",
        tab=f"{config['paths']['promoter_analysis_dir']}/shared_peak_scaled.tab"
    params:
        before=config["params"]["computematrix"]["before"],
        after=config["params"]["computematrix"]["after"],
        bin_size=config["params"]["computematrix"]["bin_size"],
        processors=config["params"]["computematrix"]["processors"]
    run:
        bw_files = " ".join(input.bw)
        shell_cmd = f"computeMatrix reference-point -R {input.bedfile} -S {bw_files} -b {params.before} -a {params.after} --binSize {params.bin_size} -p {params.processors} -o {output.gz} --outFileNameMatrix {output.tab}"
        shell(shell_cmd)

rule plot_correlation:
    input:
        DEgene=f"{config['input_files']['deg_csv']}",
        peakanno=f"{config['paths']['promoter_analysis_dir']}/shared_peak_INpromote.bed",
        Peaktab=f"{config['paths']['promoter_analysis_dir']}/shared_peak_scaled.tab",
        Shared_UPgene=f"{config['paths']['promoter_analysis_dir']}/shared_up_genes_INpromote.txt",
        Shared_DOWNgene=f"{config['paths']['promoter_analysis_dir']}/shared_down_genes_INpromote.txt"
    output:
        pdf=f"{config['paths']['promoter_analysis_dir']}/Expr_ccess_corr_plot.pdf"
    params:
        Log2FC=config["params"]["log2fc_up"],
        padj=config["params"]["padj"],
        col_ranges=",".join(config["params"]["group_analysis"]["column_ranges"]),
        sample_names=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["sample_labels"]]),
        Agroup=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["treatment_group"]]),
        Ngroup=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["control_group"]])
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/plot_correlation.R \
            {input.DEgene} \
            {params.Log2FC} \
            {params.padj} \
            {input.peakanno} \
            {input.Peaktab} \
            {params.col_ranges} \
            {params.sample_names} \
            {params.Agroup} \
            {params.Ngroup} \
            {input.Shared_UPgene} \
            {input.Shared_DOWNgene} \
            {output.pdf}
        """

rule sharedUpGene_INPromote_GSEA:
    input:
        UPgene=f"{config['paths']['promoter_analysis_dir']}/shared_up_genes_INpromote.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        GOenrich=f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_GOenrich.csv",
        GOrds=f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_GOenrich.rds",
        KEGGenrich=f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_KEGGenrich.csv",
        KEGGrds=f"{config['params']['enrichment_analysis']['promoter_up_dir']}/sharedUpGene_InPromote_gsea_KEGGenrich.rds"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/enrichment_analysis/GSEA_enrichment_shared_UPgene.R \
            {input.UPgene} \
            {input.DEgene} \
            {output.GOenrich} \
            {output.KEGGenrich} \
            {output.GOrds} \
            {output.KEGGrds}
        """

rule sharedDownGene_INPromote_GSEA:
    input:
        DOWNgene=f"{config['paths']['promoter_analysis_dir']}/shared_down_genes_INpromote.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        GOenrich=f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_GOenrich.csv",
        GOrds==f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_GOenrich.rds",
        KEGGenrich=f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_KEGGenrich.csv",
        KEGGrds=f"{config['params']['enrichment_analysis']['promoter_down_dir']}/sharedDownGene_InPromote_gsea_KEGGenrich.rds"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/enrichment_analysis/GSEA_enrichment_shared_DOWNgene.R \
            {input.DOWNgene} \
            {input.DEgene} \
            {output.GOenrich} \
            {output.KEGGenrich} \
            {output.GOrds} \
            {output.KEGGrds}
        """

rule format_TF_motifs:
    input:
        Raw=config["paths"]["jaspar_raw"]
    output:
        JASPAR=f"{config['paths']['jaspar_dir']}/Alx1.meme"
    params:
        DIR=config["paths"]["jaspar_dir"]
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/format_TF_motifs.pl \
            --motiffile {input.Raw} \
            --outdir {params.DIR}
        """

rule extr_shared_UpGene_promoterSeq:
    input:
        UPgene=f"{config['paths']['promoter_analysis_dir']}/shared_up_genes_INpromote.txt"
    output:
        fa=f"{config['params']['regulation_analysis']['promoter_up_dir']}/promoter_seqs.fa"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/extr_shared_UpGene_promoterSeq.R \
            {input.UPgene} \
            {output.fa}
        """

rule motif_enrich_for_sharedUpGene_INpromoter:
    input:
        fa=f"{config['params']['regulation_analysis']['promoter_up_dir']}/promoter_seqs.fa",
        motif_dir=config["paths"]["jaspar_dir"]
    output:
        complete=f"{config['params']['regulation_analysis']['promoter_up_dir']}/motif_enrich_result/analysis_complete.txt",
        flag=touch("flags/motif_enrich_for_sharedUpGene_INpromoter.done")
    params:
        output_dir=f"{config['params']['regulation_analysis']['promoter_up_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/tf_motif_enrich_analysis.pl \
            --fastafile {input.fa} \
            --outputdir {params.output_dir} \
            --motifdir {input.motif_dir}
        
        echo "Analysis completed on $(date)" > {output.complete}
        """

rule TF_Gene_network_for_sharedUpGene_INpromoter:
    input:
        motif_analysis=f"{config['params']['regulation_analysis']['promoter_up_dir']}/motif_enrich_result/analysis_complete.txt"
    output:
        network=f"{config['params']['regulation_analysis']['promoter_up_dir']}/SharedUPgene_in_promoter_TF_gene_network.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['promoter_up_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/build_TF_Gene_motif_network.pl \
            --motifdir {params.motif_dir} \
            --outfile {output.network}
        """

rule extract_enriched_motifs_for_sharedUpGene_INpromoter:
    input:
        flag="flags/motif_enrich_for_sharedUpGene_INpromoter.done",
        tfmapfile=config["paths"]["tf_mapping_file"]
    output:
        enriched_tfs=f"{config['params']['regulation_analysis']['promoter_up_dir']}/enriched_tfs_list.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['promoter_up_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/extr_enriched_motifs.pl \
            --tfmapfile {input.tfmapfile} \
            --motifdir {params.motif_dir} \
            --outfile {output.enriched_tfs}
        """

rule extr_shared_DownGene_promoterSeq:
    input:
        DOWNgene=f"{config['paths']['promoter_analysis_dir']}/shared_down_genes_INpromote.txt"
    output:
        fa=f"{config['params']['regulation_analysis']['promoter_down_dir']}/promoter_seqs.fa"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/extr_shared_DownGene_promoterSeq.R \
            {input.DOWNgene} \
            {output.fa}
        """

rule motif_enrich_for_sharedDownGene_INpromoter:
    input:
        fa=f"{config['params']['regulation_analysis']['promoter_down_dir']}/promoter_seqs.fa",
        motif_dir=config["paths"]["jaspar_dir"]
    output:
        complete=f"{config['params']['regulation_analysis']['promoter_down_dir']}/motif_enrich_result/analysis_complete.txt",
        flag=touch("flags/motif_enrich_for_sharedDownGene_INpromoter.done")
    params:
        output_dir=f"{config['params']['regulation_analysis']['promoter_down_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/tf_motif_enrich_analysis.pl \
            --fastafile {input.fa} \
            --outputdir {params.output_dir} \
            --motifdir {input.motif_dir}
        
        echo "Analysis completed on $(date)" > {output.complete}
        """

rule TF_Gene_network_for_sharedDownGene_INpromoter:
    input:
        motif_analysis=f"{config['params']['regulation_analysis']['promoter_down_dir']}/motif_enrich_result/analysis_complete.txt"
    output:
        network=f"{config['params']['regulation_analysis']['promoter_down_dir']}/SharedDOWNgene_in_promoter_TF_gene_network.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['promoter_down_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/build_TF_Gene_motif_network.pl \
            --motifdir {params.motif_dir} \
            --outfile {output.network}
        """

rule extract_enriched_motifs_for_sharedDownGene_INpromoter:
    input:
        flag="flags/motif_enrich_for_sharedDownGene_INpromoter.done",
        tfmapfile=config["paths"]["tf_mapping_file"]
    output:
        enriched_tfs=f"{config['params']['regulation_analysis']['promoter_down_dir']}/enriched_tfs_list.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['promoter_down_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/extr_enriched_motifs.pl \
            --tfmapfile {input.tfmapfile} \
            --motifdir {params.motif_dir} \
            --outfile {output.enriched_tfs}
        """

rule TF_heatmap_for_sharedUpGene_INpromoter: 
    input:
        enriched_tfs=f"{config['params']['regulation_analysis']['promoter_up_dir']}/enriched_tfs_list.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        heatmap=f"{config['params']['regulation_analysis']['promoter_up_dir']}/shared_UpGene_InPromote_tfs_heatmap_RNAseqExpr.pdf"
    params:
        up_log2fc=config["params"]["log2fc_up"],
        up_padj=config["params"]["padj"],
        down_log2fc=config["params"]["log2fc_down"],
        down_padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/TFheatmap.R \
            {input.enriched_tfs} \
            {input.DEgene} \
            {params.up_log2fc} \
            {params.up_padj} \
            {params.down_log2fc} \
            {params.down_padj} \
            {output.heatmap}
        """

rule TF_heatmap_for_sharedDownGene_INpromoter:
    input:
        enriched_tfs=f"{config['params']['regulation_analysis']['promoter_down_dir']}/enriched_tfs_list.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        heatmap=f"{config['params']['regulation_analysis']['promoter_down_dir']}/shared_DownGene_InPromote_tfs_heatmap_RNAseqExpr.pdf"
    params:
        up_log2fc=config["params"]["log2fc_up"],
        up_padj=config["params"]["padj"],
        down_log2fc=config["params"]["log2fc_down"],
        down_padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/TFheatmap.R \
            {input.enriched_tfs} \
            {input.DEgene} \
            {params.up_log2fc} \
            {params.up_padj} \
            {params.down_log2fc} \
            {params.down_padj} \
            {output.heatmap}
        """

rule get_DEGene_anno:
    input:
        gene_anno=config["paths"]["geneanno_bed"],
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        UPgene_anno=f"{config['paths']['enhancer_analysis_dir']}/UpGene_anno.bed",
        DOWNgene_anno=f"{config['paths']['enhancer_analysis_dir']}/DownGene_anno.bed"
    params:
        up_log2fc=config["params"]["log2fc_up"],
        up_padj=config["params"]["padj"],
        down_log2fc=config["params"]["log2fc_down"],
        down_padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/get_DEGene_anno.R \
            {input.gene_anno} \
            {input.DEgene} \
            {params.up_log2fc} \
            {params.up_padj} \
            {output.UPgene_anno} \
            {params.down_log2fc} \
            {params.down_padj} \
            {output.DOWNgene_anno}
        """

rule calculate_eRNA_UpGene_distance:
    input:
        enhfile=config["paths"]["enhancer_bed"],
        genefile=f"{config['paths']['enhancer_analysis_dir']}/UpGene_anno.bed"
    output:
        outfile=f"{config['paths']['enhancer_analysis_dir']}/enh_UpGene_within_1mb.txt"
    shell:
        """
        perl {config[paths][scripts_dir]}/identificate_shared_genes/calculate_eRNA_gene_distance.pl \
            --enhfile {input.enhfile} \
            --genefile {input.genefile} \
            --outfile {output.outfile}
        """

rule calculate_eRNA_DownGene_distance:
    input:
        enhfile=config["paths"]["enhancer_bed"],
        genefile=f"{config['paths']['enhancer_analysis_dir']}/DownGene_anno.bed"
    output:
        outfile=f"{config['paths']['enhancer_analysis_dir']}/enh_DownGene_within_1mb.txt"
    shell:
        """
        perl {config[paths][scripts_dir]}/identificate_shared_genes/calculate_eRNA_gene_distance.pl \
            --enhfile {input.enhfile} \
            --genefile {input.genefile} \
            --outfile {output.outfile}
        """

rule get_enhancer_anno:
    input:
        enh_UpGene=f"{config['paths']['enhancer_analysis_dir']}/enh_UpGene_within_1mb.txt",
        enh_DownGene=f"{config['paths']['enhancer_analysis_dir']}/enh_DownGene_within_1mb.txt",
        enhancer_ref=config["paths"]["enhancer_bed"]
    output:
        enh_anno_UpGene=f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_UpGene.bed",
        enh_anno_DownGene=f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_DownGene.bed"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/get_enhancer_anno.R \
            {input.enh_UpGene} \
            {input.enhancer_ref} \
            {output.enh_anno_UpGene} \
            {input.enh_DownGene} \
            {output.enh_anno_DownGene}
        """

rule calculate_enhancer_UpGene_signal:
    input:
        bw=expand(f"{config['paths']['atacseq_bowtie2_dir']}/{{sample}}.bw", sample=SAMPLES),
        bedfile=f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_UpGene.bed"
    output:
        gz=f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_UPgene.gz",
        tab=f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_UPgene.tab",
        sorted_regions=f"{config['paths']['enhancer_analysis_dir']}/enh_sortedRegions_for_UPgene.bed"
    params:
        before=config["params"]["computematrix"]["before"],
        after=config["params"]["computematrix"]["after"],
        bin_size=config["params"]["computematrix"]["bin_size"],
        processors=config["params"]["computematrix"]["processors"]
    run:
        bw_files = " ".join(input.bw)
        shell_cmd = f"computeMatrix reference-point -R {input.bedfile} -S {bw_files} -b {params.before} -a {params.after} --binSize {params.bin_size} -p {params.processors} -o {output.gz} --outFileNameMatrix {output.tab} --outFileSortedRegions {output.sorted_regions}"
        shell(shell_cmd)

rule calculate_enhancer_DownGene_signal:
    input:
        bw=expand(f"{config['paths']['atacseq_bowtie2_dir']}/{{sample}}.bw", sample=SAMPLES),
        bedfile=f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_DownGene.bed"
    output:
        gz=f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_DOWNgene.gz",
        tab=f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_DOWNgene.tab",
        sorted_regions=f"{config['paths']['enhancer_analysis_dir']}/enh_sortedRegions_for_DOWNgene.bed"
    params:
        before=config["params"]["computematrix"]["before"],
        after=config["params"]["computematrix"]["after"],
        bin_size=config["params"]["computematrix"]["bin_size"],
        processors=config["params"]["computematrix"]["processors"]
    run:
        bw_files = " ".join(input.bw)
        shell_cmd = f"computeMatrix reference-point -R {input.bedfile} -S {bw_files} -b {params.before} -a {params.after} --binSize {params.bin_size} -p {params.processors} -o {output.gz} --outFileNameMatrix {output.tab} --outFileSortedRegions {output.sorted_regions}"
        shell(shell_cmd)

rule get_sharedUp_enhancer_and_gene:
    input:
        bed=f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_UpGene.bed",
        sortedRegions=f"{config['paths']['enhancer_analysis_dir']}/enh_sortedRegions_for_UPgene.bed",
        tab=f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_UPgene.tab",
        txt=f"{config['paths']['enhancer_analysis_dir']}/enh_UpGene_within_1mb.txt"
    output:
        enhancer=f"{config['paths']['enhancer_analysis_dir']}/SharedUPenhancer.txt",
        gene=f"{config['paths']['enhancer_analysis_dir']}/SharedUPgene_INenhancer.txt"
    params:
        col_ranges=",".join(config["params"]["group_analysis"]["column_ranges"]),
        sample_names=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["sample_labels"]]),
        Agroup=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["treatment_group"]]),
        Ngroup=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["control_group"]])
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/get_sharedUp_enhancer_and_gene.R \
            {input.bed} \
            {input.sortedRegions} \
            {input.tab} \
            {params.col_ranges} \
            {params.sample_names} \
            {params.Agroup} \
            {params.Ngroup} \
            {input.txt} \
            {output.enhancer} \
            {output.gene}
        """

rule get_sharedDown_enhancer_and_gene:
    input:
        bed=f"{config['paths']['enhancer_analysis_dir']}/annotation_of_enh_for_potential_regulate_DownGene.bed",
        sortedRegions=f"{config['paths']['enhancer_analysis_dir']}/enh_sortedRegions_for_DOWNgene.bed",
        tab=f"{config['paths']['enhancer_analysis_dir']}/enh_expr_signal_for_DOWNgene.tab",
        txt=f"{config['paths']['enhancer_analysis_dir']}/enh_DownGene_within_1mb.txt"
    output:
        enhancer=f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNenhancer.txt",
        gene=f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNgene_INenhancer.txt"
    params:
        col_ranges=",".join(config["params"]["group_analysis"]["column_ranges"]),
        sample_names=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["sample_labels"]]),
        Agroup=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["treatment_group"]]),
        Ngroup=",".join([f'"{name}"' for name in config["params"]["group_analysis"]["control_group"]])
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/get_sharedDown_enhancer_and_gene.R \
            {input.bed} \
            {input.sortedRegions} \
            {input.tab} \
            {params.col_ranges} \
            {params.sample_names} \
            {params.Agroup} \
            {params.Ngroup} \
            {input.txt} \
            {output.enhancer} \
            {output.gene}
        """

rule shared_genes_INenhancer_volcano_plot:
    input:
        UPgene=f"{config['paths']['enhancer_analysis_dir']}/SharedUPgene_INenhancer.txt",
        DOWNgene=f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNgene_INenhancer.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        volcano_plot=f"{config['paths']['enhancer_analysis_dir']}/shared_genes_INenhancer_volcano_plot.pdf"
    params:
        Log2FC=config["params"]["log2fc_up"],
        padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/identificate_shared_genes/shared_gene_volcano_plot.R \
            {input.DEgene} \
            {input.UPgene} \
            {input.DOWNgene} \
            {params.Log2FC} \
            {params.padj} \
            {output.volcano_plot}
        """

rule sharedUpGene_INenhancer_GSEA:
    input:
        UPgene=f"{config['paths']['enhancer_analysis_dir']}/SharedUPgene_INenhancer.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        GOenrich=f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_GOenrich.csv",
        GOrds=f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_GOenrich.rds",
        KEGGenrich=f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_KEGGenrich.csv",
        KEGGrds=f"{config['params']['enrichment_analysis']['enhancer_up_dir']}/sharedUpGene_InEnhancer_gsea_KEGGenrich.rds"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/enrichment_analysis/GSEA_enrichment_shared_UPgene.R \
            {input.UPgene} \
            {input.DEgene} \
            {output.GOenrich} \
            {output.KEGGenrich} \
            {output.GOrds} \
            {output.KEGGrds}
        """

rule sharedDownGene_INenhancer_GSEA:
    input:
        DOWNgene=f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNgene_INenhancer.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        GOenrich=f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_GOenrich.csv",
        GOrds=f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_GOenrich.rds",
        KEGGenrich=f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_KEGGenrich.csv",
        KEGGrds=f"{config['params']['enrichment_analysis']['enhancer_down_dir']}/sharedDownGene_InEnhancer_gsea_KEGGenrich.rds"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/enrichment_analysis/GSEA_enrichment_shared_DOWNgene.R \
            {input.DOWNgene} \
            {input.DEgene} \
            {output.GOenrich} \
            {output.KEGGenrich} \
            {output.GOrds} \
            {output.KEGGrds}
        """

rule extract_enhancer_sequences:
    input:
        shared_up_enh=f"{config['paths']['enhancer_analysis_dir']}/SharedUPenhancer.txt",
        shared_down_enh=f"{config['paths']['enhancer_analysis_dir']}/SharedDOWNenhancer.txt",
        enhancer_ref=config["paths"]["enhancer_bed"]
    output:
        up_sequences=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/SharedUP_enhancer_sequences.fa",
        down_sequences=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/SharedDOWN_enhancer_sequences.fa"
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/extrEnhancerSeq.R \
            {input.shared_up_enh} \
            {input.enhancer_ref} \
            {output.up_sequences} \
            {input.shared_down_enh} \
            {output.down_sequences}
        """

rule motif_enrich_for_sharedUpEnhancer:
    input:
        fa=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/SharedUP_enhancer_sequences.fa",
        motif_dir=config["paths"]["jaspar_dir"]
    output:
        complete=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/motif_enrich_result/analysis_complete.txt",
        flag=touch("flags/motif_enrich_for_sharedUpEnhancer.done")
    params:
        output_dir=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/tf_motif_enrich_analysis.pl \
            --fastafile {input.fa} \
            --outputdir {params.output_dir} \
            --motifdir {input.motif_dir}
        
        echo "Analysis completed on $(date)" > {output.complete}
        """

rule TF_Gene_network_for_sharedUpEnhancer:
    input:
        motif_analysis=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/motif_enrich_result/analysis_complete.txt"
    output:
        network=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/SharedUPenhancer_TF_gene_network.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/build_TF_Gene_motif_network.pl \
            --motifdir {params.motif_dir} \
            --outfile {output.network}
        """

rule extract_enriched_motifs_for_sharedUpEnhancer:
    input:
        flag="flags/motif_enrich_for_sharedUpEnhancer.done",
        tfmapfile=config["paths"]["tf_mapping_file"]
    output:
        enriched_tfs=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/enriched_tfs_list.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/extr_enriched_motifs.pl \
            --tfmapfile {input.tfmapfile} \
            --motifdir {params.motif_dir} \
            --outfile {output.enriched_tfs}
        """

rule motif_enrich_for_sharedDownEnhancer:
    input:
        fa=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/SharedDOWN_enhancer_sequences.fa",
        motif_dir=config["paths"]["jaspar_dir"]
    output:
        complete=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/motif_enrich_result/analysis_complete.txt",
        flag=touch("flags/motif_enrich_for_sharedDownEnhancer.done")
    params:
        output_dir=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/tf_motif_enrich_analysis.pl \
            --fastafile {input.fa} \
            --outputdir {params.output_dir} \
            --motifdir {input.motif_dir}
        
        echo "Analysis completed on $(date)" > {output.complete}
        """

rule TF_Gene_network_for_sharedDownEnhancer:
    input:
        motif_analysis=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/motif_enrich_result/analysis_complete.txt"
    output:
        network=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/SharedDOWNenhancer_TF_gene_network.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/build_TF_Gene_motif_network.pl \
            --motifdir {params.motif_dir} \
            --outfile {output.network}
        """

rule extract_enriched_motifs_for_sharedDownEnhancer:
    input:
        flag="flags/motif_enrich_for_sharedDownEnhancer.done",
        tfmapfile=config["paths"]["tf_mapping_file"]
    output:
        enriched_tfs=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/enriched_tfs_list.txt"
    params:
        motif_dir=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/motif_enrich_result/"
    shell:
        """
        perl {config[paths][scripts_dir]}/regulation_analysis/extr_enriched_motifs.pl \
            --tfmapfile {input.tfmapfile} \
            --motifdir {params.motif_dir} \
            --outfile {output.enriched_tfs}
        """

rule TF_heatmap_for_sharedUpEnhancer:  
    input:
        enriched_tfs=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/enriched_tfs_list.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        heatmap=f"{config['params']['regulation_analysis']['enhancer_up_dir']}/shared_UpEnhancer_tfs_heatmap_RNAseqExpr.pdf"
    params:
        up_log2fc=config["params"]["log2fc_up"],
        up_padj=config["params"]["padj"],
        down_log2fc=config["params"]["log2fc_down"],
        down_padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/TFheatmap.R \
            {input.enriched_tfs} \
            {input.DEgene} \
            {params.up_log2fc} \
            {params.up_padj} \
            {params.down_log2fc} \
            {params.down_padj} \
            {output.heatmap}
        """

rule TF_heatmap_for_sharedDownEnhancer:    
    input:
        enriched_tfs=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/enriched_tfs_list.txt",
        DEgene=f"{config['input_files']['deg_csv']}"
    output:
        heatmap=f"{config['params']['regulation_analysis']['enhancer_down_dir']}/shared_DownEnhancer_tfs_heatmap_RNAseqExpr.pdf"
    params:
        up_log2fc=config["params"]["log2fc_up"],
        up_padj=config["params"]["padj"],
        down_log2fc=config["params"]["log2fc_down"],
        down_padj=config["params"]["padj"]
    shell:
        """
        Rscript --vanilla {config[paths][scripts_dir]}/regulation_analysis/TFheatmap.R \
            {input.enriched_tfs} \
            {input.DEgene} \
            {params.up_log2fc} \
            {params.up_padj} \
            {params.down_log2fc} \
            {params.down_padj} \
            {output.heatmap}
        """
