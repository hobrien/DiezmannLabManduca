#snakemake --cluster "qsub -pe smp {params.num_cores}" -j 20

import os
import re

# read config info into this namespace
include: "config.py"

rule all:
    input:
        "R/BamQC.pdf",
        expand("Mappings/{sample}.counts.txt", sample=samples)       
        
rule download_ref:
    output: 
        [os.path.join('Reference', os.path.basename(i).replace('.gz', '')) for i in lookup.values()]
    params:
        num_cores = 1
    run:
        rev_lookup = dict((os.path.basename(i), i) for i in lookup.values())
        for o in output:
            if os.path.basename(o) == 'ms_ogs.gff':
                ro = rev_lookup[os.path.basename(o)]
                shell("wget -P Reference {ro}")
            elif os.path.basename(o) == 'ms_scaffolds.fa':
                ro = rev_lookup[os.path.basename('.'.join([o, 'gz']))]
                shell("wget -P Reference {ro}")
                ro = os.path.basename(ro)
                shell("gunzip Reference/{ro}")
                shell("perl -pi -e 's/gi.*gb\|(\w+\.\d)\|/$1/' Reference/ms_scaffolds.fa")
            else:    
                ro = rev_lookup[os.path.basename('.'.join([o, 'gz']))]
                shell("wget -P Reference {ro}")
                ro = os.path.basename(ro)
                shell("gunzip Reference/{ro}")

rule index_ref:
    input:
        [os.path.join('Reference', os.path.basename(i).replace('.gz', '')) for i in (lookup['candida_genome'], lookup['manduca_genome'])]
    output:
        index = [os.path.join("Reference", '.'.join(["combined", str(i+1), "ht2"])) for i in range(8)],
        fasta = "Reference/combined.fa"
    params:
        outprefix = "Reference/combined",
        num_cores = 8
    log:
        "Logs/Hisat/hisat_build.txt"
    shell:
        "cat {input} > {output.fasta}; (hisat2-build -p 8 {output.fasta} {params.outprefix}) 2> {log}"

rule gff_to_gtf:
    input:
        [os.path.join('Reference', os.path.basename(i).replace('.gz', '')) for i in (lookup['candida_gff'], lookup['manduca_gff'])]
    output:
        os.path.join('Reference', os.path.basename(lookup['candida_gff']).replace('.gz', '').replace('gff3', 'gtf')),
        os.path.join('Reference', os.path.basename(lookup['manduca_gff']).replace('.gz', '').replace('gff', 'gtf'))
    params:
        num_cores = 1
    log:
        "Logs/Hisat/gff_to_gtf.txt"
    run:
        shell("(gffread {input[0]} -T -o {output[0]} ) 2> {log}")
        shell("(gffread {input[1]} -T -o {output[1]} ) 2> {log}")

rule combine_gtf:        
    input:
        os.path.join('Reference', os.path.basename(lookup['candida_gff']).replace('.gz', '').replace('gff3', 'gtf')),
        os.path.join('Reference', os.path.basename(lookup['manduca_gff']).replace('.gz', '').replace('gff', 'gtf'))
    output:
        "Reference/combined.gtf"
    params:
        num_cores = 1
    shell:
        "cat {input} | grep gene_id > {output}"

rule extract_splice_sites:
    input:
        rules.combine_gtf.output
    output:
        "Reference/splicesites.txt"
    params:
        num_cores = 1
    log:
        "Logs/Hisat/hisat_extract_splice.txt"
    shell:
        "(cat {input} | hisat2_extract_splice_sites.py '-' > {output} ) 2> {log}"

        
rule hisat:
    input:
        forward = "raw_data/{sample}_1.fq.gz",
        reverse = "raw_data/{sample}_2.fq.gz",
        splice_sites = rules.extract_splice_sites.output,
        index = rules.index_ref.output.index
    output:
        "Mappings/{sample}.bam"
    params:
        num_cores = 8,
        index_pref = "Reference/combined",
        index = [os.path.join("Reference", '.'.join(["combined", str(i+1), "ht2"])) for i in range(8)]
    log:
        "Logs/Hisat/{sample}_hisat_map.txt"
    shell:
        "(hisat2 --fr --threads 8 -x {params.index_pref} --known-splicesite-infile "
        "{input.splice_sites} -1 {input.forward} -2 {input.reverse} | samtools view "
        "-S -bo {output} -) 2> {log}"

rule sort_bam:
    input:
        rules.hisat.output
    output:
        "Mappings/{sample}.sort.bam"
    params:
         num_cores = 1
    shell:
        "samtools sort -o {output} {input}"

rule samtools_index:
    input:
        rules.sort_bam.output
    output:
        "Mappings/{sample}.sort.bam.bai"
    params:
         num_cores = 1
    shell:
        "samtools index {input}"

rule make_bed:
    input:
        rules.combine_gtf.output
    output:
        "Reference/combined.bed"
    log:
        "Logs/BamQC/make_bed.txt"
    params:
         num_cores = 1
    shell:
        "(gtf2bed.pl {input} > {output}) 2> {log}"

rule transcript_seqs:
    input:
        bed = rules.make_bed.output,
        fasta = rules.index_ref.output.fasta
    output:
        "Reference/combined_transcripts.fa"
    params:
         num_cores = 1
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -split -name -s "
        "| fold -w 60 | perl -pe 's/\([+-]\)//' > {output}"

rule exclude_candida:
    input:
        rules.transcript_seqs.output
    output:
        "Reference/Candida_genes.txt"
    params:
         num_cores = 1
    shell:
        "grep gene {input} | perl -pe 's/>//' > {output}"
    
rule bam_stats:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output
    output:
        "Mappings/{sample}.stats.txt"
    params:
         num_cores = 1
    shell:
        "bam_stat.py -i {input.bam} > {output}"

rule bam_dist:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output,
        bed=rules.make_bed.output
    output:
        "Mappings/{sample}.dist.txt"
    params:
         num_cores = 1
    log:
        "Logs/BamQC/{sample}_dist.txt"
    shell:
        "(read_distribution.py -r {input.bed} -i {input.bam} > {output}) 2> {log}"

rule bam_strand:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output,
        bed=rules.make_bed.output
    output:
        "Mappings/{sample}.expt.txt"
    params:
         num_cores = 1
    log:
        "Logs/BamQC/{sample}_expt.txt"
    shell:
        "(infer_experiment.py -r {input.bed} -i {input.bam} > {output}) 2> {log}"

rule bam_inner_dist:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output,
        bed=rules.make_bed.output
    output:
        "Mappings/{sample}.inner_distance_freq.txt"
    params:
         num_cores = 1,
         prefix = "Mappings/{sample}"
    log:
        "Logs/BamQC/{sample}_inner_dist.txt"
    shell:
        "(inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} -u 1000 -s 10) 2>{log}"

rule bam_junction_sat:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output,
        bed=rules.make_bed.output
    output:
        "Mappings/{sample}.junctionSaturation_plot.r"
    params:
         num_cores = 1,
         prefix = "Mappings/{sample}"
    log:
        "Logs/BamQC/{sample}_junction_sat.txt"
    shell:
        "(junction_saturation.py -r {input.bed} -i {input.bam} -o {params.prefix}) 2>{log}"

rule bam_cov:
    input:
        bam=rules.sort_bam.output,
        index=rules.samtools_index.output
    output:
        "Mappings/{sample}.cov.txt"
    params:
         num_cores = 1,
    shell:
        "/c8000xd3/rnaseq-heath/bin/samtools idxstats {input.bam} > {output}"
    
rule summarise_qc:
    input:
        expand("Mappings/{sample}.stats.txt", sample=samples),
        expand("Mappings/{sample}.dist.txt", sample=samples),
        expand("Mappings/{sample}.expt.txt", sample=samples),
        expand("Mappings/{sample}.inner_distance_freq.txt", sample=samples),
        expand("Mappings/{sample}.junctionSaturation_plot.r", sample=samples),
        expand("Mappings/{sample}.cov.txt", sample=samples)
    output:
        "Tables/read_numbers.txt",
        "Tables/read_distribution.txt",
        "Tables/read_strand.txt",
        "Tables/read_distance.txt",
        "Tables/junction_sat.txt"
    params:
        num_cores=1,
        stats_folder="Mappings"
    log:
        "Logs/BamQC/summarise_qc.txt"
    shell:
        "(Rscript R/SummariseBamQC.R {params.stats_folder}) 2> {log}"

rule qc_report:
    input:
        "Tables/read_numbers.txt",
        "Tables/read_distribution.txt",
        "Tables/read_strand.txt",
        "Tables/read_distance.txt",
        "Tables/junction_sat.txt"
    output:
        "R/BamQC.pdf"
    params:
         num_cores = 1
    shell:
        "Rscript -e 'library(rmarkdown); render(\"R/BamQC.Rmd\")'"

rule count_reads:
    input:
        bam=rules.sort_bam.output,
        gtf=rules.combine_gtf.output
    output:
        "Mappings/{sample}.counts.txt"
    params:
         num_cores = 1
    shell:
        "htseq-count -f bam -s no -t exon -i gene_id -m intersection-strict {input.bam} {input.gtf} > {output}"
        
