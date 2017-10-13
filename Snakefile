#snakemake --cluster "qsub -pe smp {params.num_cores}" -j 20

import os
import re

# read config info into this namespace
include: "config.py"

rule all:
    input:
        expand("Mappings/{sample}.stats.txt", sample=samples),
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
            else:    
                ro = rev_lookup[os.path.basename('.'.join([o, 'gz']))]
                shell("wget -P Reference {ro}")
                ro = os.path.basename(ro)
                shell("gunzip Reference/{ro}")

rule index_ref:
    input:
        [os.path.join('Reference', os.path.basename(i).replace('.gz', '')) for i in (lookup['candida_genome'], lookup['manduca_genome'])]
    output:
        [os.path.join("Reference", '.'.join(["combined", str(i+1), "ht2"])) for i in range(8)]
    params:
        infiles = ','.join([os.path.join('Reference', os.path.basename(i).replace('.gz', '')) for i in (lookup['candida_genome'], lookup['manduca_genome'])]),
        outprefix = "Reference/combined",
        num_cores = 8
    log:
        "Logs/Hisat/hisat_build.txt"
    shell:
        "(hisat2-build -p 8 {params.infiles} {params.outprefix}) 2> {log}"

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
        splice_sites = rules.extract_splice_sites.output
    output:
        "Mappings/{sample}.bam"
    params:
        num_cores = 8,
        index_pref = "Reference/combined",
        index = [os.path.join("Reference", '.'.join(["combined", str(i+1), "ht2"])) for i in range(8)]
    log:
        "Logs/Hisat/hisat_map.txt"
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

rule bam_stats:
    input:
        rules.sort_bam.output,
        rules.samtools_index.output
    output:
        "Mappings/{sample}.stats.txt"
    params:
         num_cores = 1
    shell:
        "bam_stat.py -i {input} > {output}"
        
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
