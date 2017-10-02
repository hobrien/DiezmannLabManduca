import os

# read config info into this namespace
include: "config.py"

rule all:
    input:
        [os.path.join("Reference", '.'.join(["combined", str(i+1), "ht2"])) for i in range(8)],
        "Reference/splicesites.txt"
        
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

rule extract_splice_sites:
    input:
        [os.path.join('Reference', os.path.basename(i).replace('.gz', '')) for i in (lookup['candida_gff'], lookup['manduca_gff'])]
    output:
        "Reference/splicesites.txt"
    params:
        num_cores = 8
    log:
        "Logs/Hisat/hisat_extract_splice.txt"
    shell:
        "(cat {input} | hisat2_extract_splice_sites.py '-' > {output} ) 2> {log}"