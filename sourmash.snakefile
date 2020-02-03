rule download_plass:
    output: "inputs/plass/hu-s1-plass-hardtrim-clean-jan08.2019.tar.gz"
    shell:'''
    curl -L -o {output} https://osf.io/9hg85/download
    '''

# the "nbhd" expansion var is born from the output of this checkpoint.
# It's not included in this rule because it doesn't exist yet.
# It exists after this rule is run, because function aggregate_decompress_plass_for_cdhit
# looks into that directory, lobs off the suffix of the file specified in the
# function, and uses what's left to create nbhd.
checkpoint decompress_plass:
    output: directory("inputs/plass/hu-s1_k31_r1_search_oh0")
    input: "inputs/plass/hu-s1-plass-hardtrim-clean-jan08.2019.tar.gz"
    params: folder = "inputs/plass"
    shell:'''
    mkdir -p {params.folder}
    tar xvf {input} -C {params.folder}
    '''

rule cdhit_plass:
    output: "outputs/cd-hit95/{nbhd}.cdhit95.faa"
    input: "inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup"
    benchmark: "benchmarks/{nbhd}.cdhit95.benchmark.txt"
    conda: ENV
    shell:'''
    cd-hit -i {input} -o {output} -c .95 -d 0
    '''

checkpoint parse_cdhit_clusters:
    """
    cd-hit produces clusters of proteins, but only the representative sequence
    is output in fasta format. Using the amino acid sequence names, construct
    one amino acid output for each each cluster, where each file contains the
    amino acid sequences of all proteins in the cluster. This rule makes one
    file per cluster, and records the names of the sequences in the cluster in
    the file. 
    """
    input: "outputs/cd-hit95/{nbhd}.cdhit95.faa"
    output: directory("outputs/cd-hit95_clusters")
    conda: "tidyverse.yml"
    script: "scripts/parse_cdhit.R"
 
rule extract_cdhit_clusters:
    """
    Extract amino acid sequences into clusters
    """
    input: "outputs/cd-hit95-clusters"
    output:
    conda:
    script:


rule download_reads:
    output: "inputs/reads/hu-s1-reads-hardtrim-jan08.2019.tar.gz"
    shell:'''
    curl -L -o {output} https://osf.io/x5ezk/download
    '''

def aggregate_decompress_plass(wildcards):
    # checkpoint_output produces the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.decompress_plass.get(**wildcards).output[0]
    file_names = expand("outputs/paladin/{nbhd}.sam",
                        nbhd = glob_wildcards(os.path.join(checkpoint_output, "{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup")).nbhd)
    return file_names

rule fake_out:
    output: "fake_out_plass.txt"
    input: aggregate_decompress_plass
    shell:'''
    touch fake_out_plass.txt
    '''
