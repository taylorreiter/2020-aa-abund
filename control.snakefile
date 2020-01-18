ENV = 'env.yml'
#nbhd = ['hu-genome19', 'hu-genome20', 'hu-genome21', 'hu-genome22', 
#        'hu-genome23', 'hu-genome24', 'hu-genome25', 'hu-genome26',
#        'hu-genome27', 'hu-genome28', 'hu-genome29', 'hu-genome30',
#        'hu-genome31', 'hu-genome32', 'hu-genome33', 'hu-genome34',
#        'hu-genome35', 'hu-genome36', 'hu-genome37', 'hu-genome38',
#        'hu-genome39', 'hu-genome40', 'hu-genome41']

rule all: 
    input:
        "fake_out_plass.txt",
        #"fake_out_reads.txt"

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
    cd-hit -i {input} -o {output} -c .95
    '''

rule download_reads:
    output: "inputs/reads/hu-s1-reads-hardtrim-jan08.2019.tar.gz"
    shell:'''
    curl -L -o {output} https://osf.io/x5ezk/download
    '''

#rule decompress_reads:
#    output: "inputs/reads/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz"
#    input: "inputs/reads/hu-s1-reads-hardtrim-jan08.2019.tar.gz"
#    params:
#        folder = "inputs/reads"
#    shell:'''
#    tar xvf {input} -C {params.folder} 
#    '''

#def aggregate_decompress_reads(wildcards):
#    # checkpoint_output produces the output dir from the checkpoint rule.
#    checkpoint_output = checkpoints.decompress_reads.get(**wildcards).output[0]    
#    file_names = expand("inputs/reads/hu-s1_k31_r1_search_oh0/{nbhd_reads}.fa.cdbg_ids.reads.hardtrim.fa.gz", 
#                        nbhd_reads = glob_wildcards(os.path.join(checkpoint_output, "{nbhd_reads}.fa.cdbg_ids.reads.hardtrim.fa.gz")).nbhd_reads)
#    return file_names

rule paladin_index_plass:
    output: "outputs/cd-hit95/{nbhd}.faa.bwt"
    input: "outputs/cd-hit95/{nbhd}.faa"
    benchmark: "benchmarks/{nbhd}.paladin_index.txt"
    conda: ENV
    shell:'''
    paladin index -r3 {input}
    '''

rule paladin_align:
    output: "outputs/paladin/{nbhd}.sam"
    input:
        indx="outputs/cd-hit95/{nbhd}.faa.bwt"
        reads="inputs/reads/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz"
        #reads = aggregate_decompress_reads
    params:
        indx= "outputs/cd-hit95/{nbhd}.faa"
    benchmark: "benchmarks/{nbhd}.paladin_align.txt"
    conda: ENV
    shell:'''
    paladin align -f 125 -t 2 {params.indx} {input.reads} > {output}
    '''

rule samtools_flagstat_paladin:
    output: "outputs/paladin/{nbhd}.sam.flagstat"
    input: "outputs/paladin/{nbhd}.sam"
    conda: ENV
    shell:'''
    samtools flagstat {input} > {output}
    '''

rule samtools_view_paladin:
    output: "outputs/paladin/{nbhd}.bam"
    input: "outputs/paladin/{nbhd}.sam"
    conda: ENV
    shell:'''
    samtools view -b {input} > {output}
    '''

rule salmon_paladin:
    output: "outputs/salmon/{nbhd}_quant/quant.sf"
    input:"outputs/cd-hit95/{nbhd}.faa"
    conda: ENV
    shell:'''
    salmon quant -t {input.cdhit} -l A -a {input.bam} -o {params.out}
    '''

def aggregate_decompress_plass(wildcards):
    # checkpoint_output produces the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.decompress_plass.get(**wildcards).output[0]    
    file_names = expand("outputs/salmon/{nbhd}_quant/quant.sf", 
                        nbhd = glob_wildcards(os.path.join(checkpoint_output, "{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup")).nbhd)
    return file_names

rule fake_out:
    output: "fake_out_plass.txt"
    input: aggregate_decompress_plass
    shell:'''
    touch fake_out_plass.txt
    '''

#rule fake_out_reads:
#    output: "fake_out_reads.txt"
#    input: aggregate_decompress_plass
#    shell:'''
#    touch fake_reads_plass.txt
#    '''

