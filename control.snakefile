ENV = 'env.yml'
nbhd = ['hu-genome19', 'hu-genome20', 'hu-genome21', 'hu-genome22', 
        'hu-genome23', 'hu-genome24', 'hu-genome25', 'hu-genome26',
        'hu-genome27', 'hu-genome28', 'hu-genome29', 'hu-genome30',
        'hu-genome31', 'hu-genome32', 'hu-genome33', 'hu-genome34',
        'hu-genome35', 'hu-genome36', 'hu-genome37', 'hu-genome38',
        'hu-genome39', 'hu-genome40', 'hu-genome41']

rule all: 
    input:
        "fake_out.txt"

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

def aggregate_decompress_plass(wildcards):
    checkpoint_output = checkpoints.decompress_plass.get(**wildcards).output[0]    
    file_names = expand("outputs/cd-hit95/{nbhd}.cdhit95.faa", 
                        nbhd = glob_wildcards(os.path.join(checkpoint_output, "{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup")).nbhd)
    print(file_names)
    return file_names


rule cdhit_plass:
    output: "outputs/cd-hit95/{nbhd}.cdhit95.faa"
    input: "inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup"
    benchmark: "benchmarks/{nbhd}.cdhit95.benchmark.txt"
    conda: ENV
    shell:'''
    cd-hit -i {input} -o {output} -c .95
    '''

rule fake_out:
    output: "fake_out.txt"
    input: aggregate_decompress_plass
    shell:'''
    touch fake_out.txt
    '''

rule download_reads:
    output: "inputs/reads/hu-s1-reads-hardtrim-jan08.2019.tar.gz"
    shell:'''
    curl -L -o {output} https://osf.io/x5ezk/download
    '''

rule decompress_reads:
    output: "inputs/reads/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz"
    input: "inputs/reads/hu-s1-reads-hardtrim-jan08.2019.tar.gz"
    params:
        folder = "inputs/reads",
        out = "hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz"
    shell:'''
    tar xvf {input} -C {params.folder} {params.out}
    '''

rule paladin_index_plass:
    output: "inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup.bwt"
    input: "inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup"
    conda: ENV
    shell:'''
    paladin index -r3 {input}
    '''

rule paladin_align:
    output: "outputs/{nbhd}/paladin/{nbhd}.hardtrim.sam"
    input:
        indx="inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup.bwt",
        reads="inputs/reads/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz"
    params:
        indx="inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup"
    conda: ENV
    shell:'''
    paladin align -f 125 -t 2 {params.indx} {input.reads} > {output}
    '''
