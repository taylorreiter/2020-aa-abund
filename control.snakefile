ENV = 'env.yml'

rule download_plass:
    output: "inputs/plass/hu-s1-plass-hardtrim-clean-jan08.2019.tar.gz"
    shell:'''
    curl -L -o {output} https://osf.io/9hg85/download
    '''

rule decompress_plass:
    output: dynamic("inputs/plass/hu-s1_k31_r1_search_oh0/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup")
    input: "inputs/plass/hu-s1-plass-hardtrim-clean-jan08.2019.tar.gz"
    params: folder = "inputs/plass"
    shell:'''
    tar xvf {input} -C {params.folder}
    '''
    
rule cdhit_plass:
    output: "outputs/cd-hit/{nbhd}.fa.cdbg_ids.reads.hardtrim.fa.gz.plass.cdhit.fa.clean.cut.dup-cdhit95.faa"
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
