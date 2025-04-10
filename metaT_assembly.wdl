version 1.0

import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/-/raw/main1.0/common/mapping.wdl?ref=0e589f4dfbb4285089c4c99b422e2eec79185ba6" as mapping
import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/-/raw/main1.0/metatranscriptome/metatranscriptome_assy_rnaspades.wdl?ref=855279a58daccf298bca0372c034f29cf95792d7" as http_rnaspades

workflow metatranscriptome_assy {
    input{
        Array[File] input_files # fastq.gz
        String proj_id
        String prefix=sub(proj_id, ":", "_")
        String bbtools_container = "microbiomedata/bbtools:38.96"
        String spades_container_prod = "bryce911/spades:3.15.2"
        String workflowmeta_container="microbiomedata/workflowmeta:1.1.1"
        Int assy_thr = 8 # half of defaults
        Int assy_mem = 120 # half of defaults
    }

    call stage {
        input:
        input_files = input_files,
        container = workflowmeta_container
    }

    call http_rnaspades.readstats_raw {
        input:
        reads_files = stage.reads_files, 
        container = bbtools_container
    }

    call http_rnaspades.assy {
        input:
        reads_files = stage.reads_files,
        container = spades_container_prod,
        threads = assy_thr,
        memory = assy_mem
    }
    call http_rnaspades.create_agp {
        input:
        contigs_in = assy.out,
        container = bbtools_container
    }
    call rename_contig {
        input:
        contigs = create_agp.outcontigs,
        scaffolds = create_agp.outscaffolds,
        agp = create_agp.outagp,
        legend = create_agp.outlegend,
        proj_id = proj_id,
        prefix = prefix,
        container = bbtools_container
    }

    call mapping.mappingtask as single_run {
      input:
        reads = stage.reads_files[0],
        reference = rename_contig.outcontigs,
        container = bbtools_container
    }

    call mapping.tar_bams as tar_bams {
            input:
            insing = single_run.outbamfile,
            container = bbtools_container
    }

    call mapping.finalize_bams as finalize_bams{
            input:
            insing = single_run.outbamfile,
            container = bbtools_container
    }


    call finish_asm {
        input:
        prefix = prefix,
        tar_bam = tar_bams.outtarbam,
        contigs = rename_contig.outcontigs,
        scaffolds = rename_contig.outscaffolds,
        log = assy.log,
        readlen = readstats_raw.outreadlen,
        sam = finalize_bams.outsam,
        bam = finalize_bams.outbam,
        bamidx = finalize_bams.outbamidx,
        cov = finalize_bams.outcov,
        asmstats = rename_contig.asmstats,
        container = workflowmeta_container
    }

    call make_info_file {
        input:
        bbtools_info = rename_contig.outlog,
        spades_info = assy.log,
        prefix = prefix,
        bbtools_container = bbtools_container,
        spades_container = spades_container_prod
    }


    output {
        File final_tar_bam = finish_asm.final_tar_bam
        File final_contigs = finish_asm.final_contigs  # annotation.input_file
        File final_scaffolds = finish_asm.final_scaffolds
        File final_log = finish_asm.final_log
        File final_readlen = finish_asm.final_readlen
        File final_sam = finish_asm.final_sam
        File final_bam = finish_asm.final_bam
        File final_bamidx = finish_asm.final_bamidx
        File final_cov = finish_asm.final_cov
        File asmstats = finish_asm.final_asmstats
        File info_file = make_info_file.assyinfo
        
    }

        meta {
        author: "Migun Shakya, B-GEN, LANL"
        email: "migun@lanl.gov"
        version: "0.0.1"
    }
}

task stage {
    input {
        Array[File] input_files
        String container
        String single = if (length(input_files) == 1) then "1" else "0"
        String reads_input = "reads.input.fastq.gz"
    }
    command <<<
        if [ ~{single} == 0 ]
        then
            cat ~{sep=" "  input_files} > ~{reads_input}
        else
            ln -s ~{input_files[0]} ./~{reads_input} || ln ~{input_files[0]} ./~{reads_input}
        fi
    >>>
    output {
        Array[File] reads_files = [reads_input]
    }
}

task rename_contig{
    input{
        File contigs
        File scaffolds
        File agp
        File legend
        String proj_id
        String prefix
        String container
    }
    command <<<
        set -oeu pipefail
        grep "Version" /bbmap/README.md | sed 's/#//' 

        if [ "~{proj_id}" != "scaffold" ]; then
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{contigs} > "~{prefix}_contigs.fna"
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{scaffolds} > "~{prefix}_scaffolds.fna"
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{agp} > "~{prefix}.agp"
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{legend} > "~{prefix}_scaffolds.legend"
        fi

        bbstats.sh format=8 in=~{scaffolds} out=stats.json
    >>>

    output{
        File outcontigs = "~{prefix}_contigs.fna"
        File outscaffolds = "~{prefix}_scaffolds.fna"
        File outagp = "~{prefix}.agp"
        File outlegend = "~{prefix}_scaffolds.legend"
        File asmstats = "stats.json"
        File outlog = stdout()
    }
    runtime {
        memory: "10G"
        cpu:  4
        maxRetries: 1
        docker: container
    }
}

task finish_asm {
    input{
        String prefix
        File tar_bam
        File contigs
        File scaffolds
        File log 
        File readlen 
        File sam 
        File bam 
        File bamidx
        File cov
        File asmstats
        String container
    }

    command <<<
        set -oeu pipefail
        ln ~{tar_bam} ~{prefix}_bamfiles.tar || ln -s ~{tar_bam} ~{prefix}_bamfiles.tar
        ln ~{contigs} ~{prefix}_contigs.fna || ln -s ~{contigs} ~{prefix}_contigs.fna
        ln ~{scaffolds} ~{prefix}_scaffolds.fna || ln -s ~{scaffolds} ~{prefix}_scaffolds.fna
        ln ~{log} ~{prefix}_spades.log || ln -s ~{log} ~{prefix}_spades.log
        ln ~{readlen} ~{prefix}_readlen.txt || ln -s ~{readlen} ~{prefix}_readlen.txt
        ln ~{sam} ~{prefix}_pairedMapped.sam.gz || ln -s ~{sam} ~{prefix}_pairedMapped.sam.gz
        ln ~{bam} ~{prefix}_pairedMapped_sorted.bam || ln -s ~{bam} ~{prefix}_pairedMapped_sorted.bam
        ln ~{bamidx} ~{prefix}_pairedMapped_sorted.bam.bai|| ln -s ~{bamidx} ~{prefix}_pairedMapped_sorted.bam.bai
        ln ~{cov} ~{prefix}_pairedMapped_sorted.bam.cov || ln -s ~{cov} ~{prefix}_pairedMapped_sorted.bam.cov

        sed -i 's/l_gt50k/l_gt50K/g' ~{asmstats}
        cat ~{asmstats} | jq 'del(.filename)' > scaffold_stats.json

    >>>

    output{
        File final_tar_bam = "~{prefix}_bamfiles.tar"
        File final_contigs = "~{prefix}_contigs.fna"
        File final_scaffolds = "~{prefix}_scaffolds.fna"
        File final_log = "~{prefix}_spades.log"
        File final_readlen = "~{prefix}_readlen.txt"
        File final_sam = "~{prefix}_pairedMapped.sam.gz"
        File final_bam = "~{prefix}_pairedMapped_sorted.bam"
        File final_bamidx = "~{prefix}_pairedMapped_sorted.bam.bai"
        File final_cov = "~{prefix}_pairedMapped_sorted.bam.cov"
        File final_asmstats = "scaffold_stats.json"
        
    }
    runtime{
        memory: "10G"
        cpu:  4
        maxRetries: 1
        docker: container
    }
}


task make_info_file{
    input{
    File bbtools_info
    File spades_info
    String prefix
    String bbtools_container
    String spades_container
    }

    command <<<
    set -oeu pipefail
    bbtools_version=`grep Version ~{bbtools_info}`

    echo -e "Metatranscriptomic Assembly Workflow - Info File" > ~{prefix}_metaT_assy.info
    echo -e "This workflow assembles metatranscriptomic reads using a workflow developed by Brian Foster at JGI." >> ~{prefix}_metaT_assy.info
    echo -e "The reads are assembled using SPAdes(1):" >> ~{prefix}_metaT_assy.info
    echo -e "`head -6 ~{spades_info} | tail -4`" >> ~{prefix}_metaT_assy.info
    echo -e "An AGP file is created using fungalrelease.sh (BBTools(2)${bbtools_version})." >> ~{prefix}_metaT_assy.info
    echo -e "Assembled reads are mapped using bbmap.sh (BBTools(2)${bbtools_version})." >> ~{prefix}_metaT_assy.info

    echo -e "\nThe following are the Docker images used in this workflow:" >> ~{prefix}_metaT_assy.info
    echo -e "   ~{bbtools_container}" >> ~{prefix}_metaT_assy.info
    echo -e "   ~{spades_container}" >> ~{prefix}_metaT_assy.info

    echo -e "\n(1) Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). Spades: A new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology, 19(5), 455-477. https://doi.org/10.1089/cmb.2012.0021" >> ~{prefix}_metaT_assy.info
    echo -e "(2) B. Bushnell: BBTools software package, http://bbtools.jgi.doe.gov/" >> ~{prefix}_metaT_assy.info

    >>>

    output{
        File assyinfo = "~{prefix}_metaT_assy.info"
    }
    runtime{
        memory: "2G"
        cpu:  4
        maxRetries: 1
        docker: bbtools_container
    }
}
