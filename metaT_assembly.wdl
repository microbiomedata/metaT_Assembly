version 1.0

import "https://code.jgi.doe.gov/BFoster/jgi_meta/-/raw/nmdc_runtimes/jgi_meta_wdl_sets/metatranscriptome_assembly_and_alignment/mapping.wdl?ref=heads" as mapping # replace with correct URLs once branch nmdc_runtimes created on https://code.jgi.doe.gov/BFoster/jgi_meta/
import "https://code.jgi.doe.gov/BFoster/jgi_meta/-/raw/nmdc_runtimes/jgi_meta_wdl_sets/metatranscriptome_assembly_and_alignment/metatranscriptome_asm_rnaspades.wdl?ref=heads" as http_rnaspades

workflow metatranscriptome_assy {
    input{
        Array[String] input_files
        String proj_id
        String prefix=sub(proj_id, ":", "_")
        String bbtools_container = "bryce911/bbtools:39.10"
        String bbtools_map_container = "bryce911/bbtools:38.86"
        String spades_container = "staphb/spades:4.0.0"
        String workflowmeta_container="microbiomedata/workflowmeta:1.1.1"
        Int asm_thr = 16 # half of defaults
        Int asm_mem = 200
        Int map_thr = 18
        Int map_mem = 120
        Boolean dotar_bams = true
        Boolean dofinalize_bams = true
    }

    call http_rnaspades.metatranscriptome_asm as mtasm{
        input: 
            input_files = input_files,
            bbtools_container = bbtools_container,
            spades_container_prod = spades_container,
            asm_thr = asm_thr,
            asm_mem = asm_mem
    }

    call rename_contig {
        input:
        contigs = mtasm.final_contigs,
        scaffolds = mtasm.final_scaffolds,
        agp = mtasm.final_agp,
        legend = mtasm.final_legend,
        proj_id = proj_id,
        prefix = prefix,
        container = bbtools_container
    }

    call mapping.mapping as map {
        input:
            input_files = input_files,
            input_reference = rename_contig.outcontigs,
            dotar_bams = dotar_bams,
            dofinalize_bams = dofinalize_bams,
            bbtools_container = bbtools_map_container,
            map_thr = map_thr,
            map_mem = map_mem
    }

    call finish_asm {
        input:
        prefix = prefix,
        readlen = mtasm.final_readlen,
        contigs = rename_contig.outcontigs,
        scaffolds = rename_contig.outscaffolds,
        asmstats = rename_contig.asmstats,
        log = mtasm.final_log,
        sam = select_first([map.outsam]),
        bam = select_first([map.outbam]),
        bamidx = select_first([map.outbamidx]),
        cov = select_first([map.outcov]),
        tar_bam = select_first([map.outtarbam]),
        container = workflowmeta_container
    }

    call make_info_file {
        input:
        bbtools_info = rename_contig.outlog,
        spades_info = mtasm.final_log,
        prefix = prefix,
        bbtools_container = bbtools_container,
        spades_container = spades_container
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
        cpu:  2
        maxRetries: 1
        docker: container
        runtime_minutes: 30
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
        memory: "2G"
        cpu:  1
        maxRetries: 1
        docker: container
        runtime_minutes: 30
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
        cpu:  1
        maxRetries: 1
        docker: bbtools_container
        runtime_minutes: 30
    }
}
