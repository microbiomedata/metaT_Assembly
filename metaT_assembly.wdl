version 1.0

import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/raw/main1.0/metatranscriptome/metatranscriptome_assy_rnaspades.wdl" as http_rnaspades
import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/raw/main1.0/common/mapping.wdl" as mapping

workflow metatranscriptome_assy {
    input{
        Array[File] input_files
        String proj_id
        String prefix=sub(proj_id, ":", "_")
        String bbtools_container = "bryce911/bbtools:38.86"
        String spades_container_prod = "bryce911/spades:3.15.4"
    }

    call http_rnaspades.readstats_raw {
        input:
        reads_files = input_files, 
        container = bbtools_container
    }

    call http_rnaspades.assy {
        input:
        reads_files = input_files,
        container = spades_container_prod
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
        rename_contig_prefix = prefix,
        container = bbtools_container
    }
    call mapping.mappingtask as single_run {
        input:
        reads = input_files[0],
        reference = rename_contig.outcontigs,
        container = bbtools_container
    }
    call mapping.finalize_bams as finalize_bams{
            input:
            insing = single_run.outbamfile,
            container = bbtools_container
    }

    call mapping.tar_bams as tar_bams {
            input:
            insing = single_run.outbamfile,
            container = bbtools_container
    }

    call finish_asm {
        input:
        prefix = prefix,
        tar_bam = tar_bams.outtarbam,
        contigs = create_agp.outcontigs,
        scaffolds = create_agp.outscaffolds,
        log = assy.log,
        readlen = readstats_raw.outreadlen,
        sam = finalize_bams.outsam,
        bam = finalize_bams.outbam,
        container = bbtools_container
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
        File info_file = make_info_file.assyinfo
    }

        meta {
        author: "Migun Shakya, B-GEN, LANL"
        email: "migun@lanl.gov"
        version: "0.0.1"
    }
}

task rename_contig{
    input{
        File contigs
        File scaffolds
        File agp
        File legend
        String proj_id
        String rename_contig_prefix
        String container
    }
    command <<<
        set -oeu pipefail
        grep "Version" /bbmap/README.md | sed 's/#//' 

        if [ "~{proj_id}" != "scaffold" ]; then
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{contigs} > "~{rename_contig_prefix}_contigs.fna"
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{scaffolds} > "~{rename_contig_prefix}_scaffolds.fna"
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{agp} > "~{rename_contig_prefix}.agp"
            sed -e 's/scaffold/~{proj_id}_scf/g' ~{legend} > "~{rename_contig_prefix}_scaffolds.legend"
        fi

    >>>

    output{
        File outcontigs = "~{rename_contig_prefix}_contigs.fna"
        File outscaffolds = "~{rename_contig_prefix}_scaffolds.fna"
        File outagp = "~{rename_contig_prefix}.agp"
        File outlegend = "~{rename_contig_prefix}_scaffolds.legend"
        File outlog = stdout()
    }
    runtime {
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
        String container
    }

    command <<<
        set -oeu pipefail
        ln -s ~{tar_bam} ~{prefix}_bamfiles.tar
        ln -s ~{contigs} ~{prefix}_contigs.fna
        ln -s ~{scaffolds} ~{prefix}_scaffolds.fna
        ln -s ~{log} ~{prefix}_spades.log
        ln -s ~{readlen} ~{prefix}_readlen.txt
        ln -s ~{sam} ~{prefix}_pairedMapped.sam.gz
        ln -s ~{bam} ~{prefix}_pairedMapped_sorted.bam
    >>>

    output{
        File final_tar_bam = "~{prefix}_bamfiles.tar"
        File final_contigs = "~{prefix}_contigs.fna"
        File final_scaffolds = "~{prefix}_scaffolds.fna"
        File final_log = "~{prefix}_spades.log"
        File final_readlen = "~{prefix}_readlen.txt"
        File final_sam = "~{prefix}_pairedMapped.sam.gz"
        File final_bam = "~{prefix}_pairedMapped_sorted.bam"
    }
    runtime{
        docker: container
    }
}


task make_info_file{
    input{
    File bbtools_info
    String spades_info
    String prefix
    String bbtools_container
    String spades_container
    }

    command <<<
    set -oeu pipefail
    bbtools_version=`grep Version ~{bbtools_info}`

    echo -e "Metatranscriptomic Assembly Workflow - Info File" > ~{prefix}metaT_assy.info
    echo -e "This workflow assembles metatranscriptomic reads using a workflow developed by Brian Foster at JGI." >> ~{prefix}metaT_assy.info
    echo -e "The reads are assembled using SPAdes(1):" >> ~{prefix}metaT_assy.info
    echo -e "`head -6 ~{spades_info} | tail -4`" >> ~{prefix}metaT_assy.info
    echo -e "An AGP file is created using fungalrelease.sh (BBTools(2)${bbtools_version})." >> ~{prefix}metaT_assy.info
    echo -e "Assembled reads are mapped using bbmap.sh (BBTools(2)${bbtools_version})." >> ~{prefix}metaT_assy.info

    echo -e "\nThe following are the Docker images used in this workflow:" >> ~{prefix}metaT_assy.info
    echo -e "   ~{bbtools_container}" >> ~{prefix}metaT_assy.info
    echo -e "   ~{spades_container}" >> ~{prefix}metaT_assy.info

    echo -e "\n(1) Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). Spades: A new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology, 19(5), 455-477. https://doi.org/10.1089/cmb.2012.0021" >> ~{prefix}metaT_assy.info
    echo -e "(2) B. Bushnell: BBTools software package, http://bbtools.jgi.doe.gov/" >> ~{prefix}metaT_assy.info

    >>>

    output{
        File assyinfo = "~{prefix}metaT_assy.info"
    }
    runtime{
        docker: bbtools_container
    }
}
