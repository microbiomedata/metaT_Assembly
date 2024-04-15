version 1.0

import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/raw/main1.0/metatranscriptome/metatranscriptome_assy_rnaspades.wdl" as http_rnaspades
import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/raw/main1.0/common/mapping.wdl" as mapping

workflow metatranscriptome_assy {
    input{
        Array[File] input_files
        String rename_contig_prefix="scaffold"
        String bbtools_container = "bryce911/bbtools:38.86"
        String spades_container_prod = "bryce911/spades:3.15.4"
        String prefix = ""
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
        rename_contig_prefix = rename_contig_prefix,
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

    call make_info_file {
        input:
        bbtools_info = rename_contig.outlog,
        spades_info = assy.log,
        prefix = prefix,
        bbtools_container = bbtools_container,
        spades_container = spades_container_prod
    }


    output {
        File final_tar_bam = tar_bams.outtarbam
        File final_contigs = create_agp.outcontigs  # annotation.input_file
        File final_scaffolds = create_agp.outscaffolds
        File final_log = assy.log
        File final_readlen = readstats_raw.outreadlen
        File final_sam = finalize_bams.outsam
        File final_bam = finalize_bams.outbam
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
        String rename_contig_prefix
        String container
        String filename_contigs="~{rename_contig_prefix}_contigs.fna"
        String filename_scaffolds="~{rename_contig_prefix}_scaffolds.fna"
        String filename_agp="~{rename_contig_prefix}.agp"
        String filename_legend="~{rename_contig_prefix}_scaffolds.legend"
    }
    command <<<
        set -oeu pipefail
        grep "Version" /bbmap/README.md | sed 's/#//' 

        if [ "~{rename_contig_prefix}" != "scaffold" ]; then
            sed -e 's/scaffold/~{rename_contig_prefix}_scf/g' ~{contigs} > ~{filename_contigs}
            sed -e 's/scaffold/~{rename_contig_prefix}_scf/g' ~{scaffolds} > ~{filename_scaffolds}
            sed -e 's/scaffold/~{rename_contig_prefix}_scf/g' ~{agp} > ~{filename_agp}
            sed -e 's/scaffold/~{rename_contig_prefix}_scf/g' ~{legend} > ~{filename_legend}
        fi

    >>>

    output{
        File outcontigs = filename_contigs
        File outscaffolds = filename_scaffolds
        File outagp = filename_agp
        File outlegend = filename_legend
        File outlog = stdout()
    }
    runtime {
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
    echo -e "The reads are assembled using SPAdes(2):" >> ~{prefix}metaT_assy.info
    echo -e "`head -6 ~{spades_info} | tail -4`" >> ~{prefix}metaT_assy.info
    echo -e "An AGP file is created using fungalrelease.sh (BBTools(2)${bbtools_version})." >> ~{prefix}metaT_assy.info
    echo -e "Assembled reads are mapped using bbmap.sh (BBTools(4)${bbtools_version})." >> ~{prefix}metaT_assy.info

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
