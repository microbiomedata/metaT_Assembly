version 1.0

import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/raw/main1.0/metatranscriptome/metatranscriptome_assy_rnaspades.wdl" as http_rnaspades
import "https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/raw/main1.0/common/mapping.wdl" as mapping


workflow metatranscriptome_assy {
    input{
        Array[File] input_files
        String rename_contig_prefix="scaffold"
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



    output {
        File final_tar_bam = tar_bams.outtarbam
        File final_contigs = create_agp.outcontigs  # annotation.input_file
        File final_scaffolds = create_agp.outscaffolds
        File final_log = assy.log
        File final_readlen = readstats_raw.outreadlen
        File final_sam = finalize_bams.outsam
        File final_bam = finalize_bams.outbam
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
    }
    command <<<

        if [ "~{rename_contig_prefix}" != "scaffold" ]; then
            sed -i 's/scaffold/~{rename_contig_prefix}_scf/g' ~{contigs} ~{scaffolds} ~{agp} ~{legend}
        fi

    >>>

    output{
        File outcontigs = contigs
        File outscaffolds = scaffolds
        File outagp = agp
        File outlegend = legend
    }
    runtime {
        docker: container
    }
}

