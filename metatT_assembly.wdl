# version 1.0

import "jgi_meta_wdl/metatranscriptome/metatranscriptome_assy_rnaspades.wdl" as http_rnaspades
import "jgi_meta_wdl/common/mapping.wdl" as mapping

# import "https://raw.githubusercontent.com/microsoft/five-dollar-genome-analysis-pipeline-azure/az1.1.0/structs/GermlineStructs.wdl" as http_rnaspades
#HELPDOC
# https://cromwell.readthedocs.io/en/stable/Imports/ for import

# example of version 1.0
# https://github.com/vgteam/vg_wdl/blob/master/workflows/giraffe.wdl

# testing using miniwdl
# miniwdl run -i input_test2.json jgi_git_assembly.wdl

#TODO
# there is an issue with the repo not being downloaded from the link also it seems when you do 
#  `miniwdl check jgi_git_assembly.wdl`
# it seem to throw warnings for the script in http.
# A tool called mgiht have to do git submodule first https://git-scm.com/book/en/v2/Git-Tools-Submodules#:~:text=Git%20addresses%20this%20issue%20using,and%20keep%20your%20commits%20separate.

workflow metatranscriptome_assy {
        Array[File] input_files
        String bbtools_container = "bryce911/bbtools:38.86"
        String spades_container_prod = "bryce911/spades:3.15.4"
    

    call http_rnaspades.readstats_raw {
    	 input: reads_files=input_files, container=bbtools_container
    }

    call http_rnaspades.assy {
         input: reads_files=input_files, container=spades_container_prod
    }
    call http_rnaspades.create_agp {
         input: contigs_in=assy.out, container=bbtools_container
    }

    call mapping.mappingtask as single_run {
           input: reads=input_files[0], reference=create_agp.outcontigs, container=bbtools_container
       }
    call mapping.finalize_bams{
        	input: insing=single_run.outbamfile, container=bbtools_container
    	}

    call mapping.tar_bams as tar_bams {
    		input: insing=single_run.outbamfile,  container=bbtools_container
    	}



    output {
        File final_tar_bam = tar_bams.filename_tarbam
        File final_contigs = create_agp.outcontigs
        File final_scaffolds = create_agp.outscaffolds
        File final_log = assy.log
	    File final_readlen = readstats_raw.outreadlen
    }

        meta {
        author: "Migun Shakya, B-GEN, LANL"
        email: "migun@lanl.gov"
        version: "0.0.1"
    }
}

