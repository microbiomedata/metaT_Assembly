workflow metat_interleave_asm {
    String? memory
    String? threads
    File input_file
    String proj
    String rename_contig_prefix="transcript"
    String bbtools_container="microbiomedata/bbtools:38.96"
    # String spades_container="staphb/spades:3.15.5"
    String spades_container="microbiomedata/spades:3.15.0"

    call readstats_raw {
    	 input: reads_files=[input_file], container=bbtools_container
    }

    call assy {
         input: infile=input_file, container=spades_container, threads=threads
    }

    call mappingtask {
         input: reads=input_file, reference=assy.out, container=bbtools_container
    }

    call create_agp {
         input: contigs_in=assy.out, container=bbtools_container}


    output {
        File contig=create_agp.outcontigs
        File agp=create_agp.outscaffolds
        File final_readlen=readstats_raw.outreadlen 
    }
 
    meta {
        author: "Migun Shakya, B-GEN, LANL"
        email: "migun@lanl.gov"
        version: "1.0.0"
    }

}



task readstats_raw {
     Array[File] reads_files
     String container
     String single = if (length(reads_files) == 1 ) then "1" else "0"
     String reads_input="reads.input.fastq.gz"
     String outfile="readlen.txt"
    runtime {
    docker: container
    memory: "120 GiB"
    cpu:  16
    }
     command {
        if [ ${single} == 0 ]
	then
	    cat ${sep = " " reads_files } > ${reads_input}
	else
	    ln -s ${reads_files[0]} ./${reads_input}
	fi

        readlength.sh in=${reads_input} 1>| ${outfile}
     }
     output {
         File outreadlen = outfile
      }
}


task create_agp {
    File contigs_in
    String container
    String java="-Xmx48g"
    String prefix="assembly"
    String filename_contigs="${prefix}.transcripts.fasta"
    String filename_scaffolds="${prefix}.scaffolds.fasta"
    String filename_agp="${prefix}.agp"
    String filename_legend="${prefix}.scaffolds.legend"
    runtime {
    docker: container
    memory: "120 GiB"
    cpu:  16
    }

    command{
        fungalrelease.sh ${java} in=${contigs_in} out=${filename_scaffolds} outc=${filename_contigs} agp=${filename_agp} legend=${filename_legend} mincontig=200 minscaf=200 sortscaffolds=t sortcontigs=t overwrite=t
  }
    output{
	File outcontigs = filename_contigs
	File outscaffolds = filename_scaffolds
	File outagp = filename_agp
    	File outlegend = filename_legend
    }
}



task mappingtask {
    File reads
    File reference
    String container
    String? threads
    String java="-Xmx100g"
    String filename_unsorted="pairedMapped.bam"
    String filename_sorted="pairedMapped_sorted.bam"
    String dollar="$"
    String system_cpu="$(grep \"model name\" /proc/cpuinfo | wc -l)"
    String jvm_threads=select_first([threads,system_cpu])
    runtime { docker: container}

    command{
    bbmap.sh ${java} threads=${jvm_threads}  nodisk=true \
    interleaved=true ambiguous=random rgid=filename \
    in=${reads} ref=${reference} out=${filename_unsorted}
    samtools sort -m6G -@ ${jvm_threads} ${filename_unsorted} -o ${filename_sorted};
	touch ${filename_sorted};
  }
  output{
      File outbamfile = filename_sorted
   }
}




task assy {
     File infile
     String container
     String? threads
     String outprefix="rna_spades"
     String filename_outfile="${outprefix}/transcripts.fasta"
     String filename_outfile_opts="${outprefix}/params.txt"
     String filename_spadeslog ="${outprefix}/spades.log"
     String system_cpu="$(grep \"model name\" /proc/cpuinfo | wc -l)"
     String spades_cpu=select_first([threads,system_cpu])
     runtime {
            docker: container
            memory: "120 GiB"
        cpu:  16
     }
     command{
        rnaspades.py -o ${outprefix} -t ${spades_cpu} -s ${infile}
        
     }
     output {
            File out = filename_outfile
            File opt = filename_outfile_opts
            File outlog = filename_spadeslog
     }
}
