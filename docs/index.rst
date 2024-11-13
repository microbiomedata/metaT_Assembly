MetaT Assembly Workflow (v0.0.2)
=============================

.. image:: mt_assy_workflow2024.png
   :align: center
   :scale: 50%


Workflow Overview
-----------------

This workflow was developed by Brian Foster at JGI. Original repo can be found (here)[https://code.jgi.doe.gov/BFoster/jgi_meta_wdl/-/tree/master/metatranscriptome]. This workflow uses SPAdes and `bbmap` to assemble and map QC'ed transcriptomic reads, with an AGP file created using `fungalrelease.sh`. 

Workflow Availability
---------------------

The workflow from GitHub uses all the listed docker images to run all third-party tools.
The workflow is available in GitHub: https://github.com/microbiomedata/metaT_Assembly; the corresponding
Docker images are available in DockerHub: 
- `microbiomedata/bbtools:38.96 <https://hub.docker.com/r/microbiomedata/bbtools>`_
- `bryce911/spades:3.15.2 <https://hub.docker.com/r/bryce911/spades>`_
- `microbiomedata/workflowmeta:1.1.1 <https://hub.docker.com/r/microbiomedata/workflowmeta>`_


Requirements for Execution 
--------------------------

(recommendations are in **bold**) 

- WDL-capable Workflow Execution Tool (**Cromwell**)
- Container Runtime that can load Docker images (**Docker v2.1.0.3 or higher**) 

Hardware Requirements
---------------------

- Memory: >120 GB RAM


Workflow Dependencies
---------------------

Third party software (This is included in the Docker image.)  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `BBTools v38.96 <https://jgi.doe.gov/data-and-tools/bbtools/>`_ (License: `BSD-3-Clause-LBNL <https://bitbucket.org/berkeleylab/jgi-bbtools/src/master/license.txt>`_)
- `SPAdes v3.15.4 <https://github.com/ablab/spades>`_ (License: `SPAdes team <https://github.com/ablab/spades?tab=License-1-ov-file#License-1-ov-file>`_)



Inputs
------

A JSON file containing the following information: 

1.	the path to the cleaned fastq file 
2.  input_interleaved (boolean)
3.  output file prefix
4.	(optional) parameters for memory 
5.	(optional) number of threads requested

An example input JSON file is shown below:

.. code-block:: JSON

    {
      "metatranscriptome_assy.input_files":["https://portal.nersc.gov/cfs/m3408/test_data/metaT/SRR11678315/readsqc_output/SRR11678315-int-0.1_filtered.fastq.gz"],
      "metatranscriptome_assy.proj_id":"SRR11678315-int-0.1"
    }


Output
------

In the outputs directory will be assembled contigs and scaffolds in fasta format from SPAdes. From `bbmap` will be mapped BAM and SAM files, including coverage, index, statistics, and a `.tar` collection. The log files, run information, and data statistics will also be included. 


An example output JSON file (scaffold_stats.json) is shown below:
   
.. code-block:: JSON 
    
    {
      "scaffolds": 240209,
      "contigs": 240566,
      "scaf_bp": 75371840,
      "contig_bp": 75367910,
      "gap_pct": 0.00521,
      "scaf_N50": 86540,
      "scaf_L50": 294,
      "ctg_N50": 87194,
      "ctg_L50": 293,
      "scaf_N90": 205805,
      "scaf_L90": 227,
      "ctg_N90": 205578,
      "ctg_L90": 227,
      "scaf_logsum": 30754,
      "scaf_powsum": 3220.307,
      "ctg_logsum": 30715,
      "ctg_powsum": 3216.241,
      "asm_score": 2.568,
      "scaf_max": 4382,
      "ctg_max": 4382,
      "scaf_n_gt50K": 0,
      "scaf_l_gt50K": 0,
      "scaf_pct_gt50K": 0,
      "gc_avg": 0.59142,
      "gc_std": 0.08812
    }



Below is an example of all the output directory files with descriptions to the right.

==================================== ============================================================================
FileName                              Description
==================================== ============================================================================
prefix_contigs.fna                    assembled FASTA contigs
prefix_scaffolds.fna                  assembled FASTA scaffolds
prefix_pairedMapped.bam               reads mapping back to the final assembly bam file
prefix_pairedMapped.sam.gz            reads mapping back to the final assembly sam.gz file
prefix_pairedMapped_sorted.bam.bai    reads mapping back to the final assembly sorted bam index file
prefix_pairedMapped_sorted.bam.cov    reads mapping back to the final assembly sorted bam coverage file 
prefix_bamfiles.tar                   collection of bam files 
prefix_scaffold_stats.json            scaffold coverage information 
prefix_readlen.txt                    read length information 
prefix_assy.info                      assembly workflow information
prefix_spades.log                     SPAdes workflow log 
==================================== ============================================================================


Version History
---------------

- 0.0.2 (release date **07/25/2024**)


Point of contact
----------------

- Original author: Brian Foster <bfoster@lbl.gov>

- Package maintainers: Chienchi Lo <chienchi@lanl.gov>
