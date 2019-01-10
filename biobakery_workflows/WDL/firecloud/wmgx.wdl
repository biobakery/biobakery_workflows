workflow workflowWMGX {

  File InputFileSampleList
  String FastqPair1Extension
  String FastqPair2Extension

  # set working directory (for output and temp files) to location
  # with extra local disk space
  String WorkingDirectory = "/cromwell_root/"

  Map[String, String] SamplesPaths = read_map(InputFileSampleList)

  # loop through each sample listed in input file sample list
  scatter (row in SamplesPaths){

    String SampleDir = row.right
    String SampleName = row.left

    File RawFastqPair1File = SampleDir + SampleName + FastqPair1Extension
    File RawFastqPair2File = SampleDir + SampleName + FastqPair2Extension
    
    call kneaddata {
      input: 
      sample=SampleName, 
      file1=RawFastqPair1File, 
      file2=RawFastqPair2File,
      workingdir=WorkingDirectory
    }

    call metaphlan {
      input: 
      sample=SampleName, 
      QCFastqFile=kneaddata.QCFastqFile,
      workingdir=WorkingDirectory
    }
    
    call humann2 {
      input: 
      sample=SampleName, 
      QCFastqFile=kneaddata.QCFastqFile, 
      TaxonomicProfileFile=metaphlan.TaxonomicProfileFile,
      workingdir=WorkingDirectory
    }
  }
}


task kneaddata {
  File file1
  File file2
  String sample
  String workingdir

  String human_database = workingdir + "/databases/kneaddata_human/"
  String rrna_database = workingdir + "/databases/kneaddata_rrna/"

  # download the two reference databases and then run kneaddata
  command {    
    kneaddata_database --download human_genome bowtie2 ${human_database}

    kneaddata_database --download ribosomal_RNA bowtie2 ${rrna_database}

    kneaddata --input ${file1} --input ${file2} -o ${workingdir} --serial \
    --reference-db ${human_database} --reference-db ${rrna_database} -t 8 \
    --remove-intermediate-output --output-prefix ${sample} --cat-final-output
  }
  
  output {
    File QCFastqFile = workingdir + "${sample}.fastq"
  }
  
  runtime {
    docker: "biobakery/kneaddata:latest"
    cpu: 8
      memory: "24GB"
      preemptible: 2
      disks: "local-disk 501 SSD"
  }
}

task metaphlan {
  File QCFastqFile
  String sample
  String workingdir

  String tmpdir = workingdir + "/tmp/"

  command {
    metaphlan2.py ${QCFastqFile} --input_type multifastq --nproc 8 --no_map --tmp_dir ${tmpdir} \
    --output_file ${workingdir}/${sample}.relative_abundance.txt
  }
    
  output {
    File TaxonomicProfileFile = "${workingdir}/${sample}.relative_abundance.txt"
  }
  
  runtime {
    docker: "biobakery/metaphlan2:latest"
    cpu: 8
      memory: "8GB"
      preemptible: 2
      disks: "local-disk 50 SSD"
  }

}

task humann2 {
  File QCFastqFile
  File TaxonomicProfileFile
  String sample
  String workingdir
  
  String databases = workingdir + "/databases/"

  # download the two reference databases and run humann2
  command {      
    humann2_databases --download chocophlan full ${databases}
    humann2_databases --download uniref uniref90_diamond ${databases}

    humann2 --input ${QCFastqFile} --output ${workingdir} --taxonomic-profile ${TaxonomicProfileFile} --threads 8 --remove-temp-output --o-log ${sample}.log
    }
    
    output {
      File GeneFamiliesFile = "${workingdir}/${sample}_genefamilies.tsv"
      File PathwayAbundanceFile = "$${workingdir}/{sample}_pathabundance.tsv"
      File PathwayCoverageFile = "${workingdir}/${sample}_pathcoverage.tsv"
      File LogFile = "${workingdir}/${sample}.log"
    }
  
  runtime {
    docker:"biobakery/humann2:latest"
    cpu: 8
      memory: "24GB"
      disks: "local-disk 120 SSD"
      preemptible: 2
  }
}

