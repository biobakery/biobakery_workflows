workflow workflowWMGX {
  String SampleName
  File FastqPair1
  File FastqPair2
    
  call kneaddata {
    input: 
    sample=SampleName, 
    file1=FastqPair1, 
    file2=FastqPair2
  }

  call metaphlan2 {
    input: 
    sample=SampleName, 
    QCFastqFile=kneaddata.QCFastqFile
  }
    
  call humann2 {
    input: 
    sample=SampleName, 
    QCFastqFile=kneaddata.QCFastqFile, 
    TaxonomicProfileFile=metaphlan2.TaxonomicProfileFile
  }
}


task kneaddata {
  File file1
  File file2
  String sample
  Int? mem_gb
  Int mem = select_first([mem_gb, 24])
  Int? preemptible_attempts_override
  Int preemptible_attempts = select_first([preemptible_attempts_override, 2])
  
  String humanDatabase = "/cromwell_root/databases/kneaddata_human/"
  String rrnaDatabase = "/cromwell_root/databases/kneaddata_rrna/"

  # download the two reference databases and then run kneaddata. We should have version control over these. 
  command {    
  	#create databases
    mkdir -p ${humanDatabase} ${rrnaDatabase}
    kneaddata_database --download human_genome bowtie2 ${humanDatabase}
    kneaddata_database --download ribosomal_RNA bowtie2 ${rrnaDatabase}
    
    #run kneaddata
    kneaddata --input ${file1} --input ${file2} --output ./ --serial \
    --reference-db ${humanDatabase} --reference-db ${rrnaDatabase} --threads 8 \
    --remove-intermediate-output --output-prefix ${sample} --cat-final-output
  }
  
  output {
    File QCFastqFile = "${sample}.fastq"
  }
  
  runtime {
    docker: "biobakery/kneaddata:latest" #This will be changed to a select first later.
    cpu: 8
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk 501 SSD"
  }
}

task metaphlan2 {
  File QCFastqFile
  String sample
  Int? mem_gb
  Int mem = select_first([mem_gb, 8])
  Int? preemptible_attempts_override
  Int preemptible_attempts = select_first([preemptible_attempts_override, 2])
  
  String tmpdir = "/cromwell_root/tmp/"

  command {
  	mkdir -p ${tmpdir}
    metaphlan2.py ${QCFastqFile} --input_type multifastq --nproc 8 --no_map --tmp_dir ${tmpdir} \
    --output_file ${sample}.relative_abundance.txt
  }
    
  output {
    File TaxonomicProfileFile = "${sample}.relative_abundance.txt"
  }
  
  runtime {
    docker: "biobakery/metaphlan2:latest"
    cpu: 8
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk 50 SSD"
  }
}

task humann2 {
  File QCFastqFile
  File TaxonomicProfileFile
  String sample
  Int? mem_gb
  Int mem = select_first([mem_gb, 24])
  Int? preemptible_attempts_override
  Int preemptible_attempts = select_first([preemptible_attempts_override, 2])
  
  String databases = "/cromwell_root/databases/"

  # download the two reference databases and run humann2
  command {
  	mkdir -p ${databases}
    humann2_databases --download chocophlan full ${databases}
    humann2_databases --download uniref uniref90_diamond ${databases}

    humann2 --input ${QCFastqFile} --output ./ --taxonomic-profile ${TaxonomicProfileFile} --threads 8 --remove-temp-output --o-log ${sample}.log
    }
    
    output {
      File GeneFamiliesFile = "${sample}_genefamilies.tsv"
      File PathwayAbundanceFile = "${sample}_pathabundance.tsv"
      File PathwayCoverageFile = "${sample}_pathcoverage.tsv"
      File LogFile = "${sample}.log"
    }
  
  runtime {
    docker:"biobakery/humann2:latest"
    cpu: 8
      memory: mem + " GB"
      disks: "local-disk 120 SSD"
      preemptible: preemptible_attempts
  }
}
