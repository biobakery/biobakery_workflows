version 1.0

workflow workflowMTX {

  # Required input variables
  input {
    File inputRead1Files
    String inputRead1Identifier
    String inputRead2Identifier
    String inputExtension
    String ProjectName
    String AdapterType
    
    # Database locations
    File versionSpecifichumanDB
    File versionSpecifictrancriptDB
    File versionSpecificrrnaDB
    File versionSpecificChocophlan
    File versionSpecificUniRef90
    File versionSpecificUtilityMapping
    
    # Optional input variables
    Boolean? bypassFunctionalProfiling
    String? dataType
    File? inputMetadataFile
    Int? preemptibleAttemptsOverride
    Int? MaxMemGB_QualityControlTasks
    Int? MaxMemGB_TaxonomicProfileTasks
    Int? MaxMemGB_FunctionalProfileTasks
    
    File? customQCDB1
    File? customQCDB2
    File? customQCDB3
  }
  
  # Set the docker tags
  String kneaddataDockerImage = "biobakery/kneaddata:0.9.0"
  String metaphlanDockerImage = "biobakery/metaphlan:3.0.1"
  String humannDockerImage = "biobakery/humann:3.0.0.a.4"
  String workflowsDockerImage = "biobakery/workflows:3.0.0.a.6_anadama0.7.9_no_metaphlan_db"
  
  # Set the default data type
  String dataTypeSetting = select_first([dataType, "mtx"])
  
  # Set if metadata is provided
  String metadataSet = if defined(inputMetadataFile) then "yes" else "no"
  
  # Set bypass mode
  Boolean setbypassFunctionalProfiling = select_first([bypassFunctionalProfiling,false])

  # Output file names to match AnADAMA2 workflow output names
  String QCReadCountFileName = "kneaddata_read_count_table.tsv"  

  String JoinedTaxonomicProfilesFileName="metaphlan_taxonomic_profiles.tsv"
  String TaxonomicProfilesCountsFileName="metaphlan_species_counts_table.tsv"
  String JoinGeneFamilesOutFileName="genefamilies.tsv"
  String JoinECsOutFileName="ecs.tsv"
  String JoinKOsOutFileName="kos.tsv"
  String JoinPathwaysOutFileName="pathabundance.tsv"

  String JoinGeneFamilesRelabOutFileName="genefamilies_relab.tsv"
  String JoinECsRelabOutFileName="ecs_relab.tsv"
  String JoinPathwaysRelabOutFileName="pathabundance_relab.tsv"

  String CountRelabGenesFileName="humann_genefamilies_relab_counts.tsv"
  String CountRelabECsFileName="humann_ecs_relab_counts.tsv"
  String CountRelabPathwaysFileName="humann_pathabundance_relab_counts.tsv"

  String JoinedFeatureCountsFileName="humann_feature_counts.tsv"
  String FunctionalCountFileName = "humann_read_and_species_count_table.tsv"
  
  String VisualizationsFileName = ProjectName+"_visualizations"
  
  # default mem settings
  Int JoinNormMemDefault = 10
  Int JoinNormMemDefaultGenes = 50

  # read in a file of the read1 paths
  Array[Array[String]] inputRead1 = read_tsv(inputRead1Files)
  
  # get the sample name and read2 file path
  scatter (read1 in inputRead1) {
     Array[String] pairSet = [read1[0], sub(read1[0], inputRead1Identifier + inputExtension, inputRead2Identifier + inputExtension), sub(basename(read1[0]), inputRead1Identifier + inputExtension, "")]
  }

  Array[Array[String]] PairPaths = pairSet

  # run tasks for each set of read pairs    
  scatter (ReadPair in PairPaths)
  {
    # Part 1: For each sample, run quality control with KneadData
    call QualityControl {
      input:
      rawfile1=ReadPair[0],
      rawfile2=ReadPair[1],
      sample=ReadPair[2],
      adapterType=AdapterType,
      humanDB=versionSpecifichumanDB,
      transcriptDB=versionSpecifictrancriptDB,
      rrnaDB=versionSpecificrrnaDB,
      customDB1=customQCDB1,
      customDB2=customQCDB2,
      customDB3=customQCDB3,
      dataType=dataTypeSetting,
      kneaddataDockerImage=kneaddataDockerImage,
      preemptibleAttemptsOverride=preemptibleAttemptsOverride,
      MaxMemGB=MaxMemGB_QualityControlTasks
    }
    # Part 2: For each sample, run taxonomic profiling with MetaPhlAn v2
    call TaxonomicProfile {
      input:
      sample=ReadPair[2], 
      QCFastqFile=QualityControl.QCFastqFile,
      metaphlanDockerImage=metaphlanDockerImage,
      preemptibleAttemptsOverride=preemptibleAttemptsOverride,
      MaxMemGB=MaxMemGB_TaxonomicProfileTasks
    }
   }
   
   if (! setbypassFunctionalProfiling ) {
    # Part 3: For each sample, run functional profiling with HUMAnN v2
    scatter (sample_index in range(length(PairPaths))) {
      call FunctionalProfile {
        input: 
        sample=PairPaths[sample_index][2], 
        QCFastqFile=QualityControl.QCFastqFile[sample_index], 
        TaxonomicProfileFile=TaxonomicProfile.TaxonomicProfileFile[sample_index],
        versionSpecificChocophlan=versionSpecificChocophlan,
        versionSpecificUniRef90=versionSpecificUniRef90,
        humannDockerImage=humannDockerImage,
        preemptibleAttemptsOverride=preemptibleAttemptsOverride,
        MaxMemGB=MaxMemGB_FunctionalProfileTasks
      }

      # regroup gene families to ECs
      call Regroup as RegroupECs {
        input:
        GeneFamiliesFile=FunctionalProfile.GeneFamiliesFile,
        versionSpecificUtilityMapping=versionSpecificUtilityMapping,
        OutFileName=PairPaths[sample_index][2]+"_ecs.tsv",
        humannDockerImage=humannDockerImage,
        groupName="uniref90_level4ec"
      }
      
      # regroup gene families to KOs
      call Regroup as RegroupKOs {
        input:
        GeneFamiliesFile=FunctionalProfile.GeneFamiliesFile,
        versionSpecificUtilityMapping=versionSpecificUtilityMapping,
        OutFileName=PairPaths[sample_index][2]+"_kos.tsv",
        humannDockerImage=humannDockerImage,
        groupName="uniref90_ko"
      }      
   
      # compute relative abundance for gene families, ecs, and pathways
      call RenormTable as RenormTableGenes {
        input:
        InFile=FunctionalProfile.GeneFamiliesFile,
        OutFileName=PairPaths[sample_index][2]+"_genefamilies_relab.tsv",
        humannDockerImage=humannDockerImage,
        MaxMemGB=JoinNormMemDefaultGenes
      }
      call RenormTable as RenormTableECs {
        input:
        InFile=RegroupECs.OutFile,
        OutFileName=PairPaths[sample_index][2]+"_ecs_relab.tsv",
        humannDockerImage=humannDockerImage,
        MaxMemGB=JoinNormMemDefault
      }
      call RenormTable as RenormTablePathways {
        input:
        InFile=FunctionalProfile.PathwayAbundanceFile,
        OutFileName=PairPaths[sample_index][2]+"_pathabundance_relab.tsv",
        humannDockerImage=humannDockerImage,
        MaxMemGB=JoinNormMemDefault
      }
    }
  
    # count the features during each alignment step
    call FunctionalCount {
      input:
      FunctionalLogFiles=FunctionalProfile.LogFile,
      OutFileName=FunctionalCountFileName,
      workflowsDockerImage=workflowsDockerImage
    }
  
    call JoinTables as JoinGeneFamilies {
      input:
      InFiles=FunctionalProfile.GeneFamiliesFile,
      OutFileName=JoinGeneFamilesOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefaultGenes
    }
    call JoinTables as JoinECs {
      input:
      InFiles=RegroupECs.OutFile,
      OutFileName=JoinECsOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinKOs {
      input:
      InFiles=RegroupKOs.OutFile,
      OutFileName=JoinKOsOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinPathways {
      input:
      InFiles=FunctionalProfile.PathwayAbundanceFile,
      OutFileName=JoinPathwaysOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinGeneFamiliesRelab {
      input:
      InFiles=RenormTableGenes.OutFile,
     OutFileName=JoinGeneFamilesRelabOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinECsRelab {
      input:
      InFiles=RenormTableECs.OutFile,
      OutFileName=JoinECsRelabOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinPathwaysRelab {
      input:
      InFiles=RenormTablePathways.OutFile,
      OutFileName=JoinPathwaysRelabOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }

    # create a file of feature counts from functional profiling data
    call CountFeatures as CountRelabGenes {
      input:
      InFile=JoinGeneFamiliesRelab.OutFile,
      OutFileName=CountRelabGenesFileName,
      Options="--reduce-sample-name --ignore-un-features --ignore-stratification",
      workflowsDockerImage=workflowsDockerImage
    }
    call CountFeatures as CountRelabECs {
      input:
      InFile=JoinECsRelab.OutFile,
      OutFileName=CountRelabECsFileName,
      Options="--reduce-sample-name --ignore-un-features --ignore-stratification",
      workflowsDockerImage=workflowsDockerImage
    }
    call CountFeatures as CountRelabPathways {
      input:
      InFile=JoinPathwaysRelab.OutFile,
      OutFileName=CountRelabPathwaysFileName,
      Options="--reduce-sample-name --ignore-un-features --ignore-stratification",
      workflowsDockerImage=workflowsDockerImage
    }
    call JoinTables as JoinFeatureCounts {
      input:
      InFiles=[CountRelabGenes.OutFile, CountRelabECs.OutFile, CountRelabPathways.OutFile],
      OutFileName=JoinedFeatureCountsFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }   
  }

  # count the reads during each step of QC
  call QCReadCount {
    input:
    LogFiles=QualityControl.LogFile,
    OutFileName=QCReadCountFileName,
    kneaddataDockerImage=kneaddataDockerImage
  }

  # count the species from the taxonomic profiles
  call CountFeatures as TaxonomicCount {
    input:
    InFile=JoinTaxonomicProfiles.OutFile,
    OutFileName=TaxonomicProfilesCountsFileName,
    Options="--include s__ --filter t__ --reduce-sample-name",
    workflowsDockerImage=workflowsDockerImage
  } 
  # join all taxonomic profiles, gene families, ecs, and pathways (including relative abundance files) from all samples
  call JoinTaxonomicProfiles {
    input:
    InFiles=TaxonomicProfile.TaxonomicProfileFile,
    OutFileName=JoinedTaxonomicProfilesFileName,
    workflowsDockerImage=workflowsDockerImage,
    MaxMemGB=JoinNormMemDefault
  }

  call VisualizationReport {
    input:
    QCCountsFile=QCReadCount.OutFile,
    TaxonomicProfileFile=JoinTaxonomicProfiles.OutFile,
    setbypassFunctionalProfiling=setbypassFunctionalProfiling,
    PathwaysFile=JoinPathwaysRelab.OutFile,
    ECsFile=JoinECsRelab.OutFile,
    FunctionalReadSpeciesCountFile=FunctionalCount.OutFile,
    FunctionalFeatureCountsFile=JoinFeatureCounts.OutFile,
    ProjectName=ProjectName,
    OutFileName=VisualizationsFileName,
    metadataSet=metadataSet,
    MetadataFile=inputMetadataFile,
    workflowsDockerImage=workflowsDockerImage
  }

}

task QualityControl {
  input {
    File rawfile1
    File rawfile2
    String sample
    String adapterType
    File humanDB
    File transcriptDB
    File rrnaDB
    File? customDB1
    File? customDB2
    File? customDB3
    String dataType
    String kneaddataDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
  }
  Int mem = select_first([MaxMemGB, 24])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])
  
  String useCustomDB1 = if defined(customDB1) then "yes" else "no"
  String useCustomDB2 = if defined(customDB2) then "yes" else "no"
  String useCustomDB3 = if defined(customDB3) then "yes" else "no"
  
  String humanDatabase = "databases/kneaddata_human/"
  String transcriptDatabase = "databases/kneaddata_transcript/"
  String rrnaDatabase = "databases/kneaddata_rrna/"
  
  # Add additional database to run options depending on data type
  String options = if dataType == "mtx" then "--reference-db ${transcriptDatabase} --reference-db ${rrnaDatabase}" else ""
  
  String customDatabase1 = "databases/db1/"
  String customDatabase2 = "databases/db2/"
  String customDatabase3 = "databases/db3/"
  String custom_options = if defined(customDB2) then "--reference-db ${customDatabase1} --reference-db ${customDatabase2}" else "--reference-db ${customDatabase1}"
  String custom_options_add = if defined(customDB3) then "--reference-db ${customDatabase3} " else " "
  

  # download the two reference databases and then run kneaddata.
  command <<< 
  
    # download second custom database if set
    if [ ~{useCustomDB3} == 'yes' ]; then
        mkdir -p ~{customDatabase3}
        tar xzvf ~{customDB3} -C ~{customDatabase3}
    fi

    # download second custom database if set
    if [ ~{useCustomDB2} == 'yes' ]; then
        mkdir -p ~{customDatabase2}
        tar xzvf ~{customDB2} -C ~{customDatabase2}
    fi
  
    # use custom databases if provided instead of reference
    if [ ~{useCustomDB1} == 'yes' ]; then
        mkdir -p ~{customDatabase1}
        tar xzvf ~{customDB1} -C ~{customDatabase1}
        
        #run kneaddata with custom databases
        kneaddata --input ~{rawfile1} --input ~{rawfile2} --output ./ --serial \
        --threads 8 --output-prefix ~{sample} --cat-final-output --run-fastqc-start ~{custom_options} ~{custom_options_add} --sequencer-source ~{adapterType}
    fi
    
    if [ ~{useCustomDB1} == 'no' ]; then
        # download the human reference
        mkdir -p ~{humanDatabase}
        kneaddata_database --download human_genome bowtie2 ~{humanDatabase} --database-location ~{humanDB}
    
        # if data is of type mtx, then download additional database
        if [ ~{dataType} == 'mtx' ]; then
            #create databases
            mkdir -p ~{transcriptDatabase}
            kneaddata_database --download human_transcriptome bowtie2 ~{transcriptDatabase} --database-location ~{transcriptDB}

            mkdir -p ~{rrnaDatabase}
            kneaddata_database --download ribosomal_RNA bowtie2 ~{rrnaDatabase} --database-location ~{rrnaDB}
        fi
    
        #run kneaddata with two reference databases
        kneaddata --input ~{rawfile1} --input ~{rawfile2} --output ./ --serial --reference-db ~{humanDatabase} \
        --threads 8 --output-prefix ~{sample} --cat-final-output --run-fastqc-start ~{options} --sequencer-source ~{adapterType}
    fi
    
    # gzip outputs to save space
    gzip *.fastq
  >>>
  
  output {
    File QCFastqFile = "${sample}.fastq.gz"
    File LogFile = "${sample}.log"
    Array[File] ContaminateReads = glob("*contam*.fastq.gz") # Keep the intermediate contaminate sequences
    Array[File] FastQCOutputsZip = glob("fastqc/*.zip") # Keep the fastqc output files (zip)
    Array[File] FastQCOutputsHtml = glob("fastqc/*.html") # Keep the fastqc output files (html)
  }
  
  runtime {
    docker: kneaddataDockerImage
    cpu: 8
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk 501 SSD"
  }
}

task TaxonomicProfile {
  input {
    File QCFastqFile
    String sample
    String metaphlanDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
  }
  
  Int mem = select_first([MaxMemGB, 8])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])
  
  String tmpdir = "tmp/"

  # create a temp directory and then run metaphlan
  command {
    mkdir -p ${tmpdir}
    metaphlan ${QCFastqFile} --input_type fastq --nproc 8 --no_map --tmp_dir ${tmpdir} \
    --output_file ${sample}.tsv
  }
    
  output {
    File TaxonomicProfileFile = "${sample}.tsv"
  }
  
  runtime {
    docker: metaphlanDockerImage
    cpu: 8
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk 50 SSD"
  }
}

task FunctionalProfile {

  input {
    File QCFastqFile
    File TaxonomicProfileFile
    String sample
    File versionSpecificChocophlan
    File versionSpecificUniRef90
    String humannDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
  }

  Int mem = select_first([MaxMemGB, 32])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])
  
  String databases = "databases/"

  # download the two reference databases and run humann
  command {
    mkdir -p ${databases}
    humann_databases --download chocophlan full ${databases} --database-location ${versionSpecificChocophlan}
    humann_databases --download uniref uniref90_diamond ${databases} --database-location ${versionSpecificUniRef90}

    humann --input ${QCFastqFile} --output ./ --taxonomic-profile ${TaxonomicProfileFile} --threads 8 --o-log ${sample}.log
    }
    
    output {
      File GeneFamiliesFile = "${sample}_genefamilies.tsv"
      File PathwayAbundanceFile = "${sample}_pathabundance.tsv"
      File PathwayCoverageFile = "${sample}_pathcoverage.tsv"
      File LogFile = "${sample}.log"
      Array[File] UnalignedReads = glob("${sample}_humann_temp/*.fa") # Keep the unaligned reads after each mapping step
    }

  runtime {
    docker: humannDockerImage
    cpu: 8
      memory: mem + " GB"
      disks: "local-disk 120 SSD"
      preemptible: preemptible_attempts
  }
}
task Regroup {
  input {
    File GeneFamiliesFile
    File versionSpecificUtilityMapping
    String OutFileName
    String humannDockerImage
    String groupName
  }

  String databases = "databases/"

  # download the utility databases and regroup to ECs
  command {
    mkdir -p ${databases}
    humann_databases --download utility_mapping full ${databases} --database-location ${versionSpecificUtilityMapping}

    humann_regroup_table --input ${GeneFamiliesFile} --output ${OutFileName} --groups ${groupName}
  }
    
  output {
      File OutFile = "${OutFileName}"
  }

  runtime {
    docker: humannDockerImage
    cpu: 1
      memory: "10 GB"
      disks: "local-disk 20 SSD"
  }
}

task RenormTable {
  input {
    File InFile
    String OutFileName
    String humannDockerImage
    Int? MaxMemGB
  }
  
  Int mem = select_first([MaxMemGB, 10])

  String databases = "databases/"

  # download the utility databases and renorm tables to relative abundance
  command {
    humann_renorm_table --input ${InFile} --output ${OutFileName} --units relab --special n
  }

  output {
      File OutFile = "${OutFileName}"
  }

  runtime {
    docker: humannDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 20 SSD"
  }
}

task QCReadCount {
  input {
    Array[File] LogFiles
    String OutFileName
    String kneaddataDockerImage
  }
  
  # symlink input files to working directory
  # create a table of read counts during each qc step
  command {
    for logfile in ${sep=' ' LogFiles}; do ln -s $logfile; done
    
    kneaddata_read_count_table --input ./ --output ${OutFileName}
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: kneaddataDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}

task JoinTaxonomicProfiles {
  input {
    Array[File] InFiles
    String OutFileName
    String workflowsDockerImage
    Int? MaxMemGB
  }
  
  Int mem = select_first([MaxMemGB, 10])

  # symlink input files to working directory
  # join all files into a single file
  command {
    for infile in ${sep=' ' InFiles}; do ln -s $infile; done
    
    join_taxonomic_profiles.py --input ./ --output ${OutFileName} --file_name .tsv
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 10 SSD"
  }
}

task JoinTables {
  input {
    Array[File] InFiles
    String OutFileName
    String humannDockerImage
    Int? MaxMemGB
  }
  
  Int mem = select_first([MaxMemGB, 10])

  # symlink input files to working directory
  # join all files into a single file
  command {
    for infile in ${sep=' ' InFiles}; do ln -s $infile; done
    
    humann_join_tables --input ./ --output ${OutFileName} --file_name .tsv
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: humannDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 10 SSD"
  }
}

task FunctionalCount {
  input {
    Array[File] FunctionalLogFiles
    String OutFileName
    String workflowsDockerImage
  }

  # symlink logs to working directory
  # compile alignment counts from humann logs
  command {
    for infile in ${sep=' ' FunctionalLogFiles}; do ln -s $infile; done
  
    get_counts_from_humann_logs.py --input ./ --output ${OutFileName}
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}

task CountFeatures {
  input {
    File InFile
    String OutFileName
    String Options
    String workflowsDockerImage
  }

  # count features in a table
  command {
    count_features.py --input ${InFile} --output ${OutFileName} ${Options}
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}

task VisualizationReport {
  input {
    File QCCountsFile
    File TaxonomicProfileFile
    Boolean setbypassFunctionalProfiling
    File? PathwaysFile
    File? FunctionalReadSpeciesCountFile
    File? FunctionalFeatureCountsFile
    File? ECsFile
    String ProjectName
    String OutFileName
    String metadataSet
    File? MetadataFile
    String workflowsDockerImage
  }
  
  String QCCountsFolder = "input/kneaddata/merged/"
  String TaxonomyFolder = "input/metaphlan/merged/"
  String FunctionalMergedFolder = "input/humann/merged/"
  String FunctionalCountsFolder = "input/humann/counts/"
  
  # symlink files to the expected folder locations
  # run visualizations
  command <<<
    mkdir -p ~{QCCountsFolder}
    (cd ~{QCCountsFolder} && ln -s ~{QCCountsFile})
    
    mkdir -p ~{TaxonomyFolder}
    (cd ~{TaxonomyFolder} && ln -s ~{TaxonomicProfileFile})
    
    if [ ~{setbypassFunctionalProfiling} == false ]; then
      mkdir -p ~{FunctionalMergedFolder}
      (cd ~{FunctionalMergedFolder} && ln -s ~{PathwaysFile} && ln -s ~{ECsFile})
      mkdir -p ~{FunctionalCountsFolder}
      (cd ~{FunctionalCountsFolder} && ln -s ~{FunctionalReadSpeciesCountFile} && ln -s ~{FunctionalFeatureCountsFile})
    fi

    if [ ~{metadataSet} == 'yes' ]; then
      biobakery_workflows wmgx_vis --input input --output ~{OutFileName} --project-name ~{ProjectName} --exclude-workflow-info --input-metadata ~{MetadataFile}
    else
      biobakery_workflows wmgx_vis --input input --output ~{OutFileName} --project-name ~{ProjectName} --exclude-workflow-info 
    fi
    
    >>>
  
  output {
    File OutFile = "${OutFileName}.zip"
  }
  
  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}
