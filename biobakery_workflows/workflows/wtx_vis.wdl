version 1.0

workflow workflowMTXVis {

  # Required input variables
  input {
    File QCCountsFile
    File TaxonomicProfileFile
    File PathwaysFile
    File FunctionalReadSpeciesCountFile
    File FunctionalFeatureCountsFile
    String ProjectName
    
    # Optional input variables
    File? inputMetadataFile
  }
  
  # Set if metadata is provided
  String metadataSet = if defined(inputMetadataFile) then "yes" else "no"
  
  String VisualizationsFileName = ProjectName+"_visualizations"
  
  # generate a visualization report from the joined output files
  call VisualizationReport {
    input:
    QCCountsFile=QCCountsFile,
    TaxonomicProfileFile=TaxonomicProfileFile,
    PathwaysFile=PathwaysFile,
    FunctionalReadSpeciesCountFile=FunctionalReadSpeciesCountFile,
    FunctionalFeatureCountsFile=FunctionalFeatureCountsFile,
    ProjectName=ProjectName,
    OutFileName=VisualizationsFileName,
    metadataSet=metadataSet,
    MetadataFile=inputMetadataFile
  }
}

task VisualizationReport {
  input {
    File QCCountsFile
    File TaxonomicProfileFile
    File PathwaysFile
    File FunctionalReadSpeciesCountFile
    File FunctionalFeatureCountsFile
    String ProjectName
    String OutFileName
    String metadataSet
    File? MetadataFile
  }
  
  String QCCountsFolder = "input/kneaddata/merged/"
  String TaxonomyFolder = "input/metaphlan2/merged/"
  String FunctionalMergedFolder = "input/humann2/merged/"
  String FunctionalCountsFolder = "input/humann2/counts/"
  
  # symlink files to the expected folder locations
  # run visualizations
  command <<<
    mkdir -p ~{QCCountsFolder}
    (cd ~{QCCountsFolder} && ln -s ~{QCCountsFile})
    
    mkdir -p ~{TaxonomyFolder}
    (cd ~{TaxonomyFolder} && ln -s ~{TaxonomicProfileFile})
    
    mkdir -p ~{FunctionalMergedFolder}
    (cd ~{FunctionalMergedFolder} && ln -s ~{PathwaysFile})
    
    mkdir -p ~{FunctionalCountsFolder}
    (cd ~{FunctionalCountsFolder} && ln -s ~{FunctionalReadSpeciesCountFile} && ln -s ~{FunctionalFeatureCountsFile})

    if [ ~{metadataSet} == 'yes' ]; then
      echo "Running with metadata"
      biobakery_workflows wmgx_vis --input input --output ~{OutFileName} --project-name ~{ProjectName} --exclude-workflow-info --input-metadata ~{MetadataFile}
    else
      echo "Running without metadata"
      echo ~{metadataSet}
      biobakery_workflows wmgx_vis --input input --output ~{OutFileName} --project-name ~{ProjectName} --exclude-workflow-info 
    fi
    
    >>>
  
  output {
    File OutFile = "${OutFileName}.zip"
  }
  
  runtime {
    docker:"biobakery/workflows:0.13.5_cloud_r2"
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}
