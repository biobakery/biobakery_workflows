workflow workflowBiobakery {

	File Sample_Path_List
	String Fastq1_Extension
	String Fastq2_Extension
	Map[String, String] SamplesPaths = read_map(Sample_Path_List)
	scatter (pair in SamplesPaths){

		String SampleDir = pair.right
		String sample = pair.left

		File F1 = SampleDir + sample + Fastq1_Extension
		File F2 = SampleDir + sample + Fastq2_Extension
		call qcAdapters {
			input: 
			sample=sample, 
			file1=F1, 
			file2=F2
		}
		
		call qcQualityHuman {
			input: 
			sample=sample, 
			file1=qcAdapters.fileR1, 
			file2=qcAdapters.fileR2
		}

		call metaphlan {
			input: 
			sample=sample, 
			r1=qcQualityHuman.fileR1, 
			r2=qcQualityHuman.fileR2, 
			s1=qcQualityHuman.fileS1, 
			s2=qcQualityHuman.fileS2,
			sample=sample		
		}
		
		call humann2 {
			input: 
			sample=sample, 
			r1=qcQualityHuman.fileR1, 
			r2=qcQualityHuman.fileR2, 
			s1=qcQualityHuman.fileS1, 
			s2=qcQualityHuman.fileS2,
			taxProfile=metaphlan.fileProfile,
			sample=sample		
		}

		call regroupHumann2 {
			input:
			geneFamilies=humann2.fileGeneFamilies,
			sample=sample
		}

		call strainphlanMarkers {
			input:
			sampleBAM=metaphlan.fileBam,
			sample=sample
		}
	}

	call combineMetaphlan {
		input:
		taxProfiles=metaphlan.fileProfile
	}

	call combineHumann2 {
		input:
		filesGeneFamilies=humann2.fileGeneFamilies,
		filesPathways=humann2.filePathwayAbundamce,
		filesEC=regroupHumann2.fileECAbundance,
		filesKO=regroupHumann2.fileKOAbundance,
		filesEggnog=regroupHumann2.fileEggnogAbundance,
		filesGO1000=regroupHumann2.fileGO1000Abundance,
		filesPfam=regroupHumann2.filePfamAbundance

	}

	call strainphlanClades {
		input:
		filesSampleMarkers=strainphlanMarkers.sampleMarker
	}

	Array[String] Clades = read_lines(strainphlanClades.cladeList)
	scatter (clade in Clades){
		
		call strainphlanTree {
			input:
			filesSampleMarkers=strainphlanMarkers.sampleMarker,
			clade=clade
		}
	}
}


task qcAdapters {
	File file1
	File file2
	String sample

	command {

		# move file into name that fits with the sample naming strategy
		mv ${file1} ${sample}.1.fq.gz
		mv ${file2} ${sample}.2.fq.gz

		trim_galore --paired --phred33 --quality 0 --stringency 5 --length 10 \
		${sample}.1.fq.gz ${sample}.2.fq.gz

		mv ${sample}.1_val_1.fq.gz ${sample}.adapterTrimmed.1.fq.gz
		mv ${sample}.2_val_2.fq.gz ${sample}.adapterTrimmed.2.fq.gz
	}
	
	output {
		File fileR1 = "${sample}.adapterTrimmed.1.fq.gz"
		File fileR2 = "${sample}.adapterTrimmed.2.fq.gz"
	}
	
	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 1
  		memory: "1GB"
  		preemptible: 2
  		disks: "local-disk 40 SSD"
	}
}

task qcQualityHuman {
	File file1
	File file2
	String sample
	File ref1
	File ref2
	File ref3
	File ref4
	File ref5
	File ref6

	command {
		kneaddata --input ${file1} --input ${file2} -o . \
		-db tools-rx/DATABASES/HG19 --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:4:15 MINLEN:50" -t 4
		rm *trimmed*
		rm *bowtie2*
		
		gzip ${sample}.adapterTrimmed.1_kneaddata_paired_1.fastq
		gzip ${sample}.adapterTrimmed.1_kneaddata_paired_2.fastq
		gzip ${sample}.adapterTrimmed.1_kneaddata_unmatched_1.fastq
		gzip ${sample}.adapterTrimmed.1_kneaddata_unmatched_2.fastq
	}
	
	output {
		File fileR1 = "${sample}.adapterTrimmed.1_kneaddata_paired_1.fastq.gz"
		File fileR2 = "${sample}.adapterTrimmed.1_kneaddata_paired_2.fastq.gz"
		File fileS1 = "${sample}.adapterTrimmed.1_kneaddata_unmatched_1.fastq.gz"
		File fileS2 = "${sample}.adapterTrimmed.1_kneaddata_unmatched_2.fastq.gz"
	}
	
	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 4
  		memory: "24GB"
  		preemptible: 2
  		disks: "local-disk 501 SSD"
	}
}

task metaphlan {
	File r1
	File r2
	File s1
	File s2
	String sample
	File ref1
	File ref2
	File ref3
	File ref4
	File ref5
	File ref6
	File ref7
	File ref8
	File ref9

	command {
    	ls -l /cromwell_root/tools-rx/DATABASES/METAPHLAN2/db_v20/
    	
    	zcat ${r1} ${r2} ${s1} ${s2} | metaphlan2.py --bowtie2db /cromwell_root/tools-rx/DATABASES/METAPHLAN2/db_v20 --index v20_m200 --mpa_pkl ${ref1} --input_type multifastq -t rel_ab  --bt2_ps "very-sensitive" --tmp_dir /tmp --ignore_viruses -s ${sample}.sam --bowtie2out ${sample}.bowtie2.out > ${sample}.relative_abundance.txt

        # Copy bam file
        samtools view -bS ${sample}.sam -o ${sample}.bam
    }
    
    output {
    	File fileProfile = "${sample}.relative_abundance.txt"
    	File fileBam = "${sample}.bam"
    	File fileBowtie = "${sample}.bowtie2.out"
    }
	
	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 1
  		memory: "8GB"
  		preemptible: 2
  		disks: "local-disk 50 SSD"
	}

}

task humann2 {
	File r1
	File r2
	File s1
	File s2
	File taxProfile
	File refChocophlan
	File refUniref90
	String sample
	
	command {    	
    	# get the chocophlan database
    	tar -zxvf ${refChocophlan}

    	#prepare fastqs
    	zcat ${r1} ${r2} ${s1} ${s2} > ${sample}.fq
    	#prepare output folder
    	mkdir humann2_run
    	humann2 -i ${sample}.fq -o humann2_run --taxonomic-profile ${taxProfile} --threads 8 --translated-alignment diamond --search-mode uniref90 --gap-fill on --nucleotide-database chocophlan --protein-database /cromwell_root/tools-rx/DATABASES/HUMANN2/UNIREF --bowtie2 /appdownload/bowtie2-2.3.4.1-linux-x86_64 
    }
    
    output {
    	File fileGeneFamilies = "humann2_run/${sample}_genefamilies.tsv"
    	File filePathwayAbundamce = "humann2_run/${sample}_pathabundance.tsv"
    	File filePathwayCoverage = "humann2_run/${sample}_pathcoverage.tsv"
    	File fileLog = "humann2_run/${sample}_humann2_temp/${sample}.log"
    }
	
	runtime {
		docker:"gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 8
  		memory: "24GB"
  		disks: "local-disk 120 SSD"
  		preemptible: 2
	}

}

task regroupHumann2 {
	File geneFamilies
	File refEC
	File refKO
	File refPfam
	File refEggnog
	File refGO1000
	String sample

	command {
		humann2_regroup_table -i ${geneFamilies} -c ${refEC} -o ${sample}_ec.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refKO} -o ${sample}_ko.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refPfam} -o ${sample}_pfam.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refEggnog} -o ${sample}_eggnog.tsv
		humann2_regroup_table -i ${geneFamilies} -c ${refGO1000} -o ${sample}_go1000.tsv

	}

	output {
		File fileECAbundance = "${sample}_ec.tsv"
		File fileKOAbundance = "${sample}_ko.tsv"
		File filePfamAbundance = "${sample}_pfam.tsv"
		File fileEggnogAbundance = "${sample}_eggnog.tsv"
		File fileGO1000Abundance = "${sample}_go1000.tsv"

	}	

	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 1
  		memory: "4GB"
  		preemptible: 2
  		disks: "local-disk 50 SSD"
	}

}

task combineMetaphlan {
	Array[File] taxProfiles

	command {
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " taxProfiles} > metaphlan_merged_results.tsv
	}

	output {
		File metaphlanMerged = "metaphlan_merged_results.tsv"
	}

	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 1
  		memory: "4GB"
  		preemptible: 2
  		disks: "local-disk 50 SSD"
	}


}

task combineHumann2 {
	Array[File] filesGeneFamilies
	Array[File] filesPathways
	Array[File] filesEC
	Array[File] filesKO
	Array[File] filesEggnog
	Array[File] filesGO1000
	Array[File] filesPfam

	command {
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesGeneFamilies} > genefamilies_merged_results.tsv
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesPathways} > pathways_merged_results.tsv
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesEC} > ec_merged_results.tsv
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesKO} > ko_merged_results.tsv
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesEggnog} > eggnog_merged_results.tsv
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesGO1000} > go1000_merged_results.tsv
		/appdownload/metaphlan2/utils/merge_metaphlan_tables.py ${sep=" " filesPfam} > pfam_merged_results.tsv

	}

	output {
		File genefamiliesMerged = "genefamilies_merged_results.tsv"
		File pathwaysMerged = "pathways_merged_results.tsv"
		File ecMerged = "ec_merged_results.tsv"
		File koMerged = "ko_merged_results.tsv"
		File eggnogMerged = "eggnog_merged_results.tsv"
		File go1000Merged = "go1000_merged_results.tsv"
		File pfamMerged = "pfam_merged_results.tsv"

	}

	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 1
  		memory: "100GB"
  		preemptible: 2
  		disks: "local-disk 500 SSD"
	}


}

task strainphlanMarkers {
	
	File sampleBAM
	String sample
	
	command {
		samtools-0.1.19 view -h ${sampleBAM} -o ${sample}.sam
		/appdownload/metaphlan2/strainphlan_src/sample2markers.py --min_read_depth 5 --ifn_samples ${sample}.sam --input_type sam --output_dir . --nprocs 2 --verbose --samtools_exe /app/samtools-0.1.19
	}

	output {
		
		File sampleMarker = "${sample}.markers"
	}

	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 2
  		memory: "10GB"
  		preemptible: 2
  		disks: "local-disk 100 SSD"
	}

}

task strainphlanClades {

	Array [File] filesSampleMarkers
	File refPKL

	command <<<
		/appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --output_dir . --nprocs_main 1 --print_clades_only > clades.txt
		cat clades.txt | tr -d "'" | sed 's/,.*//' | sed 's/(//' > clades_edited.txt
	>>>

	output {
		File cladeList = "clades_edited.txt"
	}

	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 1
  		memory: "10GB"
  		preemptible: 2
  		disks: "local-disk 100 SSD"
	}

}

task strainphlanTree {
	String clade
	File refPKL
	File refAllMarkers
	File refStrainPhlanDB
	Array [File] filesSampleMarkers

	command <<<
		# get the strainphlan database
    	tar -zxvf ${refStrainPhlanDB} --strip-components=3
		
		/appdownload/metaphlan2/strainphlan_src/extract_markers.py --mpa_pkl ${refPKL} --ifn_markers ${refAllMarkers} --clade ${clade} --ofn_markers ${clade}.markers.fasta

		mkdir ${clade}

    	if [ -d strainphlan_ref/${clade} ]; then
			/appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --ifn_markers ${clade}.markers.fasta --ifn_ref_genomes strainphlan_ref/${clade}/* --output_dir ${clade} --nprocs_main 16 --clades ${clade}
		else
			/appdownload/metaphlan2/strainphlan.py --mpa_pkl ${refPKL} --ifn_samples ${sep=' ' filesSampleMarkers} --ifn_markers ${clade}.markers.fasta --output_dir ${clade} --nprocs_main 16 --clades ${clade}
		fi

		tar -czf ${clade}.tar.gz ${clade}
		
    >>>
	
	output {
		File strainphlanTar = "${clade}.tar.gz"
	}

	runtime {
		docker: "gcr.io/microbiome-xavier/metagenomicstools:070318"
		cpu: 16
  		memory: "104GB"
  		preemptible: 2
  		disks: "local-disk 100 SSD"
	}
}
