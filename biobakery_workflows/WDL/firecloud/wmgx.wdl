workflow workflowWMGX {

	File Input_File_Sample_List
	String Fastq_Pair1_Extension
	String Fastq_Pair2_Extension
	Map[String, String] SamplesPaths = read_map(Input_File_Sample_List)
	scatter (pair in SamplesPaths){

		String SampleDir = pair.right
		String sample = pair.left

		File F1 = SampleDir + sample + Fastq_Pair1_Extension
		File F2 = SampleDir + sample + Fastq_Pair2_Extension
		
		call qcQualityHuman {
			input: 
			sample=sample, 
			file1=F1, 
			file2=F2
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

