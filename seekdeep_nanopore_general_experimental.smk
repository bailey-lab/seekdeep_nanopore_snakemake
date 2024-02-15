configfile: 'seekdeep_nanopore_general_experimental.yaml'
rule all:
	input:
		analysis_done=config['output_folder']+'/finished_analysis.txt',
		out_snakefile=config['output_folder']+'/seekdeep_nanopore_general_experimental.smk',
		out_config_file=config['output_folder']+'/seekdeep_nanopore_general_experimental.yaml'

rule copy_files:
	'''
	copies snakemake script and config files to output folder for reproducibility.
	'''
	input:
		in_snakefile='seekdeep_nanopore_general_experimental.smk',
		in_config_file='seekdeep_nanopore_general_experimental.yaml'
	output:
		out_snakefile=config['output_folder']+'/seekdeep_nanopore_general_experimental.smk',
		out_config_file=config['output_folder']+'/seekdeep_nanopore_general_experimental.yaml'
	shell:
		'''
		cp {input.in_snakefile} {output.out_snakefile}
		cp {input.in_config_file} {output.out_config_file}
		'''

rule genTargetInfoFromGenomes:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		genome_root_folder=config['genome_binding'],
		sif_file=config['sif_file_location']
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		gff_subfolder=config['gff_subfolder'],
		genome_subfolder=config['genome_subfolder'],
		extra_args=config['extra_gen_target_info_cmds']
	output:
		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.genome_root_folder}:/genome_info \
		-B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output {input.sif_file} \
		SeekDeep genTargetInfoFromGenomes --primers /input_data/{params.primer_file} \
		--longRangeAmplicon --useBlast --genomeDir /genome_info/{params.genome_subfolder} \
		--gffDir /genome_info/{params.gff_subfolder} {params.extra_args} \
		--dout /seekdeep_output/extractedRefSeqs --overWriteDir --numThreads \
		{threads} --shortNames
		'''

rule setupTarAmpAnalysis:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		primer_info=config['output_folder']+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	params:
		output_dir=config['output_folder'],
		primer_file=config['primer_file'],
		fastq_folder=config['fastq_subfolder'],
		for_seekdeep='/seekdeep_output/extractedRefSeqs/forSeekDeep',
		softlink_fastq_binding=config['softlink_fastq_binding'],
		extra_extractor_cmds=config['extra_extractor_cmds'],
		extra_kluster_cmds=config['extra_kluster_cmds'],
		extra_process_cluster_cmds=config['extra_process_cluster_cmds']
	output:
		setup_done=config['output_folder']+'/finished_setup.txt'
	threads: config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		{params.softlink_fastq_binding} {input.sif_file} \
		SeekDeep setupTarAmpAnalysis  --outDir \
		/seekdeep_output/analysis --technology nanopore \
		--uniqueKmersPerTarget \
		{params.for_seekdeep}/uniqueKmers.tab.txt.gz \
		--inputDir /input_data/{params.fastq_folder} \
		--idFile /input_data/{params.primer_file} --lenCutOffs \
		{params.for_seekdeep}/lenCutOffs.txt \
		--doNotGuessRecFlags {params.extra_extractor_cmds} \
		{params.extra_kluster_cmds} {params.extra_process_cluster_cmds} \
		--previousPopSeqsDir {params.for_seekdeep}/refSeqs/ \
		--numThreads {threads}
		touch {output.setup_done}
		'''

rule setup_analysis:
	'''
	A kludgy shell script solution to get snakemake to change directory to the
	analysis directory, where runAnalysis.sh uses relative paths
	'''
	input:
		setup_done=config['output_folder']+'/finished_setup.txt',
	threads: config['cpus_to_use']
	output:
		final_script=config['output_folder']+'/analysis_script.sh',
	shell:
		'''
		echo cd $'/seekdeep_output/analysis\n./runAnalysis.sh {threads}' >{output.final_script}
		'''

rule runAnalysis:
	input:
		final_script=config['output_folder']+'/analysis_script.sh',
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
		genome_root_folder=config['genome_binding']
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding']
	output:
		analysis_done=config['output_folder']+'/finished_analysis.txt'
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec -B {input.data_folder}:/input_data \
		-B {params.output_dir}:/seekdeep_output \
		-B {input.genome_root_folder}:/genome_info \
		{params.softlink_fastq_binding} \
		{input.sif_file} bash /seekdeep_output/analysis_script.sh
		touch {output.analysis_done}
		'''
