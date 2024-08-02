configfile: 'seekdeep_nanopore_general_experimental.yaml'


out_folder=config['output_folder']
rule all:
	input:
		setup_done=out_folder+'/finished_setup.txt',
		out_snakefile=out_folder+'/snakemake_params/setup_run.smk',
		out_config_file=out_folder+'/snakemake_params/seekdeep_nanopore_general.yaml'

rule copy_files:
	'''
	copies snakemake script and config files to output folder for reproducibility.
	'''
	input:
		setup_file='setup_run.smk',
		extractor_file='run_extractor.smk',
		finish_file='finish_process.smk',
		all_steps = 'run_pipeline.sh',
		scripts='scripts',
		config_file='seekdeep_nanopore_general_experimental.yaml'
	output:
		setup_file=out_folder+'/snakemake_params/setup_run.smk',
		extractor_file=out_folder+'/snakemake_params/run_extractor.smk',
		finish_file=out_folder+'/snakemake_params/finish_process.smk',
		all_steps=out_folder+'/snakemake_params/run_pipeline.sh',
		scripts=directory(out_folder+'/snakemake_params/scripts'),
		config_file=out_folder+'/snakemake_params/seekdeep_nanopore_general.yaml'
	shell:
		'''
		cp {input.setup_file} {output.setup_file}
		cp {input.extractor_file} {output.extractor_file}
		cp {input.finish_file} {output.finish_file}
		cp {input.all_steps} {output.all_steps}
		cp -r {input.scripts} {output.scripts}
		cp {input.config_file} {output.config_file}
		'''

rule genTargetInfoFromGenomes:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		genome_root_folder=config['genome_binding'],
		sif_file=config['sif_file_location']
	params:
		output_dir=out_folder,
		primer_file=config['primer_file'],
		gff_subfolder=config['gff_subfolder'],
		genome_subfolder=config['genome_subfolder'],
		fake_insert_size = config['fake_insert_size'],
		extra_args=config['extra_gen_target_info_cmds']
	output:
		primer_info=out_folder+'/extractedRefSeqs/locationsByGenome/Pf3D7_infos.tab.txt'
	threads: config['cpus_to_use']
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec \
			-B {input.genome_root_folder}:/genome_info \
			-B {input.data_folder}:/input_data \
			-B {params.output_dir}:/seekdeep_output \
			{input.sif_file} SeekDeep genTargetInfoFromGenomes \
				--primers /input_data/{params.primer_file} \
				--genomeDir /genome_info/{params.genome_subfolder} \
				--gffDir /genome_info/{params.gff_subfolder} \
				--pairedEndLength {params.fake_insert_size} \
				{params.extra_args} \
				--dout /seekdeep_output/extractedRefSeqs \
				--overWriteDir \
				--numThreads {threads} \
				--shortnames
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
		setup_done=config['output_folder']+'/finished_setup.txt',
	threads: config['cpus_to_use']
	shell:
		'''
		singularity exec \
			-B {input.data_folder}:/input_data \
			-B {params.output_dir}:/seekdeep_output \
			{params.softlink_fastq_binding} \
			{input.sif_file} SeekDeep setupTarAmpAnalysis \
				--outDir /seekdeep_output/analysis \
				--technology nanopore \
				--uniqueKmersPerTarget {params.for_seekdeep}/uniqueKmers.tab.txt.gz \
				--inputDir /input_data/{params.fastq_folder} \
				--idFile /input_data/{params.primer_file} \
				--lenCutOffs {params.for_seekdeep}/lenCutOffs.txt \
				--doNotGuessRecFlags {params.extra_extractor_cmds} \
				{params.extra_kluster_cmds} \
				{params.extra_process_cluster_cmds} \
				--previousPopSeqsDir {params.for_seekdeep}/refSeqs/ \
				--numThreads {threads}
			touch {output.setup_done}
		'''
