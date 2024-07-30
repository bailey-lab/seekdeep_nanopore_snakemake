configfile: 'seekdeep_illumina_general.yaml'


sample_names=config['primer_plus_fastq_binding']+'/'+config['sample_names']
all_lines=[line.strip().split('\t') for line in open(sample_names)]
all_reps=[]
for line in all_lines[1:]:
	all_reps.extend(line[2:])
all_reps=sorted(list(set(all_reps)))
	
rule all:
	input:
#		profile=expand(config['output_folder']+'/analysis/{sample}_extraction/extractionProfile.tab.txt', sample=config['samples'])
		runnable_samples=config['output_folder']+'/non-empty_extractions.txt',

rule prep_extractor:
	input:
		extractor_commands=config['output_folder']+'/analysis/extractorCmds.txt'
	params:
		output_folder=config['output_folder']+'/extractor_shell_commands',
		analysis_dir='/home/analysis'
	output:
		all_sample_commands=expand(config['output_folder']+'/extractor_shell_commands/{sample}_extraction_command.sh', sample=all_reps)
	script:
		'scripts/prep_extractor.py'

rule run_extractor:
	input:
		data_folder=config['primer_plus_fastq_binding'],
		sif_file=config['sif_file_location'],
#		setup_done=config['output_folder']+'/finished_setup.txt',
		genome_root_folder=config['genome_binding'],
		actual_shell_script=config['output_folder']+'/extractor_shell_commands/{sample}_extraction_command.sh'	
	params:
		output_dir=config['output_folder'],
		softlink_fastq_binding=config['softlink_fastq_binding'],
		singularity_shell_script='/seekdeep_output/extractor_shell_commands/{sample}_extraction_command.sh'
	output:
		profile=config['output_folder']+'/analysis/{sample}_extraction/extractionProfile.tab.txt',
		folder=directory(config['output_folder']+'/analysis/{sample}_extraction')
	resources:
		time_min=config['max_run_time_min'],
		mem_mb=config['max_memory_mb'],
		nodes=config['cpus_to_use']
	shell:
		'''
		singularity exec \
			-B {input.data_folder}:/input_data \
			-B {params.output_dir}:/seekdeep_output \
			-B {input.genome_root_folder}:/genome_info \
			-B {params.output_dir}/analysis:/home/analysis \
			{params.softlink_fastq_binding} \
			{input.sif_file} bash {params.singularity_shell_script}
		'''

rule analyze_extractor:
	input:
		profiles=expand(config['output_folder']+'/analysis/{sample}_extraction/extractionProfile.tab.txt', sample=all_reps),
		folders=expand(config['output_folder']+'/analysis/{sample}_extraction', sample=all_reps)
	output:
		runnable_samples=config['output_folder']+'/non-empty_extractions.txt'
	script:
		'scripts/analyze_extractor.py'
