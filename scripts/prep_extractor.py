import subprocess
extractor_commands=snakemake.input['extractor_commands']
output_folder=snakemake.params['output_folder']
subprocess.call(f'mkdir -p {output_folder}', shell=True)
analysis_dir=snakemake.params.analysis_dir
for line in open(extractor_commands):
	split_line=line.strip().split()
	sample_index=split_line.index('--sampleName')+1
	sample=split_line[sample_index][3:]
	output_file=open(output_folder+f'/{sample}_extraction_command.sh', 'w')
	output_file.write(f'cd {analysis_dir}\n')
	output_file.write(line)
