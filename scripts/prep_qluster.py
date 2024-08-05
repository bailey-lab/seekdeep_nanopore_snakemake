import subprocess
qluster_commands=snakemake.input['qluster_commands']
output_folder=snakemake.params['output_folder']
subprocess.call(f'mkdir -p {output_folder}', shell=True)
analysis_dir=snakemake.params.analysis_dir
for line in open(qluster_commands):
	split_line=line.strip().split()
	sample_index=split_line.index('-f')+1
	sample=split_line[sample_index].split('.fastq')[0]
	output_file=open(output_folder+f'/{sample}_qluster_command.sh', 'w')
	output_file.write(f'cd {analysis_dir}\n')
	output_file.write(line)
