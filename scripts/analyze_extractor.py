import os
folders=snakemake.input['folders']
runnable_samples=open(snakemake.output['runnable_samples'], 'w')

# def find_fastqs(folder):
# 	good_files=[]
# 	for root, dirs, files in os.walk(folder):
# 		for file_name in files:
# 			if file_name.endswith('.fastq.gz') and file_name[:-9]!='output':
# 				good_files.append(root+'/'+file_name)
# 	return good_files

def find_fastqs(folder):
	good_files=[folder + '/' + f for f in os.listdir(folder) 
					if f.endswith('.fastq.gz') 
					and 'undetermined' not in f]
	return good_files

for folder in folders:
	read_fastq='something'
	good_files=find_fastqs(folder)
	if len(good_files)>0:
		for good_file in good_files:
			if os.path.getsize(good_file)>0:
				runnable_samples.write(good_file.split('/')[-1].split('.fastq')[0]+'\n')
