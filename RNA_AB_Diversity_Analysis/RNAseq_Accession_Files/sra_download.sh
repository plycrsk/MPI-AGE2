#!/bin/bash
#SBATCH -p hooli
#SBATCH -o /beegfs/group_dv/home/MMihaljevic/sam/test.singularity.out
#SBATCH --mem=15gb
#SBATCH --cpus-per-task 24
# set name of job
#SBATCH --job-name=exec
#SBATCH --error /beegfs/group_dv/home/MMihaljevic/sam/test.singularity.err


FILENAME="accession_RNA.txt"
LINES=$(cat $FILENAME)

echo $LINES
echo $FILENAME

cd ~/sam
for LINE in $LINES
	do
		if [ -e data_RNA/${LINE}_pass_1.fastq.gz.done ]; then
			echo data_RNA/${LINE}_pass_1.fastq.gz.done finished;
		else
	  	echo ${LINE}	
	  	sbatch --cpus-per-task=18 --mem=15gb --time=5-24 -p hooli -o /beegfs/group_dv/home/MMihaljevic/sam/logs_RNA/${LINE}.%j.out <<EOF
#!/bin/bash
singularity exec ~/bioinformatics_software.v3.0.0.sif /bin/bash << SHI
#!/bin/bash
source ~/.bashrc
module load sratoolkit rlang
echo ${LINE} "processing"
fastq-dump --outdir data_RNA --gzip --skip-technical --split-3 --readids --read-filter pass --dumpbase --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' --clip ${LINE}
SHI
touch data_RNA/data/${LINE}_pass_1.fastq.gz.done
EOF
fi
done
exit
