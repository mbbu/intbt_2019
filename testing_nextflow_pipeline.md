We had issues with compatibility of our current system and the containers provided by H3ABioNet. We overcome these by creating alternative steup for running the pipleline. Here are the steps you will take to ste up and run the pipeline:
# Intbt Set up
You will be provide with a username and password during the session.

#### 1. If nextflow folder is not in your home directory, execute the script the following script

`bash installnextflow.sh`

#### 2. Activate the anaconda by loading the module

`module load anaconda3/anaconda3`

For first time use of conda, you need to initialize conda by running"

`conda init bash`

then exit shell by typing `exit` then log back in for the changes to take effect. 


#### 3. Activate Qiime conda environment

`conda activate /opt/apps/anaconda3/envs/qiime2/`

#### 4. Run the pipeline...the data is locate in a shared folder. 

`nextflow run /opt/data/Int_BT/16S-rDNA-dada2-pipeline/main.nf -profile standard -resume --reads="/opt/data/Int_BT/test-data/*_R{1,2}.fastq.gz" --trimFor 24 --trimRev 25 --reference="/opt/data/Int_BT/ref-data/silva_nr_v132_train_set.fa.gz" --species="/opt/data/Int_BT/ref-data/silva_species_assignment_v132.fa.gz" --outdir="$HOME/Int_BT/out"`
