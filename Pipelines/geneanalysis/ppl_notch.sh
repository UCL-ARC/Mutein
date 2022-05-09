# Load the necessary python libraries
module load python3/recommended
echo "~~~~~~~~~~~~~~~~~~ Call python script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
python pipeline_qsubber.py batch_geneprot.yml qsub notch NOTCH1 ""


