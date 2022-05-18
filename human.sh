#!/bin/sh
cd /data/EmiolaLab/duann2/human_proj/human_R_proj/human_dat/getab/
for k in chr* 
do
cat <<EOT >> $k\.qsub
#!/bin/sh
#SBATCH --partition=norm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --gres=lscratch:500
#SBATCH --time=8-02:00:00
module load R
for i in \`cat $k \`
do
Rscript human.R -i /data/EmiolaLab/duann2/human_proj/human_R_proj/dataset/\$i\.tsv -o /data/EmiolaLab/duann2/human_proj/human_R_proj/results/\$i\.output.tsv -m /data/EmiolaLab/duann2/human_proj/human_R_proj/results/\$i\.pdf
done

EOT
sbatch /data/EmiolaLab/duann2/human_proj/human_R_proj/human_dat/getab/$k\.qsub
done

