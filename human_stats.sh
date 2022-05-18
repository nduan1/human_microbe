#!/bin/sh
cd /data/EmiolaLab/duann2/human_proj/human_R_proj/human_dat/
for k in `cat chrlist.txt` 
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
Rscript human_stats.R -i $k -o /data/EmiolaLab/duann2/human_proj/human_R_proj/results/${k}.sig.tsv  -m /data/EmiolaLab/duann2/human_proj/human_R_proj/results/${k}.pdf

EOT
sbatch /data/EmiolaLab/duann2/human_proj/human_R_proj/human_dat/$k\.qsub
done

