#!/bin/bash

################################ Slurm options #################################
### Job name
#SBATCH --job-name=Scores_P_plus_T_MEMENTO
### Max run time "hours:minutes:seconds"
#SBATCH --time=200:00:00
### Requirements nodes/servers (default: 1)
#SBATCH --nodes=1
### Requirements cpu/core/task (default: 1)
#SBATCH --ntasks-per-node=1
### Requirements memory (default: 12.5GB per cpu, requesting 50 GB)
#SBATCH --mem-per-cpu=50GB

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
# Detail Information:
scontrol show job $SLURM_JOB_ID
echo '########################################'

# load required software modules

# Change to working directory
cd /mnt/cbib/INSERM_U897/HEMA/PRS_cSVD/WMH_Single_trait_best_model_generate_score_in_UKBiobank/WMH_Residuals_UKBiobank

module load R/4.0.0


R --vanilla <<EOF


library(glmnet)
########################################################################################################################################################

#Reading file
pheno.prs <- read.table("GWAS_WMH_AMishra_31102023_fieldID_phenotypes.txt.tab_INT_WMH_and_INT_SCORE_covariates_required.txt", header=T)
head(pheno.prs)

#######################################################################################################

prs.result <- NULL
model <- lm(pheno.prs\$f.25781.2.0 ~ pheno.prs\$AGE_MRI + pheno.prs\$SEX + pheno.prs\$eTIV + pheno.prs\$PC1 + pheno.prs\$PC2 + pheno.prs\$PC3 + pheno.prs\$PC4 + pheno.prs\$PC5 + pheno.prs\$PC6 + pheno.prs\$PC7 + pheno.prs\$PC8 + pheno.prs\$PC9 + pheno.prs\$PC10)
summary(model)

model.r2 <- summary(model)\$r.squared
prs.r2 <- model.r2
prs.coef <- summary(model)\$coeff["pheno.prs\$SEX", ]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
prs.result <- rbind(prs.result, data.frame(PhenotypeName=paste0("Results_predictive_power_without_PGS_WMH_UKB"), R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
write.table(prs.result,"Results_predictive_power_of_without_PGS_WMH_in_UKB", col.names=T, quote=F, sep="\\t", row.names=F)


#######################################################################################################
prs.result <- NULL
model <- lm(pheno.prs\$f.25781.2.0 ~ pheno.prs\$AGE_MRI + pheno.prs\$SEX + pheno.prs\$eTIV + pheno.prs\$PC1 + pheno.prs\$PC2 + pheno.prs\$PC3 + pheno.prs\$PC4 + pheno.prs\$PC5 + pheno.prs\$PC6 + pheno.prs\$PC7 + pheno.prs\$PC8 + pheno.prs\$PC9 + pheno.prs\$PC10 + pheno.prs\$SCORE1_AVG)
summary(model)

model.r2 <- summary(model)\$r.squared
prs.r2 <- model.r2
prs.coef <- summary(model)\$coeff["pheno.prs\$SCORE1_AVG", ]
prs.beta <- as.numeric(prs.coef[1])
prs.se <- as.numeric(prs.coef[2])
prs.p <- as.numeric(prs.coef[4])
prs.result <- rbind(prs.result, data.frame(PhenotypeName=paste0("Results_predictive_power_PGS_WMH_UKB"), R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
write.table(prs.result,"Results_predictive_power_of_PGS_WMH_in_UKB", col.names=T, quote=F, sep="\\t", row.names=F)
#######################################################################################################

EOF

