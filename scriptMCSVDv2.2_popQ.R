

# MCSVD_script.R by Dr. Aniket Mishra and Pr. Stephanie Debette
# To create MRI defined CSVD extreme phenotypes
# For queries or bugs please contact aniket.mishra@u-bordeaux.fr

# Usage= Rscript MCSVD_script.R MCSVD.input
# Make sure input file contains exactly 20 columns in following order 1)Sample ID 2)Age 3)Sex 4)log(WMHV+1) 5)Mask Volume or Total Intracranial Volume 6)Presence of any brain infarct (coded 0/1) 7)Presence of lacunar brain infarct (coded 0/1)
# 8)WMH, 9)DBP, 10)SBP, 11)HT_status, 12)Fast_Glucose, 13)Diabetes, 14)HDL, 15)LDL, 16)TG, 17)BMI, 18)CVD, 19)AntiHT_drug, 20)lipid_low_drug
# For descrete variables 0 = absent and 1 = present

# Run already avaialble "Hmisc" library or run the "cutt" function
#install.packages("Hmisc")
#library(Hmisc)
cat("\n")
cat("\t\tThis R script categorises individulas into the extremes of MRI defined CSVD (MCSVD) and reports sample characteristics\n")
cat("\t\tIt was developed by Dr. Aniket Mishra and Pr. Stephanie Debette\n")
cat("\t\tIf you encounter any issues please contact Aniket Mishra at aniket.mishra@u-bordeaux.fr\n")
cat("\n")

# if you encounter difficulties installing the "Hmisc" library run the cutt function available from the Hmisc package as follows
cutt<- function (x, cuts, m = 150, g, levels.mean = FALSE, digits, minmax = TRUE, 
          oneval = TRUE, onlycuts = FALSE) 
{
  method <- 1
  x.unique <- sort(unique(c(x[!is.na(x)], if (!missing(cuts)) cuts)))
  min.dif <- min(diff(x.unique))/2
  min.dif.factor <- 1
  if (missing(digits)) 
    digits <- if (levels.mean) 
      5
  else 3
  oldopt <- options(digits = digits)
  on.exit(options(oldopt))
  xlab <- attr(x, "label")
  if (missing(cuts)) {
    nnm <- sum(!is.na(x))
    if (missing(g)) 
      g <- max(1, floor(nnm/m))
    if (g < 1) 
      stop("g must be >=1, m must be positive")
    options(digits = 15)
    n <- table(x)
    xx <- as.double(names(n))
    options(digits = digits)
    cum <- cumsum(n)
    m <- length(xx)
    y <- as.integer(ifelse(is.na(x), NA, 1))
    labs <- character(g)
    cuts <- approx(cum, xx, xout = (1:g) * nnm/g, method = "constant", 
                   rule = 2, f = 1)$y
    cuts[length(cuts)] <- max(xx)
    lower <- xx[1]
    upper <- 1e+45
    up <- low <- double(g)
    i <- 0
    for (j in 1:g) {
      cj <- if (method == 1 || j == 1) 
        cuts[j]
      else {
        if (i == 0) 
          stop("program logic error")
        s <- if (is.na(lower)) 
          FALSE
        else xx >= lower
        cum.used <- if (all(s)) 
          0
        else max(cum[!s])
        if (j == m) 
          max(xx)
        else if (sum(s) < 2) 
          max(xx)
        else approx(cum[s] - cum.used, xx[s], xout = (nnm - 
                                                        cum.used)/(g - j + 1), method = "constant", 
                    rule = 2, f = 1)$y
      }
      if (cj == upper) 
        next
      i <- i + 1
      upper <- cj
      y[x >= (lower - min.dif.factor * min.dif)] <- i
      low[i] <- lower
      lower <- if (j == g) 
        upper
      else min(xx[xx > upper])
      if (is.na(lower)) 
        lower <- upper
      up[i] <- lower
    }
    low <- low[1:i]
    up <- up[1:i]
    variation <- logical(i)
    for (ii in 1:i) {
      r <- range(x[y == ii], na.rm = TRUE)
      variation[ii] <- diff(r) > 0
    }
    if (onlycuts) 
      return(unique(c(low, max(xx))))
    flow <- format(low)
    fup <- format(up)
    bb <- c(rep(")", i - 1), "]")
    labs <- ifelse(low == up | (oneval & !variation), flow, 
                   paste("[", flow, ",", fup, bb, sep = ""))
    ss <- y == 0 & !is.na(y)
    if (any(ss)) 
      stop(paste("categorization error in cut2.  Values of x not appearing in any interval:\n", 
                 paste(format(x[ss], digits = 12), collapse = " "), 
                 "\nLower endpoints:", paste(format(low, digits = 12), 
                                             collapse = " "), "\nUpper endpoints:", paste(format(up, 
                                                                                                 digits = 12), collapse = " ")))
    y <- structure(y, class = "factor", levels = labs)
  }
  else {
    if (minmax) {
      r <- range(x, na.rm = TRUE)
      if (r[1] < cuts[1]) 
        cuts <- c(r[1], cuts)
      if (r[2] > max(cuts)) 
        cuts <- c(cuts, r[2])
    }
    l <- length(cuts)
    k2 <- cuts - min.dif
    k2[l] <- cuts[l]
    y <- cut(x, k2)
    if (!levels.mean) {
      brack <- rep(")", l - 1)
      brack[l - 1] <- "]"
      fmt <- format(cuts)
      labs <- paste("[", fmt[1:(l - 1)], ",", fmt[2:l], 
                    brack, sep = "")
      if (oneval) {
        nu <- table(cut(x.unique, k2))
        if (length(nu) != length(levels(y))) 
          stop("program logic error")
        levels(y) <- ifelse(nu == 1, c(fmt[1:(l - 2)], 
                                       fmt[l]), labs)
      }
      else levels(y) <- labs
    }
  }
  if (levels.mean) {
    means <- tapply(x, y, function(w) mean(w, na.rm = TRUE))
    levels(y) <- format(means)
  }
  attr(y, "class") <- "factor"
  if (length(xlab)) 
    label(y) <- xlab
  y
}



# Arguments
args                  <- commandArgs(TRUE)

# Read input file
#DATA                  <- read.delim(args[1])
# DATA<- read.delim("MCSVD.input.example")
DATA <- read.table("GWAS_WMH_AMishra_31102023_fieldID_phenotypes.txt.tab_INT_WMH_and_INT_SCORE_covariates_required_20columns_data2.txt", header=T)
cat(paste("\tNumber of rows in provided input file: ", nrow(DATA), "\n",sep=""))
# Extract first seven columns
cat("\tExtracting first seven columns to crete S1trans file")
DATA2                 <- DATA[,c(1:7)] 

# Input file check
# Remove missing values
INPUT                 <- na.omit(DATA2)
cat(paste("\tNumber of rows in S1trans file after NA filter: ", nrow(INPUT), "\n", sep=""))
if(nrow(INPUT) < 100){
    print("ATTENTION: NA filter removed large sample")
}

#check whether the number of columns are equal to 7
if(ncol(INPUT) != 7){
    stop("ERROR 1: Number of columns is not equal to 7")
}

names(INPUT)          <- c("ID", "Age", "Sex", "l1WMHV", "MaskVolume", "AnyBI", "LBI")

#check column 6 and 7 are coded 0/1
if(min(INPUT$AnyBI) < 0){
    stop("ERROR 2: The minimum Any Brain Infarct (column 6 in Input file) value is less than 0")
}
if(max(INPUT$AnyBI) > 1){
    stop("ERROR 3: The maximum Any Brain Infarct value (column 6 in Input file) is more than 1")
}

if(min(INPUT$LBI) < 0){
    stop("ERROR 4:The minimum Lacunar Brain Infarct value (column 7 in Input file) is less than 0")
}
if(max(INPUT$LBI) > 1){
    stop("ERROR 5:The maximum Lacunar Brain Infarct value (column 7 in Input file) is more than 1")
}

# Define N_Target
# N_TARGET is 1/6th of the total sample size, the number is rounded
N_TARGET              <- round(length(INPUT[,1]) / 6)
cat(paste("\tTarget N in the top and bottom extremes (1/6th of total N): ", N_TARGET, "\n", sep = ""))

# Calculate residuals of l1WMHV adjusted for Age, Sex and MaskVolume
WMHV_regression       <- lm(l1WMHV ~ Age + Sex + MaskVolume, data = INPUT)
INPUT$WMHV_residuals  <- resid(WMHV_regression)

# Based on the WMHV residuals categorise individuals into quantile groups
INPUT$WMHV_quantiles  <- as.numeric(cutt(INPUT$WMHV_residuals,g = 4))

# Here we divide individuals in four quantiles (quartiles here) and define MCSVD_1 and MCSVD_2 status as follows
TOP_Q                 <- subset(INPUT, INPUT$WMHV_quantiles == 4)
TOP_Q_LENGTH          <- nrow(TOP_Q)
cat(paste("\tTop quantile length: ", TOP_Q_LENGTH, "\n", sep=""))

TOP_Q_ORDERED         <- TOP_Q[order(TOP_Q$LBI, TOP_Q$WMHV_residuals, decreasing = T),]
TOP_Q_ORDERED$MCSVD_1 <- c(rep(1, N_TARGET), rep(NA, TOP_Q_LENGTH - N_TARGET))
TOP_Q_ORDERED$MCSVD_2 <- ifelse(TOP_Q_ORDERED$LBI == 1, 1, "NA")

MID_Q1                <- subset(INPUT, INPUT$WMHV_quantiles == 2)
MID_Q2                <- subset(INPUT, INPUT$WMHV_quantiles == 3)
MID_Q                 <- rbind(MID_Q1, MID_Q2)
MID_Q_LENGTH          <- nrow(MID_Q)
cat(paste("\tIndviduals not in extreme quantiles length: ", MID_Q_LENGTH, "\n", sep=""))

MID_Q$MCSVD_1         <- rep("NA", MID_Q_LENGTH)
MID_Q$MCSVD_2         <- rep("NA", MID_Q_LENGTH)

BOT_Q                 <- subset(INPUT, INPUT$WMHV_quantiles == 1)
BOT_Q_LENGTH          <- nrow(BOT_Q)
cat(paste("\tBottom quantile length: ", BOT_Q_LENGTH, "\n", sep=""))

BOT_Q_ORDERED         <- BOT_Q[order(BOT_Q$AnyBI, BOT_Q$WMHV_residuals),]
BOT_Q_ORDERED$MCSVD_1 <- c(rep(0, N_TARGET), rep(NA, BOT_Q_LENGTH - N_TARGET))
N_MCSVD2              <- length(which(TOP_Q_ORDERED$MCSVD_2 == 1))
BOT_Q_ORDERED$MCSVD_2 <- c(rep(0, N_TARGET), rep(NA, BOT_Q_LENGTH - N_TARGET))


#Double check if any control individual has a brain infarct
BOT_Q_FLAG_1_LENGTH   <- length(which(BOT_Q_ORDERED$MCSVD_1 == 0 && BOT_Q_ORDERED$AnyBI == 1))

if(BOT_Q_FLAG_1_LENGTH > 0){
    print("WARNING 1: Following MCSVD_1 controls have BI")
    BOT_Q_ORDERED[which(BOT_Q_ORDERED$MCSVD_1 == 0 && BOT_Q_ORDERED$AnyBI == 1),]
}

BOT_Q_FLAG_2_LENGTH   <- length(which(BOT_Q_ORDERED$MCSVD_2 == 0 && BOT_Q_ORDERED$AnyBI == 1))

if(BOT_Q_FLAG_1_LENGTH > 0){
    print("WARNING 2: Following MCSVD_2 controls have BI:")
    BOT_Q_ORDERED[which(BOT_Q_ORDERED$MCSVD_2 == 0 && BOT_Q_ORDERED$AnyBI == 1),]
}

# rbind Q files
FINAL                 <- rbind(TOP_Q_ORDERED, MID_Q, BOT_Q_ORDERED)

# Write the output file
write.table(FINAL, "MCSVD.out", quote=F, row.names=F, col.names=T, sep= "\t")

# Make sure file FINAL2 object has 24 columns from column 12th the order should be : WMH(12), DBP(13), SBP(14), 
# HT_status(15), Fast_Glucose(16), Diabetes(17), HDL(18), LDL(19), TG(20), BMI(21), CVD(22), AntiHT_drug(23), lipid_low_drug(24) 

#DATA3                 <- DATA[,-c(1:7)] # old bug
DATA3                 <- DATA[,-c(2:7)]

cat(paste("\tNumber of rows in S1out file: ", nrow(FINAL), "\n", sep=""))
cat(paste("\tNumber of rows in S2trans file: ", nrow(DATA3), "\n", sep=""))
cat(paste("\tNumber of columns in S1out file: ", ncol(FINAL), "\n", sep=""))
cat(paste("\tNumber of columns in S2trans file: ", ncol(DATA3), "\n", sep=""))

#FINAL2<- cbind(FINAL, DATA3) # old bug

FINAL2<- merge(FINAL, DATA3, by="ID", all.y=T)

# Write the output file
write.table(FINAL2, "MCSVD_withvariables.out", quote=F, row.names=F, col.names=T, sep= "\t")

#check whether the number of columns are equal to 24
if(ncol(FINAL2) != 24){
    stop("ERROR 6: For sample characteristics number of columns is not equal to 24, input file should have 20 columns")
}

cat("\tNumber of columns in S2input file: 24")

names(FINAL2)          <- c("ID", "Age", "Sex", "l1WMHV", "MaskVolume", "AnyBI", "LBI", "WMHV_residuals", "WMHV_quantiles", 
"MCSVD_1", "MCSVD_2", "WMH", "DBP", "SBP", "HT_status", "Fast_Glucose", "Diabetes", "HDL", "LDL", "TG", "BMI", "CVD", 
"AntiHT_drug", "lipid_low_drug")


# Report sample characteristics for MCSVD1
cat("\n")
cat("\tSample Caracteristics for MCSVD_1\n")

MCSVD_1_cases 		<- FINAL2[which(FINAL2$MCSVD_1 == 1),]
MCSVD_1_controls 	<- FINAL2[which(FINAL2$MCSVD_1 == 0),]
MCSVD1_All		<- rbind(MCSVD_1_cases, MCSVD_1_controls)

# Descriptive table would look like
# Variable					ext-MCSVD1		min-MCSVD1		p-value
# N
# Age_mean(s.d)
# l1WMH_mean(s.d)
# WMH_mean(s.d)
# WMHresiduals_mean(s.d)
# TIV_mean(s.d)
# DBP_mean(s.d)
# SBP_mean(s.d)
# Fast_Glucose_mean(s.d)
# HDL_mean(s.d)
# LDL_mean(s.d)
# TG_mean(s.d)
# BMI_mean(s.d)
# Females_N(percentage)
# LBI_status_N(percentage)
# HT_status_N(percentage)
# Diabetes_status_N(percentage)
# CVD_status_N(percentage)
# AntiHT_drug_status_N(percentage)
# Lipid_low_drug_status_N(percentage)

DISCREPTIVE_MCSVD1	<-data.frame(matrix(nrow=20, ncol=4))
DISCREPTIVE_MCSVD1[1,]	<- c("N", length(MCSVD_1_cases[,1]), length(MCSVD_1_controls[,1]), "NA")
# Quantitative variables
DISCREPTIVE_MCSVD1[2,]	<- c("Age_mean(s.d)", paste(mean(MCSVD_1_cases$Age, na.rm=T),"(", sd(MCSVD_1_cases$Age, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$Age, na.rm=T),"(", sd(MCSVD_1_controls$Age, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$Age, MCSVD_1_controls$Age, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[3,]	<- c("l1WMH_mean(s.d)", paste(mean(MCSVD_1_cases$l1WMHV, na.rm=T),"(", sd(MCSVD_1_cases$l1WMHV, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$l1WMHV, na.rm=T),"(", sd(MCSVD_1_controls$l1WMHV, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$l1WMHV, MCSVD_1_controls$l1WMHV, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[4,]	<- c("WMH_mean(s.d)", paste(mean(MCSVD_1_cases$WMH, na.rm=T),"(", sd(MCSVD_1_cases$WMH, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$WMH, na.rm=T),"(", sd(MCSVD_1_controls$WMH, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$WMH, MCSVD_1_controls$WMH, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[5,]	<- c("WMHresiduals_mean(s.d)", paste(mean(MCSVD_1_cases$WMHV_residuals, na.rm=T),"(", sd(MCSVD_1_cases$WMHV_residuals, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$WMHV_residuals, na.rm=T),"(", sd(MCSVD_1_controls$WMHV_residuals, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$WMHV_residuals, MCSVD_1_controls$WMHV_residuals, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[6,]	<- c("TIV_mean(s.d)", paste(mean(MCSVD_1_cases$MaskVolume, na.rm=T),"(", sd(MCSVD_1_cases$MaskVolume, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$MaskVolume, na.rm=T),"(", sd(MCSVD_1_controls$MaskVolume, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$MaskVolume, MCSVD_1_controls$MaskVolume, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[7,]	<- c("DBP_mean(s.d)", paste(mean(MCSVD_1_cases$DBP, na.rm=T),"(", sd(MCSVD_1_cases$DBP, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$DBP, na.rm=T),"(", sd(MCSVD_1_controls$DBP, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$DBP, MCSVD_1_controls$DBP, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[8,]	<- c("SBP_mean(s.d)", paste(mean(MCSVD_1_cases$SBP, na.rm=T),"(", sd(MCSVD_1_cases$SBP, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$SBP, na.rm=T),"(", sd(MCSVD_1_controls$SBP, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$SBP, MCSVD_1_controls$SBP, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[9,]	<- c("Fast_Glucose_mean(s.d)", paste(mean(MCSVD_1_cases$Fast_Glucose, na.rm=T),"(", sd(MCSVD_1_cases$Fast_Glucose, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$Fast_Glucose, na.rm=T),"(", sd(MCSVD_1_controls$Fast_Glucose, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$Fast_Glucose, MCSVD_1_controls$Fast_Glucose, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[10,]	<- c("HDL_mean(s.d)", paste(mean(MCSVD_1_cases$HDL, na.rm=T),"(", sd(MCSVD_1_cases$HDL, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$HDL, na.rm=T),"(", sd(MCSVD_1_controls$HDL, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$HDL, MCSVD_1_controls$HDL, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[11,]	<- c("LDL_mean(s.d)", paste(mean(MCSVD_1_cases$LDL, na.rm=T),"(", sd(MCSVD_1_cases$LDL, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$LDL, na.rm=T),"(", sd(MCSVD_1_controls$LDL, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$LDL, MCSVD_1_controls$LDL, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[12,]	<- c("TG_mean(s.d)", paste(mean(MCSVD_1_cases$TG, na.rm=T),"(", sd(MCSVD_1_cases$TG, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$TG, na.rm=T),"(", sd(MCSVD_1_controls$TG, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$TG, MCSVD_1_controls$TG, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD1[13,]	<- c("BMI_mean(s.d)", paste(mean(MCSVD_1_cases$BMI, na.rm=T),"(", sd(MCSVD_1_cases$BMI, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_1_controls$BMI, na.rm=T),"(", sd(MCSVD_1_controls$BMI, na.rm=T), ")", sep=""), 
	t.test(MCSVD_1_cases$BMI, MCSVD_1_controls$BMI, var.equal = F, paired=F)$p.value)
# Now categorical variables
DISCREPTIVE_MCSVD1[14,]	<- c("Females_N(percentage)", paste(length(which(MCSVD_1_cases$Sex == 1)),"(", length(which(MCSVD_1_cases$Sex == 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$Sex == 1)),"(", length(which(MCSVD_1_controls$Sex == 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD1_All$Sex, MCSVD1_All$MCSVD_1))$p.value)
DISCREPTIVE_MCSVD1[15,]	<- c("LBI_status_N(percentage)", paste(length(which(MCSVD_1_cases$LBI== 1)),"(", length(which(MCSVD_1_cases$LBI== 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$LBI== 1)),"(", length(which(MCSVD_1_controls$LBI== 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	"NA") # It doesnt make sense to test association LBI, since we already know it will be 0 in min-CSVD
DISCREPTIVE_MCSVD1[16,]	<- c("HT_status_N(percentage)", paste(length(which(MCSVD_1_cases$HT_status== 1)),"(", length(which(MCSVD_1_cases$HT_status== 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$HT_status== 1)),"(", length(which(MCSVD_1_controls$HT_status== 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD1_All$HT_status, MCSVD1_All$MCSVD_1))$p.value)
DISCREPTIVE_MCSVD1[17,]	<- c("Diabetes_status_N(percentage)", paste(length(which(MCSVD_1_cases$Diabetes== 1)),"(", length(which(MCSVD_1_cases$Diabetes== 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$Diabetes== 1)),"(", length(which(MCSVD_1_controls$Diabetes== 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD1_All$Diabetes, MCSVD1_All$MCSVD_1))$p.value)
DISCREPTIVE_MCSVD1[18,]	<- c("CVD_status_N(percentage)", paste(length(which(MCSVD_1_cases$CVD== 1)),"(", length(which(MCSVD_1_cases$CVD== 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$CVD== 1)),"(", length(which(MCSVD_1_controls$CVD== 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD1_All$CVD, MCSVD1_All$MCSVD_1))$p.value)
DISCREPTIVE_MCSVD1[19,]	<- c("AntiHT_drug_status_N(percentage)", paste(length(which(MCSVD_1_cases$AntiHT_drug== 1)),"(", length(which(MCSVD_1_cases$AntiHT_drug== 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$AntiHT_drug== 1)),"(", length(which(MCSVD_1_controls$AntiHT_drug== 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD1_All$AntiHT_drug, MCSVD1_All$MCSVD_1))$p.value)
DISCREPTIVE_MCSVD1[20,]	<- c("Lipid_low_drug_status_N(percentage)", paste(length(which(MCSVD_1_cases$lipid_low_drug== 1)),"(", length(which(MCSVD_1_cases$lipid_low_drug== 1))/length(MCSVD_1_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_1_controls$lipid_low_drug== 1)),"(", length(which(MCSVD_1_controls$lipid_low_drug== 1))/length(MCSVD_1_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD1_All$lipid_low_drug, MCSVD1_All$MCSVD_1))$p.value)

names(DISCREPTIVE_MCSVD1) <-c("Variables","Extensive_MCSVD", "Minimal_MCSVD", "pvalue") 
write.table(DISCREPTIVE_MCSVD1, "MCSVD1_Descriptives.out", quote=F, row.names=F, col.names=T, sep= "\t")



# Report sample characteristics for MCSVD2
cat("\n")
cat("\tSample Caracteristics for MCSVD_2\n")

MCSVD_2_cases 		<- FINAL2[which(FINAL2$MCSVD_2 == 1),]
MCSVD_2_controls 	<- FINAL2[which(FINAL2$MCSVD_2 == 0),]
MCSVD2_All		<- rbind(MCSVD_2_cases, MCSVD_2_controls)

DISCREPTIVE_MCSVD2	<- data.frame(matrix(nrow=20, ncol=4))
DISCREPTIVE_MCSVD2[1,]	<- c("N", length(MCSVD_2_cases[,1]), length(MCSVD_2_controls[,1]), "NA")
# Quantitative variables
DISCREPTIVE_MCSVD2[2,]	<- c("Age_mean(s.d)", paste(mean(MCSVD_2_cases$Age, na.rm=T),"(", sd(MCSVD_2_cases$Age, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$Age, na.rm=T),"(", sd(MCSVD_2_controls$Age, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$Age, MCSVD_2_controls$Age, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[3,]	<- c("l1WMH_mean(s.d)", paste(mean(MCSVD_2_cases$l1WMHV, na.rm=T),"(", sd(MCSVD_2_cases$l1WMHV, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$l1WMHV, na.rm=T),"(", sd(MCSVD_2_controls$l1WMHV, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$l1WMHV, MCSVD_2_controls$l1WMHV, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[4,]	<- c("WMH_mean(s.d)", paste(mean(MCSVD_2_cases$WMH, na.rm=T),"(", sd(MCSVD_2_cases$WMH, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$WMH, na.rm=T),"(", sd(MCSVD_2_controls$WMH, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$WMH, MCSVD_2_controls$WMH, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[5,]	<- c("WMHresiduals_mean(s.d)", paste(mean(MCSVD_2_cases$WMHV_residuals, na.rm=T),"(", sd(MCSVD_2_cases$WMHV_residuals, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$WMHV_residuals, na.rm=T),"(", sd(MCSVD_2_controls$WMHV_residuals, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$WMHV_residuals, MCSVD_2_controls$WMHV_residuals, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[6,]	<- c("TIV_mean(s.d)", paste(mean(MCSVD_2_cases$MaskVolume, na.rm=T),"(", sd(MCSVD_2_cases$MaskVolume, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$MaskVolume, na.rm=T),"(", sd(MCSVD_2_controls$MaskVolume, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$MaskVolume, MCSVD_2_controls$MaskVolume, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[7,]	<- c("DBP_mean(s.d)", paste(mean(MCSVD_2_cases$DBP, na.rm=T),"(", sd(MCSVD_2_cases$DBP, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$DBP, na.rm=T),"(", sd(MCSVD_2_controls$DBP, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$DBP, MCSVD_2_controls$DBP, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[8,]	<- c("SBP_mean(s.d)", paste(mean(MCSVD_2_cases$SBP, na.rm=T),"(", sd(MCSVD_2_cases$SBP, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$SBP, na.rm=T),"(", sd(MCSVD_2_controls$SBP, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$SBP, MCSVD_2_controls$SBP, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[9,]	<- c("Fast_Glucose_mean(s.d)", paste(mean(MCSVD_2_cases$Fast_Glucose, na.rm=T),"(", sd(MCSVD_2_cases$Fast_Glucose, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$Fast_Glucose, na.rm=T),"(", sd(MCSVD_2_controls$Fast_Glucose, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$Fast_Glucose, MCSVD_2_controls$Fast_Glucose, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[10,]	<- c("HDL_mean(s.d)", paste(mean(MCSVD_2_cases$HDL, na.rm=T),"(", sd(MCSVD_2_cases$HDL, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$HDL, na.rm=T),"(", sd(MCSVD_2_controls$HDL, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$HDL, MCSVD_2_controls$HDL, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[11,]	<- c("LDL_mean(s.d)", paste(mean(MCSVD_2_cases$LDL, na.rm=T),"(", sd(MCSVD_2_cases$LDL, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$LDL, na.rm=T),"(", sd(MCSVD_2_controls$LDL, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$LDL, MCSVD_2_controls$LDL, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[12,]	<- c("TG_mean(s.d)", paste(mean(MCSVD_2_cases$TG, na.rm=T),"(", sd(MCSVD_2_cases$TG, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$TG, na.rm=T),"(", sd(MCSVD_2_controls$TG, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$TG, MCSVD_2_controls$TG, var.equal = F, paired=F)$p.value)
DISCREPTIVE_MCSVD2[13,]	<- c("BMI_mean(s.d)", paste(mean(MCSVD_2_cases$BMI, na.rm=T),"(", sd(MCSVD_2_cases$BMI, na.rm=T), ")", sep=""),
	paste(mean(MCSVD_2_controls$BMI, na.rm=T),"(", sd(MCSVD_2_controls$BMI, na.rm=T), ")", sep=""), 
	t.test(MCSVD_2_cases$BMI, MCSVD_2_controls$BMI, var.equal = F, paired=F)$p.value)
# Now categorical variables
DISCREPTIVE_MCSVD2[14,]	<- c("Females_N(percentage)", paste(length(which(MCSVD_2_cases$Sex == 1)),"(", length(which(MCSVD_2_cases$Sex == 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$Sex == 1)),"(", length(which(MCSVD_2_controls$Sex == 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD2_All$Sex, MCSVD2_All$MCSVD_2))$p.value)
DISCREPTIVE_MCSVD2[15,]	<- c("LBI_status_N(percentage)", paste(length(which(MCSVD_2_cases$LBI== 1)),"(", length(which(MCSVD_2_cases$LBI== 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$LBI== 1)),"(", length(which(MCSVD_2_controls$LBI== 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	"NA") # It doesnt make sense to test association LBI, since we already know it will be 0 in min-CSVD
DISCREPTIVE_MCSVD2[16,]	<- c("HT_status_N(percentage)", paste(length(which(MCSVD_2_cases$HT_status== 1)),"(", length(which(MCSVD_2_cases$HT_status== 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$HT_status== 1)),"(", length(which(MCSVD_2_controls$HT_status== 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD2_All$HT_status, MCSVD2_All$MCSVD_2))$p.value)
DISCREPTIVE_MCSVD2[17,]	<- c("Diabetes_status_N(percentage)", paste(length(which(MCSVD_2_cases$Diabetes== 1)),"(", length(which(MCSVD_2_cases$Diabetes== 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$Diabetes== 1)),"(", length(which(MCSVD_2_controls$Diabetes== 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD2_All$Diabetes, MCSVD2_All$MCSVD_2))$p.value)
DISCREPTIVE_MCSVD2[18,]	<- c("CVD_status_N(percentage)", paste(length(which(MCSVD_2_cases$CVD== 1)),"(", length(which(MCSVD_2_cases$CVD== 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$CVD== 1)),"(", length(which(MCSVD_2_controls$CVD== 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD2_All$CVD, MCSVD2_All$MCSVD_1))$p.value)
DISCREPTIVE_MCSVD2[19,]	<- c("AntiHT_drug_status_N(percentage)", paste(length(which(MCSVD_2_cases$AntiHT_drug== 1)),"(", length(which(MCSVD_2_cases$AntiHT_drug== 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$AntiHT_drug== 1)),"(", length(which(MCSVD_2_controls$AntiHT_drug== 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD2_All$AntiHT_drug, MCSVD2_All$MCSVD_2))$p.value)
DISCREPTIVE_MCSVD2[20,]	<- c("Lipid_low_drug_status_N(percentage)", paste(length(which(MCSVD_2_cases$lipid_low_drug== 1)),"(", length(which(MCSVD_2_cases$lipid_low_drug== 1))/length(MCSVD_2_cases[,1]), ")", sep=""),
	paste(length(which(MCSVD_2_controls$lipid_low_drug== 1)),"(", length(which(MCSVD_2_controls$lipid_low_drug== 1))/length(MCSVD_2_controls[,1]), ")", sep=""), 
	fisher.test(table(MCSVD2_All$lipid_low_drug, MCSVD2_All$MCSVD_2))$p.value)

names(DISCREPTIVE_MCSVD2) <-c("Variables","Extensive_MCSVD", "Minimal_MCSVD", "pvalue")
write.table(DISCREPTIVE_MCSVD2, "MCSVD2_Descriptives.out", quote=F, row.names=F, col.names=T, sep= "\t")


cat("\n\t\tThank you for using our script\n\n")

