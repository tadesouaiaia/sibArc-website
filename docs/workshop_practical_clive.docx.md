# Day 2 Morning practical: Cross-ancestry PRS analyses using BridgePRS

## Table of contents

1. [Key learning outcomes](#key-learning-outcomes)
2. [BridgePRS input data](#BridgePRS-input-data)
3. [Passing arguments to BridgePRS - example run](#passing-arguments-to-bridgePRS---example-run)
4. [BridgePRS output](#bridgePRS-output)
5. [Using BridgePRS without target summary statistics](#using-bridgePRS-without-target-summary-statistics)
6. [Further analyses with BridgePRS](#further-analyses-with-bridgePRS)
7. [Using BridgePRS SNP weights](#using-bridgePRS-snp-weights)

---

## Key learning outcomes
In the previous lecture we covered the modelling used by
BridgePRS. Here we will use the BridgePRS software to apply the method.
The aim of this practical is to provide you with a basic understanding
and some experience of running BridgePRS software. After completing
this practical, you should be able to:

* perform cross-population analyses using BridgePRS
  * understand what input files are required by BridgePRS
  * set up configuration files to pass arguments to BridgePRS
* interpret output from BridgePRS
  * use BridgePRS output to calculate PRS in unseen samples using SNP weights

<a href="#top">[Back to table of contents](#table-of-contents)</a>

## BridgePRS input data
In the BridgePRS directory there is a data folder which we will use in
this practical. View the data directory
```
$ ls -l data
total 5368
drwxr-xr-x  73 hoggac01  staff     2336 26 Jul 16:00 1000G_sample
-rw-r--r--   1 hoggac01  staff      469 12 Aug 14:35 afr.config
-rw-r--r--   1 hoggac01  staff      469 12 Aug 14:35 eas.config
-rw-r--r--   1 hoggac01  staff      410 12 Aug 14:35 eur.config
drwxr-xr-x   5 hoggac01  staff      160 14 Jul 17:22 pop_AFR
drwxr-xr-x   5 hoggac01  staff      160 14 Jul 17:22 pop_EAS
drwxr-xr-x   5 hoggac01  staff      160  7 Aug 20:30 pop_EUR
-rw-r--r--   1 hoggac01  staff   200376  7 Aug 21:30 qc_snplist.txt
```

The `pop_AFR, pop_EAS, pop_EUR` folders contain simulated genotype, phenotype and GWAS
summary statistics representative of
Europeans, East Asians and Africans for input to BridgePRS. Each pop_\*
folder is split into summary statistics, and
individual level genotype and phenotype folders, e.g.
```
$ ls -l data/pop_AFR/
total 0
drwxr-xr-x  68 hoggac01  staff  2176 14 Jul 17:22 genotypes
drwxr-xr-x   4 hoggac01  staff   128 14 Jul 17:22 phenotypes
drwxr-xr-x  68 hoggac01  staff  2176 12 Aug 11:02 sumstats
```
Look at each directory e.g. `ls pop_AFR/genotypes`. There are two sets of summary
statistics in each `sumstats` folder from the analysis of the same simulated
continuous phenotype, the "half" files were generated using half the
sample size. For computation speed the summary statistics
only have a small subset of SNPs, 19k-20k genomewide, the following
commands count the number of lines in each set of summary statistics
```
zcat data/pop_EAS/sumstats/EAS.chr* | wc -l
zcat data/pop_EUR/sumstats/EUR.chr* | wc -l
zcat data/pop_AFR/sumstats/AFR.chr* | wc -l
```
or on a Mac
```
gzcat data/pop_EAS/sumstats/EAS.chr* | wc -l
gzcat data/pop_EUR/sumstats/EUR.chr* | wc -l
gzcat data/pop_AFR/sumstats/AFR.chr* | wc -l
```
`zcat`/`gzcat` are used as the files are gzipped. `|` "pipes" the
output to `wc -l` which counts the number of lines.

Results in these files are only shown for SNPs with MAF>0 and are a
random subset of HapMap SNPs.
### Questions?
* Which population has the most polymorphic SNPs?
* Which population has the least polymorphic SNPs?
* Do you think there is a bias in the SNP ascertainment?
* Which population would you expect to have the most polymorphic SNPs
  if SNPs were ascertained randomly from 1000G?

Take a look a look at the files, e.g.
```
zcat data/pop_AFR/sumstats/AFR.chr19.glm.linear.gz | head 
zcat data/pop_AFR/sumstats/AFR_half.chr19.glm.linear.gz | head 
```
again use `gzcat` on a Mac.

In the OBS_CT column you'll see can that the "_half" summary
statistics files have half the sample size, 10,000 compared to 20,000
for the EAS and AFR populations and 40,000 compared to 80,000 for the
EUR population. The same SNPs are contained in each set of summary statistics.

The `phenotypes` folders has two files: "test" and "validation" with
IDs, the outcome phenotype and covariates. "Test" data is used to
estimate the PRS and "validation" data is used only to assess model performance,
it is not used to estimate the PRS.

The `genotypes` folders are in `plink` binary format and are split by
chromosome. These folders contain the genetic data for both test and
validation individuals in the phenotypes folder.

Test data is only used for individuals with both genotype and phenotype
information. Similarly model performance metrics is only used validation samples with both
genotype and phenotype information, however, predictions are generated
for all validation samples with genotype data.

<a href="#top">[Back to table of contents](#table-of-contents)</a>

## Passing arguments to BridgePRS - example run
Example run of BridgePRS:
```
./bridgePRS pipeline go -o out/ --config_files data/eas.config data/eur.config --fst 0.11 --phenotype y --cores 4 --restart
```
Arguments can be passed to BridgePRS both directly on the command-line and in
config files. The command-line arguments used above are:
* `-o`  Output folder
* `--config_files`  Configuration files, see below
* `--fst` Fst between the base and target populations used in the analysis, used to specify a prior
distribution. The Fst between EAS and EUR is 0.11
* `--phenotype`  The column label of the phenotype in the test and validation
files, e.g. `EAS_valid.dat`
* `--cores`  The number of cores used in the analysis.

#### Fixation index Fst
Fst is a measure of genetic distance between two populations. Details on the definition and calculation of Fst can be found on [Wikipedia](https://en.wikipedia.org/wiki/Fixation_index). Given genetic data from two populations Fst can be calculated in `plink`, for details see https://www.cog-genomics.org/plink/1.9/basic_stats#fst.

### The \*.config files
config files are a neat way to store population
specific arguments. For a standard two population analysis
two config files are required, by default the first config file is for the target
population and the second is for the base population.   
.config files tell BridgePRS where to find the required input files 
and the column headers of the summary statistics files for a
population data set, take a look, e.g.
```
cat data/eas.config 
POP=EAS
LDPOP=EAS
LD_PATH=1000G_sample
SUMSTATS_PREFIX=pop_EAS/sumstats/EAS.chr
SUMSTATS_SIZE=20000
#SUMSTATS_PREFIX=pop_EAS/sumstats/EAS_half.chr
#SUMSTATS_SIZE=10000
SUMSTATS_SUFFIX=.glm.linear.gz
GENOTYPE_PREFIX=pop_EAS/genotypes/chr
PHENOTYPE_FILE=pop_EAS/phenotypes/EAS_test.dat  
VALIDATION_FILE=pop_EAS/phenotypes/EAS_valid.dat 
SNP_FILE=qc_snplist.txt
SSF-P=P
SSF-SNPID=ID
SSF-BETA=BETA
SSF-REF=REF
SSF-ALT=A1
```
This config file contains all possible arguments that can be used in
config files. config files use the same argument names as the
command-line arguments but in uppercase, and use "=" instead of a space
between the argument name and the argument being passed.

The `POP` argument simply labels the population used in this .config
file for output.

A full list of arguments can be found [here](https://www.bridgeprs.net/guide_args/).

### Estimating linkage disequilibrium (LD)
BridgePRS requires individual level genetic data in `plink` binary
format to estimate linkage disequilibrium (LD) in the populations
which produced the GWAS summary statistics. The genotype test and
validation data could be used, e.g. data here
`pop_EUR/genotypes/`. If these data are small, less than 500 samples, or
are not representative of the GWAS population, we provide 1000 Genomes
(1000G) data to estimate LD. Suitable 1000G data for this analysis
is in the 1000G_sample folder
for the small subset of SNPs used in these examples.

The `.config` files point to
the folder with reference LD data by the `LD_PATH` argument in the
config file. LD reference data is available for the five 1000G super populations (abbreviations
required to use in brackets): East Asian (EAS), South Asian (SAS),
European (Eur), African (AFR) and American (AMR).

**For real data analyses 1000G reference data for larger subsets of SNPs can be
downloaded [here](https://www.bridgeprs.net)**

### Questions?
* Can you work out what the other command line arguments are doing?
* Can you work out what the other config file arguments are doing?

<a href="#top">[Back to table of contents](#table-of-contents)</a>

## BridgePRS output
The main output is in the folder `out/prs-combined_EAS-EUR/`. First
view the output summary plot
```
evince  out/prs-combined_EAS-EUR/bridgePRS-combo.pdf
```
on a Mac use `open` instead of `evince`.

The barplot at the top shows the variance explained (R2) by the
four PRS models BridgePRS estimates. The weighted model is BridgePRS
estimated "best" PRS.


### The models estimated by BridgePRS
The three separate target population PRS estimated by BridgePRS are:
* **Stage 2 model** -- estimated using target population data with prior effect-size distribution from the base (European) population
* **Stage 1 model** -- estimated using only the target (Non-European) data
* **Stage1+2 model** -- estimated using both stage 1 and stage 2 models

Each of these three models are given weights corresponding to how well
they fit the test data. These weights are then used to combine the PRS
to give the single weighted combined PRS.

The weighted combined PRS should typically be used.  The models,
stage1, stage2 and stage1+2, should not be used unless users have a
strong prior belief that a particular model is better. The hypotheses
of the three models are:
* **Stage 2 model** reflects the belief that the target population GWAS is only
  informative when used with the base population GWAS
* **Stage 1 model** reflects the belief that the target population GWAS is
  informative and the base population GWAS gives no addition
  information.
* **Stage 1+2 model** reflects the belief both the base and target population
  GWAS contribute independent information

Look at the following output file
```
cat out/prs-combined_EAS-EUR/EAS_weighted_combined_var_explained.txt
```
### Questions?
* Which plot in the summary plot was constructed from this output
  file?
* How do the Manhattan plots of the base and target populations
  compare?

`EAS_weighted_combined_preds.dat` has PRS predictions for samples in
the validation data using all four models: stage1, stage2, stage1+2
and weighted. `EAS_weighted_combined_snp_weights.dat` has the SNP
weights for the combined to allow this model to be applied to other
samples.

<a href="#top">[Back to table of contents](#table-of-contents)</a>

## Using BridgePRS without target summary statistics
Often GWAS summary statistics are only available in one
population. BridgePRS can use these summary statistics and optimise
them to estimate a PRS for another target population given individual
level from the target population. Here is an example using the
`data/eur_eas.config` config file
```
./bridgePRS prs-single run -o out_single/ --config_file data/eur_eas.config --phenotype y --cores 4
```
Look at `data/eur_eas.config`.
### Questions?
* What GWAS summary data is used?
* What test data is used for model optimisation?
* What validation data is used?

Results of interest are written to the folder
`out_single/prs-single_EAS/quantify/`. Model performance is shown in
the file `EAS_quantify_var_explained.txt` and plotted in ....

See how these results compare with the previous analysis which
included EAS GWAS summary statistics.
```
cat out_single/prs-single_EAS/quantify/EAS_quantify_var_explained.txt
cat out/prs-combined_EAS-EUR/EAS_weighted_combined_var_explained.txt
```
This single summary statistic analysis is equivalent to the stage 2
analysis previously but with all the weight on the EUR prior.

### Question?
* How does prediction using only EUR summary statistics compare with
those which include information from the EAS summary statistics?

<a href="#top">[Back to table of contents](#table-of-contents)</a>

## Further analyses with BridgePRS
### African analysis
Run BridgePRS again to estimate PRS in Africans using
`afr.config`. Note, you should also change the command line `fst`
argument to match the Fst between Africans and Europeans, use 0.15.

### Qustions?
* Which population, EAS or AFR, has the best prediction?
* What are the reasons for the differences in prediction between the populations?
  
### Analyses with other GWAS summary statistics
For each population the config files contain commented out links to GWAS summary
statistics of the same phenotype using half the same size: 40k for EUR
and 10k for both EAS and AFR.

Edit `eas.config` to use the EAS 10k
GWAS summary statistics. To run the analysis write results to a new
output directory e.g. `out_half_target`. Run the similar analysis for
African samples by editing  `afr.config`. Compare with previous results using the
10k EAS and EAS GWAS. Compare EAS and AFR results.

Check you've run the analyses using the correct GWAS summary
statistics, e.g.
```
less less out_half_target/logs/bridgePRS.eas-eur.pipeline.go.log
less less out_half_target/logs/bridgePRS.afr-eur.pipeline.go.log
```
or
```
grep Sumstats out_half_target/logs/bridgePRS.eas-eur.pipeline.go.log
grep Sumstats out_half_target/logs/bridgePRS.afr-eur.pipeline.go.log
```

If you have made a mistake, correct and run again using the
`--restart` flag which deletes the previously generated results.

### Questions?
* How has using the less well powered EAS and AFR GWAS affected the predictive
  accuracy of the BridgePRS models?
* How do AFR and EAS results compare?

#### Analyses with smaller EUR GWAS summary statistics
Run BridgePRS analysis using using the 40k EUR GWAS (use `EUR_half`) as base and 20k EAS GWAS as target.
You will need to edit `eur.config` and `eas.config` to use the correct summary statistics and pass the
correct GWAS sample size. Write to results to a new directory e.g. `out_half_eur`.
Run the equivalent analysis for AFR.

### Questions?
* How has using the less well powered EUR GWAS affected the predictive
  accuracy of the BridgePRS models?
* How do AFR and EAS results compare?

<a href="#top">[Back to table of contents](#table-of-contents)</a>

## Using BridgePRS SNP weights
The SNP weights in `EAS_weighted_combined_snp_weights.dat` can be used to make
predictions in other samples for which we have overlapping genotype data. We demonstrate
this using `plink` and data in `data/pop_EAS/genotypes/`. The genotype data is split
by chromosome, therefore predictions are estimated for each chromosome separately and
then combined. The following bash commands estimate the per chromosome predictions using plink,
we first make a directory to write these predictions to
```
mkdir out/preds
for chr in {1..22}
do
    plink --score out/prs-combined_EAS-EUR/EAS_weighted_combined_snp_weights.dat 1 2 4 sum \
          --bfile data/pop_EAS/genotypes/chr$chr \
          --out out/preds/EAS_chr$chr
done
```
We then read these predictions into R to combine
```
R
library(data.table)
tmp <- fread(paste0('out/preds/EAS_chr1.profile'))
pred <- data.frame( tmp$IID, tmp$SCORESUM )
colnames(pred) <- c('IID','Score')
for( chr in 2:22 ){
    tmp <- fread(paste0('out/preds/EAS_chr',chr,'.profile'))
    pred[,2] <- pred[,2] + tmp$SCORESUM
}

# Check predictions are the same as those produced directly by BridgePRS
valid <- fread('out/prs-combined_EAS-EUR/EAS_weighted_combined_preds.dat')
# Match IDs
ptr <- match( valid$IID, pred$IID )
plot( pred$Score[ptr], valid$Weighted, main="Individual PRS"
      xlab="PRS estimated by BridgePRS",
      ylab="PRS re-estimated from SNP weights" )
```
Predictions are the same \+/\- rounding error and a constant.

<a href="#top">[Back to table of contents](#table-of-contents)</a>
