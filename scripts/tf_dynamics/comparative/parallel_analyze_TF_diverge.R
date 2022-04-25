### run motifDiverge to evaluate homologous sequences for enrichment/depeletion of a given TFBS 
###
###       INPUTS: 
###         - <speciesA>.fa      :  fasta file for species A; ***note: for homologs, each sequence in <speciesA>.fa should pair in order with <speciesB>.fa
###         - <speciesB>.fa      :  fasta file for species B; ***note: for homologs, each sequence in <speciesA>.fa should pair in order with <speciesB>.fa
###         - <tf_pasm>.tsv      :  tf_pspm_file: path for tab seperated PSPM for TF of interest; 
###                              :  four rows (A,C,G,T) and num columns = num bases in motif
###                              :  ignores first row when reading this file in 
###         - "/store/outputs/"  :  output_file_path: file path for output txt to be written to 
###         - "output_header"    : output_header: file path for output txt to be written to 
###
###
###       NOTE: 
###         - uses evolutionary model for primates from PHAST
### Abin Abraham 
### 2018-03-09 09:25:26


require(motifDiverge)
require(Biostrings)
require(MotifDb)
library("foreach")
library("doParallel")
library("qvalue")
# library("dpylr")


### ====== FUNCTIONS =====

get_index <- function(TF_name) { 
  geneSymbol.rows = grep (TF_name, values(MotifDb)$geneSymbol, ignore.case=TRUE)
  organism.rows = grep ('Hsapien', values(MotifDb)$organism, ignore.case=TRUE)
  source.rows = grep ('JASPAR_2014', values(MotifDb)$dataSource, ignore.case=TRUE)
  tf_name.hu.jaspar.rows = intersect(geneSymbol.rows, intersect(organism.rows, source.rows))
  
  return(tf_name.hu.jaspar.rows[1])

}

### ======= MAIN ======

## check input args 
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("no inputs passed to script...", call.=FALSE)
} else if (length(args)<5) {
  # default output file
  stop("not enough inputs passed to script...", call.=FALSE)
}

# store input args into variables 
# speciesA_file = "/Users/abin-personal/Google Drive/GD_transfer/data/hg19_MER20.fa"
# speciesB_file = "/Users/abin-personal/Google Drive/GD_transfer/data/mm9_MER20.fa"
# output_file_path = "/Users/abin-personal/Google Drive/GD_transfer/output/test"
# output_header = format(Sys.time(), "%X")

speciesA_file = args[1]
speciesB_file = args[2]
tf_pspm_file = args[3]
output_file_path = args[4]
output_header = args[5]

# set up cores 
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores, type="FORK")  
# registerDoParallel(cores=no_cores)
registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

### ====== LOAD DATA ======
seqA = readDNAStringSet(speciesA_file)
seqB = readDNAStringSet(speciesB_file)

if (length(seqA) %% 2  == 1) {
  seqA = seqA[1:length(seqA)-1]
  seqB = seqA[1:length(seqB)-1]
  print("last input sequence was removed from both species b/c there is an ODD number of sequences")
} 

# rand_seq = sample(1:10198, 10, replace=F)
# seqA = seqA[rand_seq]
# seqB = seqB[rand_seq]

seqA = seqA
seqB = seqB
## ====== LOAD TF PROFILE ======
#tf_list = c('YY1','CTCF', 'CEBPb', 'USF1','p53','ETS1') #Pol-II, 'PRMT1/4',p300, 'FOXO1A', 'HoxA-11', 'PGR', 'TFIF'   not found 

get_pspm = as.matrix(read.table(tf_pspm_file, header=FALSE, sep="\t", skip=1))
rownames(get_pspm) = c("A","C","G","T")
colnames(get_pspm) = seq(1, dim(get_pspm)[2])
pspm = get_pspm

# index = get_index('CEBPb')
# pspm = MotifDb[index][[1]]
pssm = pspmToPssm(pspm)
pspm = pspmToPssm(pspm,return.pspm=TRUE)$pspm

## ====== CALC BACKGROUND FREQUENCES  ======
bg  = colSums(alphabetFrequency(c(seqA,seqB))[,1:4])
bg  = bg/sum(bg)
cut = scoreCutFromErr(err=.05, pssm=pssm, pspm=pspm, bg=bg, type="type1")

# load evolutionary model 
evo.mod.file = system.file( "extdata", "primates.mod", package="motifDiverge")

## ====== CALC MODEL PARAMETERS  ======
require(rphast)
neutral.mod = get.neutralMod(evo.mod.file)  #- point ucscURL to model
neutral.mod = mod.backgd.tm(neutral.mod,bg) #- adjust background
pspm.mods   = get.modList(neutral.mod,pspm) #- models for pspm columns

pstart = Sys.time()   
                                   
output <- foreach (ind = seq(1,length(seqB) , by = 2) , .combine = "rbind", .packages=c('Biostrings',"motifDiverge")) %dopar% {

  this_seqA = seqA[ind:(ind+1)]
  this_seqB= seqB[ind:(ind+1)]

  pars.model  = cbernEstimateModelPars( seqs.x    = this_seqA,
                                        seqs.y    = this_seqB,
                                        pssm      = pssm     ,
                                        pspm      = pspm     ,
                                        cut.fw    = cut      ,
                                        cut.rc    = cut      ,
                                        indep     = FALSE    ,
                                        useCounts = FALSE    ,
                                        modList   = pspm.mods) #bottle neck 
  
  print(pars.model)
  
}

pend = Sys.time() - pstart

# print summary 
cat("__________________________\n")
cat("for ", length(seqB), ":\n")
cat("time w multi-threading:", pend , "\n")

pars.model = output

### ====== CALC ENRICH/DEPLETION  ======

# filter if necessary... BUT REMEMBER TO REALIGN THE HOMOLOGOUS SEQUENCES !!!!!! 
pars.model_filter <- pars.model  

# intialize store enrichments 
store_enrichments <- matrix(0, nrow = dim(pars.model_filter)[1], ncol = 3)
store_depletion <- matrix(0, nrow = dim(pars.model_filter)[1], ncol = 3)

# find enrichments and depeletion 
for (i in 1:dim(pars.model_filter)[1]){
  # print(i)
  x = pars.model_filter[i,1:6]
  
  result_enrich = tryCatch({
    pcbern(x[1],x[2],x[3],x[4],x[5],x[6], lower.tail=F)}, 
  error = function(e) {
    # print(paste("MY_ERROR:  ",e, i))
    return(NA)})
  
  result_deplete = tryCatch({
    pcbern(x[1],x[2],x[3],x[4],x[5],x[6], lower.tail=T)}, 
    error = function(e) {
      # print(paste("MY_ERROR:  ",e,i))
      return(NA)})
  
  # print(result_enrich) 
  store_enrichments[i,1] = result_enrich
  store_enrichments[i,2] = names(seqA)[i]
  store_enrichments[i,3] = names(seqB)[i]
  store_depletion[i,1] = result_deplete
  store_depletion[i,2] = names(seqA)[i]
  store_depletion[i,3] = names(seqB)[i]
}

### ====== ADD Q-VALUES  ======

store_qvalues_enrich = tryCatch({
  # qvalue(p = as.numeric(store_enrichments[,1]), lambda=seq(0.2,0.8,0.01))$qvalues},
  qvalue(p = as.numeric(store_enrichments[,1]), lambda=0, fdr.level =0.05)$qvalues},
  
  error = function(e) { 
    print(paste("error in qvalue calculation: ", e))
    # stop(call=TRUE)
    return(NA)})  

store_qvalues_deplete = tryCatch({
  # qvalue(p = as.numeric(store_depletion[,1]), lambda=seq(0.2,0.8,0.01))$qvalues},
  qvalue(p = as.numeric(store_depletion[,1]), lambda=0, fdr.level =0.05)$qvalues},
  
  error = function(e) { 
    print(paste("error in qvalue calculation: ", e))
    # stop(call=TRUE)
    return(NA)})


# store_qvalues_enrich = qvalue(p = as.numeric(store_enrichments[,1]), lambda=seq(0.2,0.8,0.01))$qvalues 
# store_qvalues_deplete = qvalue(p = as.numeric(store_depletion[,1]), lambda=seq(0.2,0.8,0.01))$qvalues 

store_depletion = cbind(store_qvalues_deplete, store_depletion)
store_enrichments = cbind(store_qvalues_enrich, store_enrichments)

colnames(store_enrichments) = c('enrich_qvalue', 'enrich_p_value','seqA','seqB')
colnames(store_depletion) = c('deplete_qvalue', 'deplete_p_value','seqA','seqB')

### ====== FORMAT OUTPUT  ======


### ====== WRITE OUTPUTS  ======
# setwd("/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/")
enrich_output_file_path = file.path(output_file_path, paste("output_", output_header, "__enrichment_motifDiverge.tsv", sep=""))
deplete_output_file_path = file.path(output_file_path, paste("output_", output_header, "__deplete_motifDiverge.tsv", sep=""))

write.table(as.matrix(store_enrichments), file = enrich_output_file_path, append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"))

write.table(as.matrix(store_depletion), file = deplete_output_file_path, append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"))


print(paste("ouputs written to: ", enrich_output_file_path))
print(paste("ouputs written to: ", deplete_output_file_path))
print("DONE!")
