

require(motifDiverge)
require(Biostrings)
require(MotifDb)
library("foreach")
library("doMPI")

## Without this code, the MPI session may fail to shut down properly in case an
## error occurs, and then the script won't terminate until it hits the walltime
options(error=quote(assign(".mpi.err", FALSE, env = .GlobalEnv)))

cl <- startMPIcluster()

## Tell foreach() to parallelize on the MPI cluster
registerDoMPI(cl)

# enh.hg.file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files/hg19_MER20.fa"
# enh.mm.file = "/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/fasta_files/mm9_MER20.fa"
enh.hg.file = "/Users/Abin/Google Drive/GD_transfer/data/hg19_MER20.fa"
enh.mm.file = "/Users/Abin/Google Drive/GD_transfer/data/mm9_MER20.fa"


## ====== FUNCTIONS =====

get_index <- function(TF_name) { 
  geneSymbol.rows = grep (TF_name, values(MotifDb)$geneSymbol, ignore.case=TRUE)
  organism.rows = grep ('Hsapien', values(MotifDb)$organism, ignore.case=TRUE)
  source.rows = grep ('JASPAR_2014', values(MotifDb)$dataSource, ignore.case=TRUE)
  tf_name.hu.jaspar.rows = intersect(geneSymbol.rows, intersect(organism.rows, source.rows))
  
  return(tf_name.hu.jaspar.rows[1])

}

## ====== LOAD DATA ======
enh.hg = readDNAStringSet(enh.hg.file)
enh.mm = readDNAStringSet(enh.mm.file)
enh.hg = enh.hg[1:4]
enh.mm = enh.mm[1:4]

## ====== LOAD TF PROFILE ======
tf_list = c('YY1','CTCF', 'CEBPb', 'USF1','p53','ETS1') #Pol-II, 'PRMT1/4',p300, 'FOXO1A', 'HoxA-11', 'PGR', 'TFIF'   not found 
sapply(tf_list, get_index)

index      = get_index('CEBPb')
pspm       = MotifDb[index][[1]]
pssm       = pspmToPssm(pspm)
pspm       = pspmToPssm(pspm,return.pspm=TRUE)$pspm

## ------------------------------------------------------------------------
bg  = colSums(alphabetFrequency(c(enh.hg,enh.mm))[,1:4])
bg  = bg/sum(bg)
cut = scoreCutFromErr(err=.05, pssm=pssm, pspm=pspm, bg=bg, type="type1")

## ------------------------------------------------------------------------
evo.mod.file = system.file( "extdata", "primates.mod", package="motifDiverge")

## ------------------------------------------------------------------------
require(rphast)
neutral.mod = get.neutralMod(evo.mod.file)  #- point ucscURL to model
neutral.mod = mod.backgd.tm(neutral.mod,bg) #- adjust background
pspm.mods   = get.modList(neutral.mod,pspm) #- models for pspm columns

start = Sys.time() 
pars.model  = cbernEstimateModelPars(  seqs.x    = enh.mm   ,
                                       seqs.y    = enh.hg   ,
                                       pssm      = pssm     ,
                                       pspm      = pspm     ,
                                       cut.fw    = cut      ,
                                       cut.rc    = cut      ,
                                       indep     = FALSE    ,
                                       useCounts = FALSE    ,
                                       modList   = pspm.mods)
end = Sys.time() - start 

## ------------------------------------------------------------------------  

pstart = Sys.time()                                      
output <- foreach (ind = seq(1,length(enh.mm) , by = 2) , .combine = "rbind", .packages=c('Biostrings',"motifDiverge")) %dopar% {

  this_mouse = enh.mm[ind:(ind+1)]
  this_human = enh.hg[ind:(ind+1)]

  pars.model  = cbernEstimateModelPars( seqs.x    = this_mouse,
                                        seqs.y    = this_human,
                                        pssm      = pssm     ,
                                        pspm      = pspm     ,
                                        cut.fw    = cut      ,
                                        cut.rc    = cut      ,
                                        indep     = FALSE    ,
                                        useCounts = FALSE    ,
                                        modList   = pspm.mods) #bottle neck 
  
  print(pars.model )
  
}

pend = Sys.time() - pstart

cat("for ", length(enh.mm), ":")
cat("time w/o MPI:", end , "\n")
cat("time w MPI:", pend , "\n")

pars.model = output 
## ------------------------------------------------------------------------

# filter if necessary... BUT REMEMBER TO REALIGN THE HOMOLOGOUS SEQUENCES !!!!!! 
pars.model_filter <- pars.model  

#store enrichments 
store_enrichments <- matrix(0, nrow = dim(pars.model_filter)[1], ncol = 3)
store_depletion <- matrix(0, nrow = dim(pars.model_filter)[1], ncol = 3)

for (i in 1:dim(pars.model_filter)[1]){
  # print(i)
  x = pars.model_filter[i,1:6]
  
  result_enrich = tryCatch({
    pcbern(x[1],x[2],x[3],x[4],x[5],x[6], lower.tail=F)}, 
  error = function(e) {
    print(paste("MY_ERROR:  ",e, i))
    return(-1)})
  
  result_deplete = tryCatch({
    pcbern(x[1],x[2],x[3],x[4],x[5],x[6], lower.tail=T)}, 
    error = function(e) {
      print(paste("MY_ERROR:  ",e,i))
      return(-1)})
  
  print(result_enrich) 
  store_enrichments[i,1] = result_enrich
  store_enrichments[i,2] = names(enh.hg)[i]
  store_enrichments[i,3] = names(enh.mm)[i]
  store_depletion[i,1] = result_deplete
  store_depletion[i,2] = names(enh.hg)[i]
  store_depletion[i,3] = names(enh.mm)[i]
}

colnames(store_enrichments) = c('p_value','hg19','mm9')
colnames(store_depletion) = c('p_value','hg19','mm9')

## ------------------------------------------------------------------------
# setwd("/dors/capra_lab/abraha1/projects/transposable_elements/data/tf_dynamics/comparative/MER20/")
write.table(as.matrix(store_enrichments), file = "enrich_CEBPb_MER20_hg19_mm9.out", append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"))

write.table(as.matrix(store_depletion), file = "depleted_CEBPb_MER20_hg19_mm9.out", append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = T, qmethod = c("escape", "double"))

## -----------------------------------------------------------------------------
## Cluster shutdown
## -----------------------------------------------------------------------------

closeCluster(cl)
mpi.quit()
