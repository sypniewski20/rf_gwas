suppressMessages(library(randomForest))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(doParallel))

### function for parallel RF classifier ###

parRF_classifier <- function(pheno,
                             gt,
                             ntree,
                             mtry,
                             ncpu,
                             sample_size) {
  foreach(ntree = rep(ntree/ncpu, ncpu), 
          mtry = mtry,
          .combine = randomForest::combine, 
          .multicombine = T,
          .packages = 'randomForest') %dopar%
    randomForest(x = gt, 
                 y = pheno, 
                 importance = TRUE,
                 proximity = TRUE, 
                 ntree = ntree, 
                 mtry = mtry, 
                 strata = pheno, 
                 sampsize = sample_size)
}

### specify cpus ###

ncpu <- 10

### specify output name

output_path <- 'bcs_obesity_panel'
dir.create(file.path(output_path))


### create dir for RDS files ###
dir.create(file.path(paste(output_path,'/RDS', sep = '')))


### path to genotypes

genotypes_path <- 'bcs_obesity_panel/gt_data_corrected.tsv'

print('loading genotypes', quote = F)
gt_data_corrected <- fread(genotypes_path, 
                 header = T,
                 nThread = 18) 

### path to phenotype 

phenotype_path <- 'bcs_pheno.tsv'

pheno <- read.table(phenotype_path, header = T) %>%
  pull() %>%
  as.factor()

### path to optimization results ###

opt_path <- 'bcs_obesity_panel/results_optimization.tsv'

results_optimization <- fread(opt_path) %>%
  filter(cor > 0.95) %>%
  filter(OOB_ER_1 == min(OOB_ER_1)) %>%
  filter(ntree == min(ntree))

##### hyperparameters ######

## set number of trees

ntree <- results_optimization$ntree

## set number of loci per tree

mtry <- results_optimization$mtry

## set sampling size

sample_size <- rep(
  round(sum(pheno == 1)*2/3), 2
)

cl <- makePSOCKcluster(ncpu)
registerDoParallel(cl)

for (i in 1:3) {
  print(paste('rf_all_', i, sep = ''), quote = F)
  rf_all <- parRF_classifier(pheno = pheno,
                             gt = gt_data_corrected,
                             ntree = ntree,
                             mtry = mtry,
                             ncpu = ncpu,
                             sample_size = sample_size)

  write_rds(rf_all, file = paste(output_path,'/RDS/rf_all_', i, '.rds', sep = ''))
}

stopCluster(cl)

rf_all_1 <- readRDS(paste(output_path,'/RDS/rf_all_1', '.rds', sep = ''))
rf_all_2 <- readRDS(paste(output_path,'/RDS/rf_all_2', '.rds', sep = ''))
rf_all_3 <- readRDS(paste(output_path,'/RDS/rf_all_3', '.rds', sep = ''))

importance_rf_all_1 <- data.frame(importance(rf_all_1,type=1)) %>%
  rename(importance = MeanDecreaseAccuracy)

importance_rf_all_2 <- data.frame(importance(rf_all_2,type=1)) %>%
  rename(importance = MeanDecreaseAccuracy)

importance_rf_all_3 <- data.frame(importance(rf_all_3,type=1)) %>%
  rename(importance = MeanDecreaseAccuracy)

# ############################################################################################################################################
# ############################################################################################################################################
# 
# The predictive ability of classification trees is measured by the out-of-bag error rate.  An error rate is calculated for each tree within a forest.
# We will use the error rate from the last tree in the forest, which takes all previous trees into account and thus represents the error rate after the model stabilizes/converges

rf_all_1_err.rate <- tail(rf_all_1$err.rate,1)[1]
rf_all_2_err.rate <- tail(rf_all_2$err.rate,1)[1]
rf_all_3_err.rate <- tail(rf_all_3$err.rate,1)[1]


# # ############################################################################################################################################
# # ############################################################################################################################################

# Now conduct RF on subsets of the data to identify a group of loci that may be predictive of disease resistance.
# For each subset, we will use mtry=p since that is the optimal setting that we previously found.


all_initial_err_rate <- data.frame(cbind(1,  ncol(gt_data_corrected), rf_all_1_err.rate, rf_all_2_err.rate, rf_all_3_err.rate))

colnames(all_initial_err_rate) <- c('best_subset', 'nloci', 'rf1_oob', 'rf2_oob', 'rf3_oob')

perc <- c(0.01,
          0.02,
          0.03,
          0.04,
          0.05,
          0.1,
          0.2,
          0.3)

registerDoParallel(cores = ncpu)
cl <- makeCluster(ncpu)

initial_err_rate <- data.frame(
  best_subset = perc,
  nloci = NA,
  rf1_oob = NA,
  rf2_oob = NA,
  rf3_oob = NA
)

for (i in 1:length(perc)) {
  print('**********', quote = F)
  print(paste('Best ', perc[i], sep = ''), quote = F)

  names_best_1 <- rownames(importance_rf_all_1)[which(importance_rf_all_1$importance > quantile(importance_rf_all_1$importance, probs = 1 - perc[i]))]
  names_best_2 <- rownames(importance_rf_all_2)[which(importance_rf_all_2$importance > quantile(importance_rf_all_2$importance, probs = 1 - perc[i]))]
  names_best_3 <- rownames(importance_rf_all_3)[which(importance_rf_all_3$importance > quantile(importance_rf_all_3$importance, probs = 1 - perc[i]))]
  names_best_unique <- unique(c(names_best_1, names_best_2, names_best_3))

  names_best_unique %>%
    data.frame() %>%
    rename(loci = '.') %>%
    write.table(paste(output_path,'/best', perc[i], '.tsv', sep = ''), quote = F, sep = '\t', row.names = F)

  # Extract genotypes

  extract_gt <- colnames(gt_data_corrected)[colnames(gt_data_corrected) %in% names_best_unique]

  genotypes_perc <- gt_data_corrected %>%
    select(all_of(extract_gt))
  
  initial_err_rate$nloci[i] <- length(names_best_unique)
  
  for (j in 1:3) {
    rf <- parRF_classifier(pheno = pheno,
                               gt = genotypes_perc,
                               ntree = ntree,
                               mtry = length(names_best_unique),
                               ncpu = ncpu,
                               sample_size = sample_size)
    
    
    write_rds(rf, file = paste(output_path,
                               '/RDS/rf', 
                               j, 
                               '_', 
                               perc[i],
                               '.rds', 
                               sep = ''))
    
    initial_err_rate[i,j + 2] <- rf$err.rate[length(rf$err.rate)]
    colnames(initial_err_rate)[j + 2] <- paste('rf', j, '_oob', sep ='')
    
  }
  
}
stopCluster(cl)

all_initial_err_rate <- rbind(all_initial_err_rate, initial_err_rate)

all_initial_err_rate <- all_initial_err_rate %>%
  pivot_longer(-c(best_subset, nloci), names_to = 'run', values_to = 'err.rate')

all_initial_err_rate %>%
  fwrite(paste(output_path,'/all_initial_err_rate.tsv', sep = ''), sep = '\t', row.names = F)

all_initial_err_rate_avg <- all_initial_err_rate %>%
  mutate(err.rate = as.numeric(err.rate)) %>%
  group_by(best_subset) %>%
  summarise(mean.err.rate = mean(err.rate)) %>%
  arrange(desc(mean.err.rate))

all_initial_err_rate_avg %>%
  arrange(mean.err.rate) %>%
  ggplot(aes(y = reorder(best_subset, -mean.err.rate), x = mean.err.rate)) +
  xlab('Mean OOB error') +
  ylab('Subset of loci') +
  theme_minimal() +
  geom_bar(stat = 'identity')

# #################### Backward purging approach

purging_data <- all_initial_err_rate_avg %>%
  arrange(mean.err.rate) %>%
  head(2) %>%
  tail(1) %>%
  select(best_subset) %>%
  pull()

names_purging <- fread(paste(output_path, '/best', purging_data, '.tsv', sep = '')) %>%
  pull()

extract_gt <- colnames(gt_data_corrected)[colnames(gt_data_corrected) %in% names_purging]
genotypes_purging <- gt_data_corrected %>%
  select(all_of(extract_gt))

ncpu <- 10
registerDoParallel(cores = ncpu)
cl <- makeCluster(ncpu)

for (i in 1:3) {
  print('**********', quote = F)
  print(paste('rf purging ', i, sep = ''), quote = F)

  rf_purging <- parRF_classifier(pheno = pheno,
                                 gt = genotypes_purging,
                                 ntree = ntree,
                                 mtry = length(names_purging),
                                 ncpu = ncpu,
                                 sample_size = sample_size)

  write_rds(rf_purging,file = paste(output_path,"/RDS/rf_all_purging_", i, '.rds', sep = ''))
}

stopCluster(cl)

rf_purging_1 <- readRDS(paste(output_path,'/RDS/rf_all_purging_1.rds', sep = ''))
rf_purging_2 <- readRDS(paste(output_path,'/RDS/rf_all_purging_2.rds', sep = ''))
rf_purging_3 <- readRDS(paste(output_path,'/RDS/rf_all_purging_3.rds', sep = ''))

names_all_iterations <- list()
names_all_iterations[[length(names_purging)]] <- names_purging

err.rate_best <- data.frame(V1 = 1:length(names_purging),
                            V2 = 1:length(names_purging),
                            V3 = 1:length(names_purging))

rownames(err.rate_best) <- 1:length(names_purging)
err.rate_best[length(names_purging),] <- c(tail(rf_purging_1$err.rate,1)[1],
                                           tail(rf_purging_2$err.rate,1)[1],
                                           tail(rf_purging_3$err.rate,1)[1])

for (i in 1:(length(names_purging)-2)){  # RF cannot be conducted with 1 locus, which is why the loop is from 1:length(names_purging)-2
  print(i)
  imp_purging_1 <- data.frame(importance(rf_purging_1,type=1))
  imp_purging_2 <- data.frame(importance(rf_purging_2,type=1))
  imp_purging_3 <- data.frame(importance(rf_purging_3,type=1))

  colnames(imp_purging_1)<-"Mean_Decrease_Accuracy1"
  colnames(imp_purging_2)<-"Mean_Decrease_Accuracy2"
  colnames(imp_purging_3)<-"Mean_Decrease_Accuracy3"

  all_imp <- cbind(imp_purging_1,imp_purging_2,imp_purging_3)
  all_imp$average <- apply(all_imp[,1:3],1,mean)

  dont_keep <- which(all_imp[,'average']==min(all_imp[,'average']))
  
  if (length(dont_keep) == 1) {
    table_keep <- all_imp[-dont_keep,]
  } else {
    table_keep <- all_imp[-dont_keep[sample(x = dont_keep, size = 1)]]
  }

  names_keep <- rownames(table_keep)
  names_all_iterations[[length(names_purging)-i]] <- names_keep

  extract_gt <- colnames(gt_data_corrected)[colnames(gt_data_corrected) %in% names_keep]
  genotypes_purging <- gt_data_corrected %>%
    select(all_of(extract_gt))

  
  for (j in 1:3) {
    rf_purging <- parRF_classifier(pheno = pheno,
                                   gt = genotypes_purging,
                                   ntree = ntree,
                                   mtry = length(names_keep), 
                                   ncpu = ncpu,
                                   sample_size = sample_size)
    
    err.rate_best[length(names_purging)-i,j] <- tail(rf_purging$err.rate,1)[1]
  }
  
  #err.rate_best[length(names_purging)-i,] <- c(tail(rf_purging_1$err.rate,1),tail(rf_purging_2$err.rate,1),tail(rf_purging_3$err.rate,1))
}


err.rate_best$Average<-apply(err.rate_best,1,mean)
write.table(err.rate_best, paste(output_path,'/backward_purging.tsv', sep = ''), quote = F, sep = '\t')

# Now plot the backward purging results. Omit the row from one locus since RF cannot be conducted with just one locus
err.rate_best[-1,] %>%
  ggplot(aes(x = 2:nrow(err.rate_best), y = Average)) +
  geom_point(stat = 'identity') +
  ylab('mean OOB error') +
  xlab('Number of loci') +
  theme_minimal()
ggsave(paste(output_path,'/purging_results.png', sep = ''), 
       bg = 'white')
  
# Which group of loci explains the most variation?
best_pred <- which(err.rate_best$Average==max(err.rate_best$Average[-c(1)]))

# Export the names of the predictor loci
write.table(names_all_iterations[[best_pred]], 
          paste(output_path,'/best_predictior_loci.tsv', 
                sep = ''),
          row.names = F,
          col.names = F,
          quote = F)
