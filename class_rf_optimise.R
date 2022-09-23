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

### path to genotypes

genotypes_path <- 'gt_obesity_panel.tsv'

print('loading genotypes', quote = F)
gt_data <- fread(genotypes_path, 
                 header = T,
                           nThread = ncpu) %>%
  select(-IID)

### path to covariates

covars_path <- 'covars.tsv'
covars <- read.table(covars_path, 
                     header = T) %>%
  mutate(sex = as.factor(sex),
         neutering = as.factor(neutering))

### path to phenotype 

phenotype_path <- 'bcs_pheno.tsv'

pheno <- read.table(phenotype_path, header = T) %>%
  pull() %>%
  as.factor()

### getting residuals of genotypes with linear regression

print('computing residuals', quote = F)
gt_data_corrected <- gt_data %>%
  map(~lm(.x ~ covars$sex + covars$age + covars$neutering)$residuals) %>%
  as.data.table()

rm(gt_data)

colnames(gt_data_corrected) <- gsub('^X', '', colnames(gt_data_corrected))

gt_data_corrected %>% 
  fwrite(paste(output_path,'/gt_data_corrected.tsv', sep = ''))

##### hyperparameters ######

## set number of trees

ntree <- c(seq(100, 5000, by = 100),
           seq(5000, 20000, by = 1000),
           seq(20000, 100000, by = 10000),
           seq(100000, 1000000, by = 10000))

## set number of loci per tree

p <- length(1:ncol(gt_data_corrected))
mtry <- round(
  c(sqrt(p), 2*sqrt(p), 0.1*(p), 0.2*(p), p/3, p)
)

## set sampling size

pheno_sample <- c(11,11)

results_optimization <- matrix(data=NA , nrow = 0, ncol = 5)

cl <- makePSOCKcluster(ncpu)
registerDoParallel(cl)

for (i in ntree) {
  print('********', quote = F)
  print(paste(i, ' trees', sep = ''), quote = F)
  print('********', quote = F)
  for (j in mtry) {
    print(paste(j, ' loci', sep = ''), quote = F)
    rf_all_1 <- parRF_classifier(pheno = pheno,
                      gt = gt_data_corrected,
                      ntree = i,
                      mtry = j,
                      ncpu = ncpu,
                      sample_size = pheno_sample)

    rf_all_2 <- parRF_classifier(pheno = pheno,
                      gt = gt_data_corrected,
                      ntree = i,
                      mtry = j,
                      ncpu = ncpu,
                      sample_size = pheno_sample)

    results_optimization <- rbind(results_optimization, 
                                  c(cor(importance(rf_all_1,type=1), importance(rf_all_2,type=1))[1], 
                                    i,
                                    j,
                                    tail(rf_all_1$err.rate,1)[1],
                                    tail(rf_all_2$err.rate,1)[1]))
  }
}
stopCluster(cl)

results_optimization<-as.data.frame(results_optimization)
colnames(results_optimization)<-c("cor", "ntree", "mtry", "OOB_ER_1", "OOB_ER_2")

results_optimization$mean_OOB <- rowMeans(results_optimization[4:5])

results_optimization %>% 
  fwrite(paste(output_path,'/results_optimization.tsv', sep = ''))

results_optimization %>%
  ggplot(aes(x = ntree, y = mean_OOB, col = factor(mtry))) +
  geom_line() +
  ylab('mean OOB Error') +
  xlab('Number of trees') +
  guides(col = guide_legend(title = 'Number of loci')) +
  theme_minimal()
ggsave(paste(output_path,'/results_optimization.png', sep = ''), 
       bg = 'white')

results_optimization %>%
  ggplot(aes(x = ntree, y = cor, col = factor(mtry))) +
  geom_line() +
  ylab(expression(
    paste("Pearson's ", R^{2})
  )) +
  xlab('Number of trees') +
  guides(col = guide_legend(title = 'Number of loci')) +
  theme_minimal()
ggsave(paste(output_path,'/cor_results_optimization.png', sep = ''), 
       bg = 'white')

