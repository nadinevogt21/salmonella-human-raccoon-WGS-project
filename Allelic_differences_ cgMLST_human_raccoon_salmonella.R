# code created by Adam Vogt. Meant for use with small datasets only, in this case, n=608 isolates. 

library(tidyr)
library(dplyr)
library(tidyverse)

ph <- read.csv("~/Dropbox/SALMONELLA ECOLI raccoon PhD/Chapter 5 Salmonella genomics chapter/manuscript files/cgMLST_dataset.csv")

#first bit of code here for comparing isolates between two different sources. 
#determine which serovars were identified in both raccoons and humans 

P1 <- ph %>% filter(sample_source=='Raccoon') %>% select(cgmlst_serovar)
P2 <-ph %>% filter(sample_source=='human')%>% select(cgmlst_serovar)

Pcommon <- intersect(  P1$cgmlst_serovar, P2$cgmlst_serovar )
rbind(
  P1$cgmlst_serovar %>% table %>% `[`(Pcommon),
  P2$cgmlst_serovar %>% table %>% `[`(Pcommon) ) %>%
  { rownames(.) <- c('Raccoon', 'human');  . }

sv <-c('Infantis', 'Newport', 'Typhimurium', 'Heidelberg', 'Enteritidis')

RH <- cgMLST_dataset %>% filter(sample_source %in% c('Raccoon','human') & cgmlst_serovar %in% sv)

sacvars <- names(ph) %>% tail(-5)

RaccoonHuman <- RH %>% group_by(cgmlst_serovar) %>% nest %>% mutate(s = map(data, ~ {
  m <- .[ .$sample_source == 'Raccoon', sacvars]
  n <- .[ .$sample_source == 'human', sacvars]
  #browser()
  s <- expand_grid(i=seq(nrow(m)), j=seq(nrow(n))) %>% rowwise %>%
    mutate(s= sum(m[i,] != n[j, ], na.rm=T)) %>% `$`('s') %>% range(na.rm=T)
}))

#code for determining allelic differences for a given serovar, for all sources. 
library(pacman)
p_load(tidyverse, ggplot2, directlabels, gt)

ph <- read.csv("~/Dropbox/SALMONELLA ECOLI raccoon PhD/Chapter 5 Salmonella genomics chapter/manuscript files/cgMLST_dataset.csv")

# the maximum number of alleles
diameters <- function(d) {
  ch <- select(d, STMMW_06041 : STMMW_02191) %>% as.matrix
  # %>% convhulln too few points
  if (nrow(ch) <= 1) NA else {
    map(2:nrow(ch), function(j)
      map_dbl(1:(j-1), function(i) {
        sum( ch[i, ] != ch[j, ] )
      })) %>% do.call(c, .)
  }
}

sacphD <- ph %>% group_by(cgmlst_serovar) %>% nest() %>%
  summarize(
    nsamp = data[[1]] %>% nrow(),
    ds = diameters(data[[1]]) %>% list,
    .groups="keep")

sacphD %>%
  summarize(maxDiff = map_dbl(ds, max),
            minDiff=map_dbl(ds, min),
            npairs = map_dbl(ds, length),
            nsamp= nsamp) %>%
            { .[order(.$maxDiff), ] } %>%
  gt


