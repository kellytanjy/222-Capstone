##read in data
library(readr)
library(tidyr)
library(dplyr)
library(reshape2)
Broad202 = read_delim("222_data/Broad202.txt","\t", escape_double = FALSE, trim_ws = TRUE)
LBL202 = read_delim("222_data/LBL202.txt","\t", escape_double = FALSE, trim_ws = TRUE)
UNC202 = read_delim("222_data/UNC202.txt","\t", escape_double = FALSE, trim_ws = TRUE)
colnames(UNC202)[1]="Gene"

firstrow=read_delim("222_data/unifiedScaled.txt", 
                    "\t", escape_double = FALSE, col_names = FALSE,
                    trim_ws = TRUE, skip = 0,n_max=1)
unifiedScaled = read_delim("222_data/unifiedScaled.txt",
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE, skip = 1) %>% 
  `colnames<-`(value=c('Gene',firstrow))
rm(list="firstrow")

merge_genes=intersect(LBL202$Gene,UNC202$Gene) %>% intersect(Broad202$Gene)
merge_ids=intersect(colnames(LBL202)[-1],colnames(UNC202)[-1]) %>% 
  intersect(colnames(Broad202)[-1])

#reshape the data
unc= UNC202 %>%
  filter(Gene %in% merge_genes) %>%
  select(Gene,match(merge_ids,colnames(UNC202))) %>%
  melt(id.vars="Gene")
broad=Broad202 %>%
  filter(Gene %in% merge_genes) %>%
  select(Gene,match(merge_ids,colnames(Broad202))) %>%
  melt(id.vars='Gene')
lbl=LBL202 %>%
  filter(Gene %in% merge_genes) %>%
  select(Gene,match(merge_ids,colnames(LBL202))) %>%
  melt(id.vars='Gene')

total=unc %>%
  merge(broad,by=c('Gene','variable')) %>%
  merge(lbl,by=c('Gene','variable')) %>%
  `colnames<-`(value=c('Gene','Case','unc','broad','lbl'))
totaltry=total %>% 
  mutate(
    unc=scale(unc),
    broad=scale(broad),
    lbl=scale(lbl))

#Factor Analysis
fit=factanal(totaltry[,3:5],1)
#print(fit)
t=as.matrix(totaltry[,3:5]) %*% fit$loadings
fa=data.frame(Gene=totaltry$Gene,Case=totaltry$Case,Measure=t)
reshapeback=dcast(fa,Gene~Case)
#firstrow=read_delim("222_data/unifiedScaled.txt", 
#                    "\t", escape_double = FALSE, col_names = FALSE, 
#                    trim_ws = TRUE, skip = 0,n_max=1)
#unified=read_delim("222_data/unifiedScaledFiltered.txt",
#                    "\t", escape_double = FALSE, col_names = FALSE, 
#                    trim_ws = TRUE, skip = 1) %>%
#  `colnames<-`(value=c('Hugo_symbol',firstrow)) %>%
#  t() %>% 
#  as_data_frame()
#  
#colnames(filtered)=filtered[1,]
#filtered=filtered[-1,] %>%
#  sapply(as.numeric) %>%
#  as_data_frame()




