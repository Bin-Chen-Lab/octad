

compute_bsabs=function(antigene_1,antigene_2=NULL,data_input,pheno_input){
#data=data_input
#colnames(data)=gsub('-','_',colnames(data))
#antigene_1=gsub('-','_',antigene_1)
#antigene_2=gsub('-','_',antigene_2)
#pheno_data=pheno_input
#ag1=antigene_1
#ag2=antigene_2
#temp_data_frame=as.data.frame(t(mapply(compute_bsabs_single,ag1=antigene_1,ag2=antigene_2)))
#colnames(temp_data_frame)=c("antigen_1","antigen_2","distance","spread","angle_cos","pair_score","case_greater")
#temp_data_frame[3:6]=apply(temp_data_frame[3:7],2,function(x)as.numeric(as.character(x)))
#return(temp_data_frame)
#plug, will replace a bit later with faster code above:
if(is.null(antigene_2)){
gene_vector=antigene_1
all_pairs=t(combn(gene_vector,2))
result_table=data.frame(antigen_1=all_pairs[,1],
                        antigen_2=all_pairs[,2],
                        distance=0,spread=0,angle_cos=0,
                        pair_score=0,p.value=0,p.adj=0,case_greater=FALSE)
}else{
gene_vector=unique(c(antigene_1,antigene_2))
result_table=data.frame(antigen_1=antigene_1,
                        antigen_2=antigene_2,
                        distance=0,spread=0,angle_cos=0,
                        pair_score=0,p.value=0,p.adj=0,case_greater=FALSE)
}
pb <- txtProgressBar(min = 1, max = nrow(result_table), style = 3)
counter=0


#replace - with _ in data input & antigenes
result_table$antigen_1=gsub('-','_',result_table$antigen_1)
result_table$antigen_2=gsub('-','_',result_table$antigen_2)
colnames(data_input)=gsub('-','_',colnames(data_input))
#check if antigenes exist in data_input
if(any(gene_vector%in%colnames(data_input)==FALSE)){
warning(paste('These antigenes do not exist in the expression matrix:','\n',
paste0(gene_vector[!gene_vector%in%colnames(data_input)],collapse=' ')),'\nomiting them')
result_table=subset(result_table,antigen_1%in%colnames(data_input)&antigen_2%in%colnames(data_input))
}


for (x in 1:nrow(result_table)){
  i=result_table[x,1]
  j=result_table[x,2]
    counter=counter+1
	#cat(paste0(round(counter/nrow(result_table)*100,3),'%'),'\n') #counter
	setTxtProgressBar(pb, counter) 
	
    if(i!=j){
##########define centers and dist
distance=sqrt(sum((pam(data_input[which(pheno_train=='case'),c(j,i)], 1)$medoids-
                       pam(data_input[which(pheno_train=='control'),c(j,i)], 1)$medoids)^2)) #distance between clusters
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'distance']=distance

#silhouette
temp_dist=dist(data_input[,c(j,i)])
silh=(mean(silhouette(as.numeric(pheno_train),temp_dist)[,3]))+1
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spread']=silh

#angle coefficient
table_temp=rbind(pam(data_input[which(pheno_train=='case'),c(j,i)], 1)$medoids,
                 pam(data_input[which(pheno_train=='control'),c(j,i)], 1)$medoids)
#Boolean logic to check the position of clusters
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'case_greater']=paste0(table_temp[1,]>table_temp[2,],collapse ='_')
 
#angle to between 45 and metoids
m1=lm(table_temp[,i]~table_temp[,j])$coefficients[2]
cos_to_45=1/(1+(abs((1-m1)/(1+m1*1)))^2)
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'angle_cos']=cos_to_45
#score
pair_score=cos_to_45*silh*distance
result_table[result_table$antigen_1==i&result_table$antigen_2==j,'pair_score']=pair_score
#pvalue
  table_temp=data_input[,c(j,i)]
  table_temp$pheno=pheno_input
 pvalue=wilcox.test(as.formula(paste(i,'+',j,'~pheno',sep='')),data=table_temp)$p.value
 result_table[result_table$antigen_1==i&result_table$antigen_2==j,'p.value']=pvalue
}
}
result_table$p.adj=p.adjust(result_table$p.value,method='fdr')
result_table=na.omit(result_table)
return(result_table)
}




#visualize
plot_bsabs=function(input_table,label=c('all','case','control'),pval_cut_off=0.01,pair_score_cut_off=NULL){
if(is.null(pair_score_cut_off)){
pair_score_cut_off=quantile(input_table$pair_score,.99)
}

input_table=na.omit(input_table)
input_table$significant='NO'
input_table[input_table$p.adj<pval_cut_off&input_table$pair_score>pair_score_cut_off&input_table$case_greater=='TRUE_TRUE','significant']='Enriched in case'
input_table[input_table$p.adj<pval_cut_off&input_table$pair_score>pair_score_cut_off&input_table$case_greater=='FALSE_FALSE','significant']='Enriched in control'
#table(input_table$significant)

#label DE genes
input_table$pair_label <- NA
input_table$pair_label[input_table$significant == 'Enriched in case'] = 
  apply(input_table[input_table$significant == 'Enriched in case',], 1, function(x) paste(x[1:2], collapse = ":"))
  


#fancy volcano plot
plot=ggplot(data=input_table, aes(y=pair_score ,x=-log10(p.adj), col=significant, label=pair_label)) +
  geom_point() + 
  theme_minimal() +
  ggrepel::geom_text_repel(max.overlaps =30) +
  scale_color_manual(values=c("red", "blue", "black")) 
print(plot)
}












compute_bsabs_single=function(ag1,ag2,data=data_input,pheno_data=pheno_input){
  distance=sqrt(sum((pam(data[which(pheno_data=='case'),c(ag2,ag1)], 1)$medoids-
                       pam(data[which(pheno_data=='control'),c(ag2,ag1)], 1)$medoids)^2)) #distance between clusters
  #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'distance']=distance
  #silhouette
  temp_dist=dist(data[,c(ag2,ag1)])
  silh=(mean(silhouette(as.numeric(pheno_data),temp_dist)[,3]))+1
  #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'spread']=silh
  #angle coefficient
  table_temp=rbind(pam(data[which(pheno_data=='case'),c(ag2,ag1)], 1)$medoids,
                   pam(data[which(pheno_data=='control'),c(ag2,ag1)], 1)$medoids)
  #Boolean logic to check the position of clusters
  boolean_vector=paste0(table_temp[1,]>table_temp[2,],collapse ='_')

  m1=lm(table_temp[,ag1]~table_temp[,ag2])$coefficients[2]
  cos_to_45=1/(1+(abs((1-m1)/(1+m1*1)))^2)
  #result_table[result_table$antigen_1==i&result_table$antigen_2==j,'angle_cos']=cos_to_45
  #score
  pair_score=cos_to_45*silh*distance
  table_temp=data[,c(ag2,ag1)]
  table_temp$pheno=pheno_data
  pvalue=wilcox.test(as.formula(paste(ag1,'+',ag2,'~pheno',sep='')),data=table_temp)$p.value
  return(c(i,j,distance,silh,cos_to_45,pair_score,pvalue,boolean_vector))
}
