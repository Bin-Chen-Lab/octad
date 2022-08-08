

data=data_input
pheno_data=pheno_input
ag1=
temp_data_frame=as.data.frame(t(mapply(compute_bsabs_score,ag1=result_table$antigen_1,ag2=result_table$antigen_2)))
colnames(temp_data_frame)=c("antigen_1","antigen_2","distance","spread","angle_cos","pair_score","case_greater")
temp_data_frame[3:6]=apply(temp_data_frame[3:7],2,function(x)as.numeric(as.character(x)))



compute_bsabs(antigene_1=result_table$antigen_1,antigene_2=result_table$antigen_2,data_input=train,pheno_input=pheno_train)
compute_bsabs_single(result_table$antigen_1[1],result_table$antigen_2[1],data=train,pheno_data=pheno_train)


small_res=compute_bsabs(antigene_1=result_table$antigen_1,antigene_2=result_table$antigen_2,data_input=train,pheno_input=pheno_train)

.

input_table=small_res
input_table=na.omit(input_table)
input_table$significant='NO'
input_table[input_table$p.adj<0.01&input_table$pair_score>quantile(input_table$pair_score,.99)&input_table$case_greater=='TRUE_TRUE','significant']='Enriched in case'
input_table[input_table$p.adj<0.01&input_table$pair_score>quantile(input_table$pair_score,.99)&input_table$case_greater=='FALSE_FALSE','significant']='Enriched in control'
table(input_table$significant)

#label DE genes
input_table$pair_label <- NA
input_table$pair_label[input_table$significant == 'Enriched in case'] = 
  apply(input_table[input_table$significant == 'Enriched in case',], 1, function(x) paste(x[1:2], collapse = ":"))
  


#fancy volcano plot
ggplot(data=input_table, aes(y=pair_score ,x=-log10(p.adj), col=significant, label=pair_label)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps =30) +
  scale_color_manual(values=c("red", "blue", "black")) 

  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
