#' @export
#' @import ggplot2
#' @import ggrepel


#visualize
plot_bsabs=function(input_table,label=c('all','case','control','none'),pval_cut_off=0.01,pair_score_cut_off=NULL){
if(is.null(pair_score_cut_off)){
pair_score_cut_off=quantile(input_table$pair_score,.99)
}

input_table=na.omit(input_table)
input_table$significant='NO'
input_table[input_table$p.adj<pval_cut_off&input_table$pair_score>pair_score_cut_off&input_table$case_greater=='TRUE_TRUE','significant']='Enriched in case'
input_table[input_table$p.adj<pval_cut_off&input_table$pair_score>pair_score_cut_off&input_table$case_greater=='FALSE_FALSE','significant']='Enriched in control'
#table(input_table$significant)

#label pairs

input_table$pair_label <- NA
if(label=='case'){
input_table$pair_label[input_table$significant == 'Enriched in case'] = 
  apply(input_table[input_table$significant == 'Enriched in case',], 1, function(x) paste(x[1:2], collapse = ":"))
}else if(label=='control'){
input_table$pair_label[input_table$significant == 'Enriched in control'] = 
  apply(input_table[input_table$significant == 'Enriched in control',], 1, function(x) paste(x[1:2], collapse = ":"))
}else if(label=='all'){
input_table$pair_label[input_table$significant != 'NO'] = 
  apply(input_table[input_table$significant != 'NO',], 1, function(x) paste(x[1:2], collapse = ":"))
}else if(label=='none'){}


#fancy volcano plot
plot=ggplot(data=input_table, aes(y=pair_score ,x=-log10(p.adj), col=significant, label=pair_label)) +
  geom_point() + 
  theme_minimal() +
  ggrepel::geom_text_repel(max.overlaps =30) +
  scale_color_manual(values=c("red", "blue", "black")) 
print(plot)
}
