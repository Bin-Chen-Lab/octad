#' @export
cmap_score_ultimate=function (sig_up, sig_down, drug_signature) 
{
    num_genes <- length(drug_signature)
    ks_up <- 0
    ks_down <- 0
    connectivity_score <- 0
	
	
	drug_signature=rank(drug_signature)
	up_tags_rank=drug_signature[as.vector(sig_up)]
	down_tags_rank=drug_signature[as.vector(sig_down)]

up_tags_position=sort(up_tags_rank)
down_tags_position=sort(down_tags_rank)

    num_tags_up <- length(up_tags_position)
    num_tags_down <- length(down_tags_position)
    if (num_tags_up > 1) {
        a_up <- 0
        b_up <- 0
        a_up <- max(sapply(1:num_tags_up, function(j) {
            j/num_tags_up - up_tags_position[j]/num_genes
        }))
        b_up <- max(sapply(1:num_tags_up, function(j) {
            up_tags_position[j]/num_genes - (j - 1)/num_tags_up
        }))
        if (a_up > b_up) {
            ks_up <- a_up
        }
        else {
            ks_up <- -b_up
        }
    }
    else {
        ks_up <- 0
    }
    if (num_tags_down > 1) {
        a_down <- 0
        b_down <- 0
        a_down <- max(sapply(1:num_tags_down, function(j) {
            j/num_tags_down - down_tags_position[j]/num_genes
        }))
        b_down <- max(sapply(1:num_tags_down, function(j) {
            down_tags_position[j]/num_genes - (j - 1)/num_tags_down
        }))
        if (a_down > b_down) {
            ks_down <- a_down
        }
        else {
            ks_down <- -b_down
        }
    }
    else {
        ks_down <- 0
    }
    if (ks_up == 0 & ks_down != 0) {
        connectivity_score <- -ks_down
    }
    else if (ks_up != 0 & ks_down == 0) {
        connectivity_score <- ks_up
    }
    else if (sum(sign(c(ks_down, ks_up))) == 0) {
        connectivity_score <- ks_up - ks_down
    }
    else {
        connectivity_score <- ks_up - ks_down
    }
    return(connectivity_score)
}
