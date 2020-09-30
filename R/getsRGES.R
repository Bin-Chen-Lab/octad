#' @export
getsRGES <- function(RGES, cor, pert_dose, pert_time, diff, max_cor){
  
  sRGES <- RGES
  pert_time <- ifelse(pert_time < 24, "short", "long")
  pert_dose <- ifelse(pert_dose < 10, "low", "high")
  if (pert_time == "short" & pert_dose == "low"){
    sRGES <- sRGES + diff[4]
  }
  if (pert_dose ==  "low" & pert_time == "long"){
    sRGES <- sRGES + diff[2]
  }
  if (pert_dose ==  "high" & pert_time == "short"){
    sRGES <- sRGES + diff[1]
  }
  return(sRGES * cor/max_cor) 
}
