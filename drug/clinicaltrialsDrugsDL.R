#place into drug/core_function when done

queryClinicalTrials <- function(dz){
  library(curl)
  library(stringr)
  #library(XML)
  #dz = 'Oligodendroglioma'
  url <- paste0("https://clinicaltrials.gov/ct2/results/download_fields?cond=",
                dz,"&down_fmt=csv")
  tmp <- "~/Documents/GitHub/OCTAD/test.csv"
  curl_download(url, tmp)
  csvfile <- read.csv(tmp)
  #gets only interventions that are drugs
  Drugs <- lapply(as.character(csvfile$Interventions), 
                      function(x){
                        unlist(strsplit(x,split='|',fixed=T))
                      })
  
  DrugsOnly <- lapply(as.character(csvfile$Interventions), 
                  function(x){
                    grep('Drug',unlist(strsplit(x,split='|',fixed=T)))
                    })
  temp <- c()
  for (i in 1:length(Drugs)) {
    #x <- Drugs[[i]][DrugsOnly[[i]]]
    temp <- c(temp,Drugs[[i]][DrugsOnly[[i]]])
  }
  temp <- tolower(temp)
  #Drug <- unlist(strsplit(temp,split='drug: ',fixed=T))
  #cleaning names as much as possible
  Drug <- temp %>% 
    str_replace("drug: ", "") %>% 
    str_replace("hydrochloride","") %>% 
    str_replace("hci","") %>% 
    str_replace("hcl","") %>% 
    trimws()
    
  Drug <- data.frame(table(Drug))
}
