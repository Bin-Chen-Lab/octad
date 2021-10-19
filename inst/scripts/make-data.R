##  Download files
download.file(url='https://chenlab-data-public.s3.amazonaws.com/octad/octad.db/CCLE.log2.read.count.matrix.rda',
              destfile='CCLE.log2.read.count.matrix.rda')

load('CCLE.log2.read.count.matrix.rda')
save(CCLE.log2.read.count.matrix, 
     file="CCLE.log2.read.count.matrix.rda")