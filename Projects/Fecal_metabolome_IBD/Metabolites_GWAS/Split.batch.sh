# split traits file into multiple batches, keep header
# adjust parameter "l" to define line number in each batch (the number of traits)

awk -v l=10 '(NR==1){header=$0;next}
                (NR%l==2) {
                   close(file); 
                   file=sprintf("%s.%0.5d.txt",FILENAME,++c)
                   sub(/csv[.]/,"",file)
                   print header > file
                }
                {print > file}' trait.file
