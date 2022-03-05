prac_names = list.files(path = "Practicals/", pattern = "*.Rmd")

for(i in 1: length(prac_names))
  rmarkdown::render(paste("Practicals/", prac_names[i], sep = ""),
                  encoding="UTF-8") 

