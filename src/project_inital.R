## 
## LOWRP Model Evaluation
## To support Draft PIR/EIS review
## 

## Code was compiled by Paul Julian
## contact info: pjulian.sccf.org

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

#Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);

#Paths
wd="C:/Julian_LaCie/_GitHub/LOWRP_ModelEval"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/GIS","/src/","/_documents/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

## Initalize ReadMe #one and done
## usethis::use_readme_rmd(open = rlang::is_interactive()); 