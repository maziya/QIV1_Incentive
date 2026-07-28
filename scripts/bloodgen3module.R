devtools::install_github("Drinchai/BloodGen3Module")
mod <- BloodGen3Module:::Module_listGen3   # triple-colon: it's internal/unexported
ann <- BloodGen3Module:::Gen3_ann

library(dplyr)
gmt_list <- split(mod$Gene, mod$Module)
gmt_lines <- mapply(function(mid, genes) {
  fn <- ann$Function[ann$Module == mid][1]   
  desc <- paste0(fn, "(", mid, ")")
  paste(c(mid, desc, genes), collapse = "\t")
}, names(gmt_list), gmt_list)
writeLines(gmt_lines, "BloodGen3_modules_v2.gmt")
