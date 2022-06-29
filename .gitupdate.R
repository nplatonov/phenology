desc <- commandArgs(TRUE)
if (!length(desc)) {
   desc <- "ongoing"
}
system("git add -A")
if (!system(paste0("git commit -m",dQuote(desc))))
   system("git push")
