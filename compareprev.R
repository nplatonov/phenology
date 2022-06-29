plutil::ursula(3)
season <- "sid2012"
#prev <- read_envi(file.path("../v43nt/phenology",season))["v"]
prev <- ursa_read(envi_list(path=c("../../v43nt/phenology"
                                ,"D:/RAS/2013/appear/versions/v43nt")
                         ,pattern=season,full.names=TRUE))["v"]
prev <- prev[prev<600]
a <- ursa_read(file.path(c("output",".")[2],season))
onset <- a[a<600]
suppl <- a[a>600]
as.table(suppl)
#display(list(a,prev,onset,suppl))
p <- c(prev=prev,dev=onset)
d <- c(diff=p["dev"]-p["prev"])
m <- !is.na(d)
p2 <- list(onset=p[m],difference=d)
p2
display(p2,height=1000,stretch="eq",ramp=FALSE,las=1)
#c(p,p[2]-p[1])
