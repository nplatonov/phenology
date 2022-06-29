plutil::ursula()
if (F) {
   ice <- ursa_read("output.v6/season_ice.tif")
   b <- local_stat(ice,cover=0.51)[c("mean","sd","slope","slopeS")]
   print(b)
   display(b,ramp=FALSE)
   q()
}
if (F) {
   ice <- ursa_read("output.v6/season_ice.tif")
   print(ice)
   q()
   iucn <- ursa:::.fasterize(file.path("D:/RAS/2022/shiny/_archive_/demography"
                                     ,"IUCN_pb_subpopulations.shp.zip"))[["POP"]]
   print(iucn)
   print(series(ice))
   b <- aggregate(ice,iucn,mean,na.rm=TRUE,table=TRUE)
   print(b)
   #display_homo(series(ice))
   q()
}
if (F)
   ursa_write(ice["2016"],"diff.tif")
ice <- ice["2016"]
ice
lice <- ice[ice<0]
less <- as.data.frame(lice)
less
xy <- less[sample(nrow(less),1),]
xy
sia <- ursa_read("output/sia2015")#["date"]
sid <- ursa_read("output/sid2016")
value_xy(sia,x=xy[,1],y=xy[,2])
value_xy(sid,x=xy[,1],y=xy[,2])
sid["date"]-sia["date"]
#print(ice>=365)
ursa_write(ice["2016"],"diff.tif")

