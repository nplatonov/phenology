plutil::ursula()
'check_flags' <- function() {
   nature <- apply(expand.grid(pheno=c("sid","sia"),year=seq(1988,2024)),1
              ,function(x) paste(x,collapse=""))
   list1 <- envi_list(pattern="si[ad]\\d{4}",path=c("./output","../v43nt/phenology")[1]
                     ,full.names=TRUE)
   if (!length(list1)) {
      list1 <- dir(pattern="si[ad]\\d{4}\\.tif",path="./output",full.names=TRUE)
   }
   list1 <- list1[na.omit(match(nature,gsub(".*(si[ad]\\d{4}).*","\\1"
                 ,basename(list1))))] |> tail(24)
   str(list1)
   q()
   res <- sapply(list1,function(fname) {
      a <- ursa_read(fname)
      if (length(a)>1)
         a <- a["kind"]
      a
   }) |> do.call(c,args=_)
   print(res[res<600])
   res[res<600] <- 700
   if (F) {
      res[res==830] <- 880
      res[res==840] <- 870
   }
   if (F) {
      res[res==794] <- 880 ## 794->820->870
      res[res==792] <- 870 ## 792->810->870
   }
  # display(res,scale=3)
   ta <- as.table(res)
   da <- array(0,dim=c(length(res),length(ta)),dimnames=list(names(res),names(ta)))
   for (i in seq(res)) {
      ta <- as.table(res[i])
      da[i,match(names(ta),colnames(da))] <- ta
   }
   print(da)
   #a
   #display(res)
   #as.table(res)
   0L
}
'inconsistent_pairs' <- function() {
   p1 <- ursa_read("sid2021")
   p2 <- ursa_read("sia2021")
   m <- p2==792 | p2==794
   p <- c(p1=p1[m],p2=p2[m])
   print(p)
   if (T)
      p |> as.data.frame() |> spatialize(crs=ursa_crs(p1)) |>
         spatial_write("inconsistent.geojson")
   0L
}
'doublecheck' <- function() {
   p1 <- ursa_read("output.v7/sid2021.tif")
   p2 <- ursa_read("output/sid2021.tif")
   print(as.table(p1["flag"]))
   print(as.table(p2["flag"]))
   0L
}
'nugget' <- function() {
   p1 <- c(z=ursa_read("output.v5/sid2021.tif")["value"])
   print(p1)
   b <- as.data.frame(p1)
   sp::coordinates(b) <- ~x+y
   spatial_crs(b) <- spatial_crs(p1)
   v1 <- gstat::variogram(z~1,data=b)
   str(v1)
   v2 <- gstat::variogram(z~x+y,data=b)
   str(v2)
   0L
}
invisible({
   ##~ check_flags()
   ##~ inconsistent_pairs()
   doublecheck()
   ##~ nugget()
})
