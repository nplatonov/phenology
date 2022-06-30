if (requireNamespace("plutil")) {
   plutil::ursula(3)
} else {
   require(ursa)
   stopifnot(packageVersion("ursa")>="3.9.7")
}
'seasonality' <- function(season) {
   if (missing(season))
      season <- "any2022"
   year <- as.integer(gsub("\\D","",season))
   pheno <- gsub("\\d","",season)
   transition <- if (pheno=="sid") c("max","mlt","min")
            else if (pheno=="sia") c("min","frz","max")
            else c("max","mlt","min","frz")
   shift <- as.integer(pheno=="sia")
   period <- list(mlt=list(from=as.Date(paste0(year,"-02-15"))
                          ,to=as.Date(paste0(year,"-09-30")))
                 ,frz=list(from=as.Date(paste0(year,"-08-15"))
                          ,to=as.Date(paste0(year+1L,"-03-31"))
                          )
                 ,min=list(from=as.Date(paste0(year,"-07-15"))
                          ,to=as.Date(paste0(year,"-11-15"))
                          )
                 ,max=list(from=as.Date(paste0(year+shift,"-01-15"))
                          ,to=as.Date(paste0(year+shift,"-05-15"))
                          )
                  )[transition]
   period
}
'iceconc_download' <- function(iceconc="./iceconc",season="sia2021"
                     ,nsidc=c(nt="0051",nrt="0081",bs="0079"),verbose=F) {
   wd <- setwd(iceconc);on.exit(setwd(wd))
   nsidc <- match.arg(nsidc)
   remote <- file.path("https://n5eil01u.ecs.nsidc.org/PM"
                      ,switch(nsidc,'0079'="NSIDC-0079.003"
                                   ,'0051'="NSIDC-0051.001"
                                   ,'0081'="NSIDC-0081.001"))
   isNRT <- grepl("0081",nsidc)
   year <- as.integer(gsub("\\D","",season))
   pheno <- gsub("\\d","",season)
   period <- seasonality(season)
   transition <- names(period)
   delta <- 3L
   d3 <- do.call(c,do.call(c,period)) |> range()
   d3 <- seq(d3[1]-delta,d3[2]+delta,by="1 days")
   isMelt <- pheno %in% "sid"
   isFreeze <- pheno %in% "sia"
   print(series(d3))
  # conc3 <- file.path(iceconc,paste0("conc2",season))
  # if (envi_exists(conc3))
  #    return(10L)
  # conc1 <- file.path(iceconc,paste0("tmp1"))
  # conc2 <- file.path(iceconc,paste0("tmp2"))
   if (!dir.exists(iceconc))
      dir.create(iceconc)
   fileout <- paste0("conc1",season,".tif")
   print(range(d3))
   d3name <- format(d3,"%y%m%d")
   msk <- c(msk=ursa("ref/dist2land")>=0)
   if (F) {
      belt <- as.data.frame(ursa(1L))
      xy <- sf::sf_project(from=spatial_crs(msk),to=spatial_crs(4326)
                          ,pts=cbind(belt$x,belt$y))
      belt <- allocate(belt[xy[,2]>86.5,])*100
      print(range(d3),quote=FALSE)
   }
   if (ursa_exists(fileout))
      conc <- ursa_read(fileout)
   else
      conc <- ursa(NA,bandname=format(d3,"%y%m%d"))
   dst <- tempfile()
   if (verbose)
      cat("read metadata")
   for (i in seq(conc) |> sample(1)) {
      d <- d3[i]
      print(d)
      if (!band_blank(conc[i]))
         next
      if (verbose)
         cat(".")
      src <- file.path(remote,format(d,"%Y.%m.%d")
                     # ,format(d,"nt_%Y%m%d_f18_nrt_n.bin.xml")
                      )
      ret <- try(download.file(src,dst,method="curl",quiet=TRUE))
      if (inherits(ret,"try-error"))
         return(ursa())
      content <- xml2::read_html(dst)
      lst <- try(rvest::html_table(content,fill=TRUE)[[1]])
      if (inherits(lst,"try-error"))
         return(ursa())
      fname <- grep("_n\\.bin$",lst$Name,value=TRUE)
      ret <- try(download.file(file.path(src,fname),dst,method="curl",quiet=TRUE))
      if (inherits(lst,"try-error"))
         return(ursa())
      fs <- file.size(dst)
      if (fs!=136492)
         return(ursa())
      a <- readBin(dst,integer(),size=1,n=fs,signed=FALSE)
      a <- tail(a,-300)
      session_grid(msk)
      b <- as_ursa(a)
      f <- b[b>250]
      b <- b[b<=250]/2.5
      b[f==251] <- 100
      conc[i] <- b
      ursa_write(conc,fileout)
   }
   0L
}
'iceconc_filter' <- function() {
   iceseason <- unique(sort(format(d3,"%Y%m")))
   suff <- ifelse(method=="bs","bt",method)
   absent <- NULL
   q()
   a <- create_envi(msk,conc1,bands=length(d3),bandname=d3name
                   ,datatype=4L,bg=-99,interleave="bil")
   future <- TRUE
   for (s in season)
   {
      d3s <- seq(as.Date(paste0(s,"01"),"%Y%m%d"),len=2,by="1 month")
      d3s <- format(seq(d3s[1],d3s[2]-1L,by="1 days"),"%y%m%d")
      zipname <- list1[mygrep(substr(s,1,4),basename(list1))]
      if (!length(zipname))
         next
      if (future)
         future <- FALSE
      zippattern <- paste0(substr(s,3,6),"\\d{2}",suff)
      d3nt <- idrisi.ziplist(zipname,zippattern)
      d3zip <- mygsub("\\D","",d3nt)
      ind <- match(d3name,d3zip)
      ind1 <- which(!is.na(ind))
      ind2 <- na.omit(ind)
      ind3 <- which(is.na(match(d3s,d3zip)))
      ind4 <- which(!is.na(match(d3name,d3s[ind3])))
      if (length(ind2))
      {
         d <- read.idrisi(zipname,d3nt[ind2])
         d[d>100] <- NA
         d[is.na(d)] <- belt
         d <- focal_special(d,"gaus",size=3,sigma=0.75,cover=0.35,fillNA=TRUE
                       ,verbose=TRUE)
         d[msk<1] <- NA
         a[ind1] <- d
         print(a[ind1])
      }
      if (length(ind4))
      {
         absent <- c(absent,ind4)
         a[ind4] <- ursa_new(value=NA,nband=length(ind4))
      }
   }
   close(a)
   print(absent)
   if (future) {
      envi_remove(conc1)
      return(20L)
   }
   a <- open_envi(conc1)
   b <- create_envi(a,conc2)
   for (i in chunk_line(a))
   {
      d <- temporal_interpolate(a[,i],win=7)
      ##~ d[msk[,i]<1] <- NA
      b[,i] <- d
   }
   close(a)
   close(b)
   a <- open_envi(conc2)
   isEmpty <- all(band_blank(head(a))) | all(band_blank(tail(a)))
   if (isEmpty)
   {
      cat("omit empty images...")
      e <- NULL
      for (i in chunk_band(a))
         e <- c(e,i[which(band_blank(a[i]))])
      b <- create_envi(a,conc3,bandname=bandname(a)[-e])
      for (i in chunk_band(b))
         b[i] <- a[i]
      close(a,b)
      cat(" done!\n")
   }
   else {
      close(a)
      envi_copy(conc2,conc3)
   }
   envi_remove(c(conc1,conc2))
   0L
}
'phenology' <- function(season="sid2021",output=c(".","output")
                       ,iceconc="../v43nt/conc"
                       ,th=15,slen=30,doy=TRUE,overwrite=FALSE,ind=NULL) {
  # require(Iso)
   devel <- length(ind)>0
   if ((!devel)&&(file.exists("stop")))
      return(90L)
   fname <- file.path(output,season)
   if ((!devel)&&(!overwrite)&&(ursa_exists(fname)))
      return(10L)
   if (dir.exists(iceconc)) {
      list1 <- envi_list(path=iceconc,full.names=TRUE)
      iceconc <- list1[grep(season,basename(list1))]
   }
   if (!devel) {
      logname <- paste0(season,".Rlog")
      if ((!overwrite)&&(file.exists(logname)))
         return(20L)
      Fout <- file(logname,"wt")
      writeLines(as.character(Sys.time()),Fout)
      setWindowTitle(title=paste("R  ",season))
   }
   pheno <- gsub("\\d","",season)
   transition <- if (pheno=="sid") c("max","mlt","min") else c("min","frz","max")
   year <- as.integer(gsub("\\D","",season))
   if (F & devel) {
      seed <- sample(100:999,1) ## 625
      cat("seed: ",seed,"\n")
      set.seed(seed)
   }
   period <- seasonality(season)
   transition <- names(period)
   if (F & devel)
      str(period)
   conc <- read_envi(iceconc,cache=TRUE)
   ice <- ursa_value(conc)
   dima <- dim(ice)[1]
   if (!devel)
      res <- array(NA_integer_,dim=c(dima,2),dimnames=list(NULL,c("value","flag")))
   interim <- array(NA_integer_,dim=c(3,length(transition))
                   ,dimnames=list(c("value","flag","qual"),transition))
   cname <- names(conc)
   dname <- as.Date(cname,"%y%m%d")
   if (!devel)
      pb <- ursaProgressBar(dima)
   if (devel)
      ind0 <- ind
   else
      ind0 <- seq_len(dima) |> sample()
   fname <- head(fname,1)
   for (i in ind0) {
      if (!devel)
         setUrsaProgressBar(pb)
      if (devel) {
         cat("------------------------------------------------\n")
         print(c(i=i))
      }
      if (F & devel) {
         loc <- ursa()
         session_grid(loc)
         ursa_value(loc)[i,] <- 1
         loc <- spatialize(loc)
         glance(loc,blank="white",border=0)
      }
      interim[] <- NA_integer_
     # evolution <- switch(pheno,sid="drop",sia="raise")
     # print(transition)
      if (T) {
         if (pheno=="sid")
            evolution <- c("maxdrop","drop","mindrop")
         else
            evolution <- c("minraise","raise","maxraise")
      }
      else {
         if (pheno=="sid")
            evolution <- c("drop","drop","drop")
         else
            evolution <- c("raise","raise","raise")
      }
      names(evolution) <- transition
      for (tr in transition) {
         if (F & tr!="min")
            next
         if (devel)
            print(c(tr=tr))
         ind3 <- which(dname>=period[[tr]]$from & dname<=period[[tr]]$to)
         iname <- cname[ind3]
         y <- ice[i,][ind3]
         x <- dname[ind3]
         if (doy)
            x <- as.integer(x-as.Date(paste0(year,"-01-01"))+1L)
        # print(data.frame(evolution=swi)
         interim[,tr] <- isotonic(x=x,y=y,evolution=evolution[tr]
                                 ,th=th,slen=slen,label=paste(season,i,tr)
                                # ,season=season,transition=tr
                                 ,devel=devel)
      }
      if (devel)
         cat("----- finalizing -----\n")
      if (devel)
         print(interim)
      pre <- c(na.omit(interim["flag",]))
      if (!length(pre))
         next
      if (TRUE) {
         if (all(pre>=931)) {
            if ((pheno=="sid")&&("mlt" %in% names(pre)))
               pre <- pre["mlt"]
            else if ((pheno=="sia")&&("frz" %in% names(pre)))
               pre <- pre["frz"]
            else if (length(pre)>1) {
               print(c(i=i))
               print(pre)
               print("who is decision maker #2 now?")
               if (!devel)
                  writeLines(sprintf("2: i=%d pre[%s]=%d pre[%s]=%d"
                            ,i,names(pre)[1],pre[1],names(pre)[2],pre[2])
                            ,Fout)
               if (F) 
                  pre <- 770L
               else 
                  pre <- switch(pheno,sid=pre["min"],sia=pre["max"],770L)
            }
         }
         else {
            if (T) {
               mqual <- interim["qual",]
               if (devel) {
                  cat("^^^^^^^^^^\n")
                  print(which(mqual==max(mqual)))
                  cat("^^^^^^^^^^\n")
               }
               pre <- interim["value",][which.max(interim["qual",])]
            }
            else {
               if (any(interim["qual",]>75)) {
                  pre <- interim["value",][pre<931 & interim["qual",]>75]
               }
               else {
                  pre <- interim["value",][!is.na(pre) & pre<931]
               }
               pre <- pre[!is.na(pre)]
            }
            if ((pheno=="sid")&&("mlt" %in% names(pre)))
               pre <- pre["mlt"]
            else if ((pheno=="sia")&&("frz" %in% names(pre)))
               pre <- pre["frz"]
            else if (length(pre)>1) {
               print(c(i=i))
               print(pre)
               print("who is decision maker #1 now?")
               if (!devel)
                  writeLines(sprintf("1: i=%d pre[%s]=%d pre[%s]=%d"
                            ,i,names(pre)[1],pre[1],names(pre)[2],pre[2])
                            ,Fout)
               if (F) 
                  pre <- 760L
               else 
                  pre <- switch(pheno,sid=pre["min"],sia=pre["max"],760L)
            }
         }
      }
      if (devel) {
         print(interim[,names(pre)])
      }
      if (devel) {
         if (isTRUE(any(na.omit(interim["value",])<600)))
            break
      }
      if (!devel) {
         qc <- ursa:::.try(res[i,1:2] <- interim[1:2,names(pre)]) ## skip 3 -- qual
         if (!qc) {
            print(c(i=i))
            stop("failed")
         }
      }
   }
   if (!devel)
      close(pb)
   if (!devel) {
      onset <- ursa(res)
      onset <- c(date=doy2date(onset["value"],year=year),onset)
      ursa_write(onset,paste0(fname,".tif"))
      writeLines(as.character(Sys.time()),Fout)
      close(Fout)
     # file.remove(logname)
   }
   0L
}
'isotonic' <- function(x,y,evolution=c("mindrop","drop","maxdrop"
                                      ,"minraise","raise","maxraise")
                      ,th=15,slen=30,label="",devel=TRUE) {
   if (anyNA(y))
      return(NA)
   evolution < match.arg(evolution)
  # transition <- match.arg(transition)
   ret <- c(value=NA_integer_,flag=NA_integer_,qual=-10L)
  # pheno <- gsub("\\d","",season)
  # year <- as.integer(gsub("\\D","",season))
   if (F & devel)
      print(summary(y))
   if (all(y>=th)) {
      ret["flag"] <- 990L
      return(ret)
   }
   if (all(y<th)) {
      ret["flag"] <- 980L
      return(ret)
   }
   if (T & devel)
      print(summary(y))
   names(y) <- as.character(x)
   if (F & devel) {
      cat("names -----------\n")
      print(series(y))
   }
   if (F & devel)
      print(y)
   if (evolution %in% c("mindrop","minraise"))
      L5 <- 100-Iso::ufit(100-y)$y
   else if (evolution %in% c("maxdrop","maxraise"))
      L5 <- Iso::ufit(y)$y
   else {
      L1 <- Iso::ufit(100-y)
      L2 <- Iso::ufit(y)
      if (F & devel)
         print(c(L1mse=L1$mse,L2mse=L2$mse))
      if (L1$mse<L2$mse) {
         L5 <- 100-L1$y
      }
      else {
         L5 <- L2$y
      }
      if (T & devel) {
         str(L1)
         str(L2)
         plot(x,y,type="l",main=sprintf("red(L1)=%.2f blue(L2)=%.2f",L1$mse,L2$mse))
         points(x,100-L1$y,col="red")
         points(x,L2$y,col="blue")
         abline(h=th,lty=2,col="grey")
      }
   }
   names(L5) <- names(y)
   isFound <- FALSE
   L5t <- table(L5>=th)
   if (length(L5t)<2) {
      if ('FALSE' %in% names(L5t))
         ret["flag"] <- 980L
      else if ('TRUE' %in% names(L5t))
         ret["flag"] <- 990L
      else
         stop("never reached")
   }
   else {
      if (L5t[2]<slen) {
         ret["flag"] <- 750L+L5t[2] ## shortterm ice
      }
      if (L5t[1]<slen) {
         ret["flag"] <- 850L+L5t[1] ## shottterm ow
      }
      seg <- rle(L5<th)
      seg$cumlen <- cumsum(seg$lengths)
      if (F & devel)
         str(seg)
      if (devel)
         print(as.data.frame(unclass(seg)))
      nseg <- length(seg[[1]])
      if (nseg==3) {
         lencheck <- min(series(seg$lengths,1))
         if (lencheck<7)
            lseg <- 0L
         else {
            lseg <- seg$lengths[2]*ifelse(!seg$values[2],-1L,1L)
         }
         if (abs(lseg)>slen)
            lseg <- 0
      }
      else
         lseg <- 0L
      if (F & devel)
         print(c(len=series(seg$lengths)))
      ind5 <- c(drop=NA,raise=NA)
      if (length(ind4 <- which(seg$values))) {
         ind4d <- tail(ind4,1)
         if (ind4d>1)
            ind5["drop"] <- seg$cumlen[ind4d-1L]
         ind4r <- head(ind4,1)
         if (ind4r<length(seg$lengths))
            ind5["raise"] <- seg$cumlen[ind4r]
         ind5 <- ind5[!is.na(ind5)]
         if (length(ind5 <- ind5[!is.na(ind5)])==2) {
            if ((evolution %in% c("minraise"))&&(ind5["drop"]<ind5["raise"]))
               ind5 <- ind5["raise"]
            else if ((evolution %in% c("maxdrop"))&&(ind5["drop"]>ind5["raise"]))
               ind5 <- ind5["drop"]
            else if ((evolution %in% c("maxraise"))&&(ind5["raise"]<ind5["drop"]))
               ind5 <- ind5["raise"]
            else if ((evolution %in% c("mindrop"))&&(ind5["raise"]>ind5["drop"]))
               ind5 <- ind5["drop"]
            else if (evolution %in% c("drop"))
               ind5 <- ind5["drop"]
            else if (evolution %in% c("raise"))
               ind5 <- ind5["raise"]
         }
         if (length(ind5)==2) {
            if (devel) {
               cat("**************\n")
               print(cbind(data.frame(evolution),t(ind5)))
               cat("**************\n")
            }
            isFound <- FALSE
            ret["flag"] <- 970L
         }
         else if (grepl("drop",evolution)) {
            if (ind4d>1) {
               lens <- seg$lengths[c(ind4d-1L,ind4d)]
              # print(c(lens=lens,lseg=lseg,ind5=ind5))
               if ((lseg!=0)&&(any(lens>=7))) {
                  isFound <- TRUE
                  ret["flag"] <- if (lseg<0) -lseg+700L else lseg+800L
               }
               else if (all(lens>=7)) {
                  isFound <- TRUE
                  ret["flag"] <- 930L
               }
            }
         }
         else if (grepl("raise",evolution)) {
            if (ind4r<length(seg$lengths)) {
               lens <- seg$lengths[c(ind4r,ind4r+1L)]
               if ((lseg!=0)&&(any(lens>=7))) {
                  isFound <- TRUE
                  ret["flag"] <- if (lseg<0) -lseg+700L else lseg+800L
               }
               else if (all(lens>=7)) {
                  isFound <- TRUE
                  ret["flag"] <- 930L
               }
            }
         }
         if (devel & isFound)
            print(c(ind5=ind5))
      }
   }
   if (!isFound) {
      if (is.na(ret["value"])) {
         ymean <- mean(y)
         if (ymean>=th)
            ret["flag"] <- 990L
         else
            ret["flag"] <- 980L
      }
   }
   if (isFound) {
      if (grepl("drop",evolution))
         ind6 <- ind5-0L ## < th%, last ice-day
      else
         ind6 <- ind5+1L ## <th% ## `+0L` last ow-day, `+1L` first ice-day
      if (F & devel)
         print(L5)
      if (devel)
         print(ind6)
      ret["qual"] <- as.integer(round((0.5-abs(0.5-ind6/length(y)))*1000))
     # ret["qual"] <- min(ind6-1L,length(y)-ind6)
   }
   if (isFound) {
      ret["value"] <- x[ind6]
      if (is.na(ret["flag"]))
         ret["flag"] <- 930L
   }
   if (T & devel) {
      xp <- x
     # print(x[ind6])
      plot(xp,y,type="l",xlab="Day of Year",ylab="Ice Concentration",main=label)
     # lines(xp,predict(loess(y~xp,span=2/12)),col="orange")
      points(xp,L5)
     # if (isFound)
     #    points(xp[ind4s],L5[ind4s],color="brown")
      abline(h=th,lty=2,col="#000000A0")
      if (isFound)
         abline(v=ret["value"],lty=2,col="#000000A0")
   }
   ret
}
'extent' <- function(season,date) {
   year <- as.integer(gsub("\\D","",basename(season)))
   pheno <- gsub("\\d","",basename(season))
   Dec31 <- as.integer(as.Date(paste0(year,"-01-01"))-1L)
   a <- ursa_read(season)
   value <- a[a<800]
   flag <- a[a>800]
   flag[flag==820] <- 880
   flag[flag==830] <- 880
   if (missing(date))
      doy <- quantile(na.omit(c(ursa_value(value))),runif(1,min=0,max=1))
   else if (is.character(date)) {
      if (is.na(d3 <- try(as.Date(date,"%Y-%m-%d")))) {
         if (is.na(d3 <- try(as.Date(date,"%Y%m%d")))) {
            if (is.na(d3 <- try(as.Date(date,"%y%m%d")))) {
               if (is.na(d3 <- try(as.Date(date,"%Y%j")))) {
                  doy <- quantile(na.omit(c(ursa_value(value))),runif(1,min=0,max=1))
               }
            }
         }
      }
      if ((!is.na(d3))&&(is.character(date))) {
         doy <- as.integer(d3)-Dec31
      }
   }
   print(value)
   print(date)
   print(doy)
   ret <- if (pheno=="sid") value>=doy else value<=doy
   ret[is.na(ret)] <- flag==880
   ursa_write(ret,"inconsistent.tif")
   display(list(ret,flag==880))
   0L
}
'doy2date' <- function(season,year) {
   if (is.ursa(season)) {
      a <- season
      if (missing(year)) {
         year <- as.integer(gsub("\\D*(\\d{4})\\D*","\\1",names(a)))
         if (!(year %in% seq(1970,2070)))
            stop("need to specify 'year=' argument")
      }
      year <- as.integer(year)
   }
   else
      a <- ursa_read(season)
   if (missing(year))
      year <- as.integer(gsub("\\D","",basename(season)))
   Dec31 <- as.integer(as.Date(paste0(year,"-01-01"))-1L)
   value <- a[a<800]
   flag <- a[a>800]
   value <- value+Dec31
   value[is.na(value)] <- flag
   value
}
'date2doy' <- function(season) {
   if (is.ursa(season))
      a <- season
   else
      a <- ursa_read(season)
   value <- a[a>900]
   flag <- a[a<900]
   mval <- as.Date(band_min(value),origin=ursa:::.origin())
   year <- min(as.integer(format(mval,"%Y")),na.rm=TRUE)
   Dec31 <- as.integer(as.Date(paste0(year,"-01-01"))-1L)
   value <- value-Dec31
   value[is.na(value)] <- flag
   value
}
'cleanse' <- function(season,summer=7,winter=21,unit=c("value","date"),flag=T) {
   unit <- match.arg(unit)
   if (is.character(season)) {
      season <- ursa(season)
      if (all(c("v","d","kind") %in% names(season))) {
         season <- season[c("v","d","kind")]
         names(season) <- c("value","date","flag")
      }
   }
   if (length(season)==1)
      return(season)
   d <- as.Date(band_min(season["date"]),origin=ursa:::.origin())
   year <- as.integer(format(d,"%Y"))
   mind <- as.integer(format(d,"%j"))
   pheno <- ifelse(mind<150,"sid","sia")
   a <- season[unit]
   f <- season["flag"]
   indOW <- f>700 & f<700+winter
   indIce <- f>800 & f<800+summer
   f[indOW] <- 980
   f[indIce] <- 990
   f[f>600 & f<900] <- 930
   a[indOW] <- NA
   a[indIce] <- NA
   if (flag)
      a[is.na(a)] <- f
   names(a) <- paste0(pheno,year)
   a
}
'season_length' <- function(phenology=c("./output",".")[1],season=c("ow","ice")) {
   season <- match.arg(season)
   nature <- interleaving()
   ind1 <- grep("sid",nature)
   ind2 <- ind1+1
   ind2 <- ind2[ind2<=length(nature)]
   ind1 <- head(ind1,length(ind2))
   ow <- data.frame(sid=nature[ind1],sia=nature[ind2]
                   ,ow=gsub("\\D+","ow",nature[ind2]))
  # print(series(ow))
   ind1 <- grep("sia",nature)
   ind2 <- ind1+1
   ind2 <- ind2[ind2<=length(nature)]
   ind1 <- head(ind1,length(ind2))
   ice <- data.frame(sia=nature[ind1],sid=nature[ind2]
                    ,ice=gsub("\\D+","ice",nature[ind2]))
  # print(series(ice))
   session_grid("ref/dist2land.tif")
   if (F) {
      res <- ursa(bandname=ice$ice)
      for (i in seq(res) |> sample()) {
         cond1 <- ursa_exists(p1 <- file.path(phenology,ice$sia[i]))
         cond2 <- ursa_exists(p2 <- file.path(phenology,ice$sid[i]))
         if (!cond1 | !cond2)
            next
         v1 <- cleanse(p1,unit="date",flag=F)
         v2 <- cleanse(p2,unit="date",flag=F)
         print(ice[i,])
        # print(c(sia=v1,sid=v2))
         res[i] <- v2-v1
      }
      res <- res[!band_blank(res)]
      print(res)
      if (F)
         ursa_write(res,file.path(phenology,"season_ice.tif"))
   }
   if (T) {
      res <- ursa(bandname=ow$ow)
      for (i in seq(res) |> sample()) {
         cond1 <- ursa_exists(p1 <- file.path(phenology,ow$sid[i]))
         cond2 <- ursa_exists(p2 <- file.path(phenology,ow$sia[i]))
         if (!cond1 | !cond2)
            next
         v1 <- cleanse(p1,unit="date",flag=F)
         v2 <- cleanse(p2,unit="date",flag=F)
         print(ow[i,])
        # print(c(sid=v1,sia=v2))
         res[i] <- v2-v1
      }
      res <- res[!band_blank(res)]
      bres <- band_stat(res)
      print(bres)
      print(sum(bres$n))
      if (F)
         ursa_write(res,file.path(phenology,"season_ow.tif"))
      if (T) {
         ow <- res[res<=0]
         ow <- as.data.frame(ow)
         if (nrow(ow)>0) {
            owd <- ow[,grep("ow",colnames(ow))]
            sn <- apply(owd,1,function(x) names(x[!is.na(x)]))
            sd <- apply(owd,1,function(x) x[!is.na(x)])
           # da <- rowSums(owd,na.rm=TRUE)
           # print(da)
            sow <- cbind(ow[,grep("ow",colnames(ow),invert=TRUE)],season=sn,value=sd)
            sow <- sow[order(sow$season,decreasing=TRUE),]
            rownames(sow) <- NULL
            print(sow)
            sink("diff.Rout")
            print(sow)
            sink()
         }
      }
   }
   0L
}
'manage_converters' <- function() {
   ##~ a <- extent(season="output/sia2012",date="2012-11-15")
   o <- cleanse("./output/sia2020")
   s1 <- ursa_read("output/sia2012")
   print(tail(as.table(s1),12))
   s2 <- doy2date(s1)
   print(head(as.table(s2),12))
   s3 <- date2doy(s2)
   print(tail(as.table(s3),12))
   print(identical(s1,s3))
   0L
}
'plot_onsets' <- function(path=".",objs) {
   if (F) {
      v1 <- round(runif(120,min=60,max=240))
      v1 <- ursa_dummy(1,min=60,max=240)
      p1 <- colorize(v1,stretch="julian",ramp=F,interval=T)
      print(p1)
      str(p1)
     # d1 <- .deintervale(p1)
     # print(p1)
     # display(p1)
      q()
   }
   c1 <- sample(seq(0,359),1)
   print(c(c1=c1))
   c2 <- c1+180
   c3 <- c2+180
   b1 <- 127
   b2 <- 223
   h <- 1.5
   palIce <- cubehelix(1,weak=c2,bright=b1,hue=h)
   palOw <- cubehelix(1,weak=c1,bright=b2,hue=h)
   palSID <- cubehelix(11,weak=c1,rich=c2,light=b1,dark=b2,hue=h) |> head(-1) |> tail(-1)
   palSIA <- cubehelix(11,weak=c2,rich=c3,light=b2,dark=b1,hue=h) |> head(-1) |> tail(-1)
   print(palIce)
   print(palOw)
   print(palSIA)
   print(palSID)
   if (F) {
      ps <- pixelsize()
      display(list(colorize(ps,pal=palSID,ramp=FALSE)
                  ,colorize(ps,pal=palSIA,ramp=FALSE)))
      q()
   }
  # display(pixelsize(),ramp=FALSE,pal=palSID)
  # q()
  # pal <- cubehelix(5,weak=120,rotate=0,hue=1.5)
  # ps <- colorize(pixelsize(),pal=c(palOw,palIce),ramp=FALSE)
  # print(ps)
  # display(ps)
   src <- strsplit(objs,split="\\s+")[[1]]
  # res <- lapply(src,function(fname) ursa_read(file.path(path,fname))) |> as.ursa()
   res <- lapply(src,function(fname) cleanse(file.path(path,fname))) |> as.ursa()
   pheno <- gsub("\\d","",names(res))
   b1 <- split(res,pheno)
   b2 <- lapply(b1,function(x) {
      value <- x[x<800]
      print(value)
      isSID <- all(band_max(value)<350)
      isSIA <- all(band_min(value)>150)
      print(c(SID=isSID,SIA=isSIA))
      flag <- x[x>800]
      if (F) {
         flag[flag==820 | flag==830] <- 880
         flag[flag==810 | flag==840] <- 870
         flag[flag==850 | flag==860] <- NA
      }
      flag <- colorize(flag,value=c(980,990),name=c("OW","Ice"),pal=c(palOw,palIce))
      ctF <- ursa_colortable(flag)
      value <- colorize(value,ramp=FALSE
                       ,pal=if (isSID) palSID else palSIA
                       ,interval=!TRUE,stretch="julian",inv=F)
      ctV <- ursa_colortable(value)
      value <- value+1L
      if (isSIA)
         flag <- 1-flag
      flag <- flag*(band_max(value)+1)
      value[is.na(value)] <- flag
      if (isSID)
         ct <- c(head(ctF,1),ctV,tail(ctF,1))
      else
         ct <- c(tail(ctF,1),ctV,head(ctF,1))
      value <- colorize(value,pal=ct,name=names(ct))
     # write_envi(value[1],"C:/tmp/sia");q()
      value
   })
  # print(b2)
   ct <- lapply(b2,ursa_colortable)
   b3 <- lapply(unname(b2),function(x) as.list(x)) |> do.call(c,args=_)
   nature <- apply(expand.grid(pheno=c("sid","sia"),year=seq(1988,2024)),1
                  ,function(x) paste(x,collapse=""))
   b3 <- b3[na.omit(match(interleaving(),names(b3)))]
   cat("-----------------\n")
   print(b3)
   compose_open(b3,layout=c(1,length(b3))
              # ,legend=list("left","right")
               ,legend=list("top","bottom")
               ,fileout="./article/assets/onsets.png"
               )
   for (i in seq_along(b3)) {
      compose_panel(b3[[i]],blank="white")
      panel_annotation(names(b3[[i]]))
   }
   compose_legend(ct)
   compose_close(bpp=8)
   0L
}
'interleaving' <- function(ref="") {
   nature <- apply(expand.grid(pheno=c("sid","sia")
                              ,year=seq(1988,2024)
                             # ,year=seq(2019,2021)
                              ),1,function(x) paste(x,collapse=""))
   if ((is.null(ref))||(!dir.exists(ref)))
      return(nature)
   if (!length(list1 <- envi_list(ref))) {
      list1 <- dir(path=ref)#,pattern="si[ad]\\d{4}")
   }
   list1 <- list1[grep("si[ad]\\d{4}",basename(list1))]
   list2 <- gsub(".*(si[ad]\\d{4}).*","\\1",basename(list1))
   list2 <- list2[na.omit(match(nature,list2))]
   list2
}
'check_types' <- function(phenology=c("output","../v43nt/phenology")[1]) {
  # list1 <- dir(path=phenology,pattern="\\.tif",full.names=TRUE)
  # list1 <- list1[match(gsub("(\\..+)$","",basename(list1)),interleaving())]
   list1 <- interleaving(ref=phenology) |> tail(8) # |> head(-1)
   list2 <- lapply(file.path(phenology,list1),function(fname) {
      a <- cleanse(fname)
      a[a<900] <- 930
      as.table(a)
   }) |> do.call(rbind,args=_)
   rownames(list2) <- list1
   print(list2)
  # print(series(list1))
   0L
}
'manage' <- function(devel=FALSE) {
   ref <- ursa_read("ref/dist2land.tif")
   iceconc <- c("../v43nt/conc","./iceconc")[1]
   if (!devel) {
      list2 <- interleaving(iceconc) |> rev() |> head(1) #|> head(4)
      sapply(list2,phenology,iceconc=iceconc,overwrite=length(list2)==1)
   }
   else {
      fname <- c("sid2022","diff","diff.Rout")[3]
      ind0 <- if (F) NULL
              else if (F) ursa:::.getIndex(ref,x=1212500,y=-12500)
              else if (F) 99911
              else if (F & grepl("si[ad]\\d{4}",fname) & ursa_exists(fname)) {
                 a <- ursa_read(fname)["flag"] #,engine="rgdal")
                 a <- a[a>600]
                 print(as.table(a))
                 a <- a[a==960] ## 950 960
                 xy <- as.data.frame(a)
                 xy <- xy[sample(seq_len(nrow(xy)),1),]
                 ursa:::.getIndex(ref,x=xy$x,y=xy$y)
              }
              else if (F & grepl("diff",fname) & ursa_exists(fname)) {
                 a <- ursa_read(fname)
                 a <- a[a<=0]
                 print(as.table(a))
                 xy <- as.data.frame(a)
                 xy <- xy[sample(seq_len(nrow(xy)),1),]
                 ursa:::.getIndex(ref,x=xy$x,y=xy$y)
              }
              else if (T & file.exists(fname)) {
                 da <- read.table(fname)
                 da <- da[sample(seq_len(nrow(da)),1),]
                 print(da)
                 ind <- ursa:::.getIndex(ref,x=da$x,y=da$y)
                 attr(ind,"fname") <- paste0(c("sid","sia"),gsub("\\D","",da$season))
                 ind
                # print(ind0)
                # fname <- paste0(c("sid","sia"),gsub("\\D","",da$season))
                # print(fname)
              }
              else sample(seq_len(prod(dim(ref)[1:2])),100) 
      if (is.null(attr(ind0,"fname"))) {
         attr(ind0,"fname") <- interleaving(iceconc) |> sample(1)
      }
      if (T & length(ind0)==1) {
         xy <- spatialize(data.frame(t(coord_cr(ref,ind0))),crs=ursa_crs(ref))
         session_grid(ref)
         glance(xy)
         session_grid(ref)
      }
     # phenology(interleaving(iceconc) |> sample(1),ind=ind0)
     # phenology("sid2014",ind=ind0)
     # phenology("sia2014",ind=ind0)
      sapply(attr(ind0,"fname"),phenology,ind=ind0)
   }
   0L
}
invisible(if (ursa:::.argv0()=="phenology.R") {
  # prepare("./iceconc","sia2012")
   manage(devel=F)
  # plot_onsets(path="output","sia2019 sid2020 sia2020 sid2021")
  # season_length()
  # check_types()
})
if (requireNamespace("plutil")) {
   plutil::ursula(3)
}
