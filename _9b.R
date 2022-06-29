plutil::ursula()
(a1 <- ursa("sid2014"))
as.table(a1["flag"])
q()
v1 <- a1["value"]
f1 <- a1["flag"]
v1[f1>700 & f1<=702] <- NA
v1[f1>800 & f1<=802] <- NA
a1["value"] <- v1
a1["value"]
a2 <- ursa("output/sid2014")
a2 <- c(value=a2[a2<600],flag=a2[a2>600])
as.table(a2["flag"])
a2["value"]
q()

f2[f2<600] <- NA
as.table(f2)
a2[a2>600] <- NA
a2
m <- !is.na(a1["value"])+!is.na(a2)
m
