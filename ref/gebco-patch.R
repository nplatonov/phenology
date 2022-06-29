plutil::ursula()
a <- ursa_read("gebco")
a
a <- focal_special(a,"gaus",size=3,verbose=TRUE,fillNA=TRUE,cover=0)
a
ursa_write(a,"gebco.tif")
#a <- spatialize(is.na(a))
#glance(a)
# ursa:::widgetize(mapview::mapview(a))


