'level1' <- function() {
  # e1 <- new.env(parent = baseenv())  # this one has enclosure package:base.
  # e2 <- new.env(parent = e1)
  # assign("a", 3, envir = e1)
   for (ext in c("min","max")) {
      label <- switch(ext,min="minimum",max="maximum")
     # assign("desc",label)
     # desc <- label
      ret <- level2(ext)
      print(ret)
   }
   0L
}
'level2' <- function(x) {
   a <- get("desc",envir=parent.frame(1))
   print(a)
   paste("argument is:",x)
}
invisible({
  level1()
})