library(devtools)
library(roxygen2)

# setwd("C:/MBender/Arbeit/Github")
setwd("Z:/Aktuell/Eigene Dateien/Eigene Dateien_Marc/R/Github")

# package name
package <- "mlekb"

# create package
# create(package)

# update package
setwd(paste("./", package, sep = ""))
document()

# check package
setwd("..")
check(package)

# install package
install(package)






