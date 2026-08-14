setwd("C:/Users/01459189/OneDrive/phd/mixed-stratification")
shinylive::export(appdir = "shiny", destdir = "docs")
httpuv::runStaticServer("docs/")