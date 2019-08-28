.onLoad <- function(libname, pkgname) {
  jdk.location <- list.files(path = "C:/Program Files/Java/", pattern = "^jdk.*", full.names = TRUE)[1]
  if (is.na(jdk.location)) {
    stop("You must install a Java JDK.")
  }

  Sys.setenv(JAVA_HOME = jdk.location)

  if (!file.exists(paste0(system.file("java", package = "dismo"), "/maxent.jar"))) {
    file.copy(
      from = "maxent.jar",
      to = paste0(system.file("java", package = "dismo"), "/maxent.jar")
    )
  }

  rJava::.jinit(parameters = "-Xmx2g")
}
