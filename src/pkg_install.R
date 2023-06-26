 env_pkgs <- .packages(all.available = T)
 attached <- (.packages())
 req_libs <- c(
   "GEOquery",
   "logger",
   "meshes",
   "RDAVIDWebService",
   "MeSH.Hsa.eg.db",
   "bookdown",
   'rmarkdown',
   'GSVA',
   'GSAR'

 )
to_install <- req_libs[!req_libs %in% env_pkgs]
to_install <- 'GSVA'

if (length(to_install)) {
  lapply(to_install, function(pkg) {
    tryCatch(
      {
        message(sprintf("Installing %s", pkg))
        BiocManager::install(pkg)

      },
      error = function(e) {
        message(e)
                install.packages(to_install,
          dependencies = TRUE,
          INSTALL_opts = c("--no-lock")
        )

      }
    )
  })
}


req_libs <- c("shinybusy", "rmarkdown")
# stable version on CRAN
lapply(req_libs, function(pkg){
  tryCatch({
      install.packages(pkg,
    dependencies = TRUE,
    INSTALL_opts = c("--no-lock")
  )
  }, error = function(e){
    warning(paste('Error installing package', pkg))
  })
})

# or development version on GitHub
# remotes::install_github('rstudio/bookdown')
