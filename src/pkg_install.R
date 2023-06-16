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
   'GSVA'

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


# stable version on CRAN
install.packages("rmarkdown",
  dependencies = TRUE,
  INSTALL_opts = c("--no-lock")
)
# or development version on GitHub
# remotes::install_github('rstudio/bookdown')
