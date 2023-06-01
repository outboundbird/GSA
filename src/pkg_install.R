 env_pkgs <- .packages(all.available = T)
 attached <- (.packages())
 req_libs <- c(
   "GEOquery",
   "logger"
 )
 to_install <- req_libs[!req_libs %in% env_pkgs]


if (length(to_install)) {
  lapply(to_install, function(pkg) {
    tryCatch(
      {
        message(sprintf("Installing %s", pkg))
        install.packages(to_install,
          dependencies = TRUE,
          INSTALL_opts = c("--no-lock")
        )
      },
      error = function(e) {
        message(e)
        BiocManager::install(pkg)
      }
    )
  })
}
