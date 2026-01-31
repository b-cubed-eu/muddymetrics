#' Save Download Key to Disk
#' @param key The GBIF download key.
#' @param path File path to save the key.
save_download_key <- function(key, path = ".gbif_download_key") {
  writeLines(as.character(key), path)
  message("Download key saved to: ", path)
}

#' Load Download Key from Disk
#' @param path File path to the saved key.
#' @return The key as character, or NULL if not found.
load_download_key <- function(path = ".gbif_download_key") {
  if (file.exists(path)) {
    return(readLines(path, warn = FALSE)[1])
  }
  return(NULL)
}

#' Robust File Download with Resume
#' @description Uses system curl to download a file with resume support.
#' @param url The download URL.
#' @param dest_path The local destination path.
download_robust <- function(url, dest_path) {
  message("Starting robust download via system curl...")
  message("Destination: ", dest_path)
  
  # -L: Follow redirects
  # -C -: Resume transfer from where it left off
  # -o: Output file
  cmd <- sprintf("curl -L -C - \"%s\" -o \"%s\"", url, dest_path)
  
  system(cmd)
}

