#' Check for GBIF Credentials
#'
#' @description
#' Verifies that the required GBIF credentials (user, pwd, email) are set
#' in the R options.
#'
#' @return Logical. TRUE if all credentials are set, stops otherwise.
#' @export
check_gbif_credentials <- function() {
  user <- getOption("gbif_user")
  pwd <- getOption("gbif_pwd")
  email <- getOption("gbif_email")
  
  if (is.null(user) || is.null(pwd) || is.null(email)) {
    stop(
      "GBIF credentials are not set. Please set them using:\n",
      "options(gbif_user = '...', gbif_pwd = '...', gbif_email = '...')\n",
      "or add them to your .Renviron file."
    )
  }
  return(TRUE)
}

