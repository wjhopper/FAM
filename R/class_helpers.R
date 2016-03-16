#' FAM classes
#'
#' FAM provides these helper functions to construct objects belonging to S3 classes
#' from the \code{\link{FAM}} package, or to quickly check if an R object
#' inhereits from one of the FAM package's classes.
#'
#' It is unlikely the end-user should need to call these functions, but this page
#' serves to document their behavior should they need to be used externally from
#' the \code{\link{FAM}} package, and to give more information on the structure
#' of the S3 classes implemented by the \code{\link{FAM}} package.
#' @name FAM_classes
#'
#' @param data An R object that is, or can be coerced to, a data frame.
#' @param x an R object
#'
#' @return The \code{as*} functions (e.g.,\code{as.LB4L_IV}, \code{as.LB4L_IV_summary})
#' will construct data frames (by coercing their inputs with \code{\link{as.data.frame}})
#' with the primary class specificed by the function (e.g., \code{as.LB4L_IV} gives the
#' data frame the primary class \code{"LB4L)_IV"}).
#'
#' The \code{is*} functions (e.g.,\code{is.LB4L_IV}, \code{is.LB4L_IV_summary})
#' return logical scalars. They will return \code{TRUE} if the object has that
#' specific class (i.e.,\code{LB4L_IV}, \code{LB4L_IV_summary}, etc.) as it's
#' first value of the \code{class} attribute.

NULL
