#' @export
calc_odds <- function(final_acc,prac_score) {
  # if we take the log of (0/1), we get -Inf
  # if we take the log of (1/0), we get Inf
  # if we take the log of (0/0), we get NaN
#   if ( length(prac_score[prac_score == 1]) != length(prac_score[prac_score == 0]) ) {
#     stop("invalid input vector, lengths of groups are not equal)

#   }
  if (!all(prac_score %in% 0:1)) {
    stop("Not all values in practice score vector are zero or one")
  }
  odds <- log(final_acc[prac_score == 1]/final_acc[prac_score == 0])
  return(odds )

}
