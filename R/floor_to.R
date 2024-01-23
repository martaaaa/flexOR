#' @title floor_to: Rounds a number down
#' @description
#' This function rounds a number down to the specified multiple.
#' 
#' @param x The number to be rounded.
#' @param to The multiple to which the number should be rounded down. Default is 1.
#' 
#' @return The number rounded down to the specified multiple.
#' 
#' @examples
#' floor_to(7, 3) # Returns 6, as 6 is the largest multiple of 3 less than or equal to 7.
#' floor_to(5, 2) # Returns 4, as 4 is the largest multiple of 2 less than or equal to 5.
#' 
#' @importFrom base floor
#' @export
#'
floor_to <- function(x, to = 1) {
  floored <- floor(x / to) * to
  return(floored)
} #floor_to