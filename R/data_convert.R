#' Convert counts to individual responses
#'
#' There are many ways to represent the same data. In problem 4.2, the data
#' is given in the form of counts. This is the inverse of
#' em.counts_to_responses, e.g.,
#' ```
#'  xs <- em.counts_to_responses(counts)
#'  em.responses_to_counts(xs) == counts
#'  em.counts_to_responses(em.responses_to_counts(xs)) == xs.
#'```
#' @param counts count data
#' @return the count data converted to response data, say counts = (379,299,...),
#'         then 379 responded 0 encounters, 299 responded 1 encounter, ...
#' @keywords counts multinomial responses
#' @examples
#'   # let counts be the count data
#'   counts[j] # denotes number of respondents with j risky sexual encounters.
#'   xs <- em.counts_to_responses(counts)
#'   xs[k] # denotes the response of the i-th person
#' @export
em.counts_to_responses <- function(counts)
{
  data <- NULL
  for (i in 1:length(counts))
  {
    data <- append(data,rep((i-1),counts[i]))
  }
  data
}

#' Convert individual response data to count data
#'
#' This is the inverse of em.counts_to_responses, e.g.,
#'   em.responses_to_counts(em.counts_to_responses(counts)) == counts
#' and
#'   em.counts_to_responses(em.responses_to_counts(data)) == data.
#' @param data response data
#' @returns response data converted to count data
#' @keywords counts multinomial responses
#' @export
em.responses_to_counts <- function(data)
{
  counts <- NULL
  for (i in 0:16)
  {
    counts <- append(counts,length(data[data == i]))
  }
  counts
}
