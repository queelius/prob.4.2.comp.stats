#' Convert counts to individual responses
#'
#' There are many ways to represent the same data. In problem 4.2, the data
#' is given in the form of counts, e.g., counts[j] represents the number
#' of people who responded with j risky sexual encounters.
#' The more natural form of this data is as a vector of indivual responses,
#' e.g., x[k] is the k-th persons response to the number of risky sexual
#' encounters.
#' @param counts Count data to convert
#' @keywords counts multinomial responses
#' @export
em.counts_to_responses <- function(counts)
{
  # say counts = (379,299,...), then 379 responded 0 encounters, 299 responded
  # 1 encounter, and so on.
  data <- NULL
  for (i in 1:length(counts))
  {
    data <- append(data,rep((i-1),counts[i]))
  }
  data
}

#' Convert individual responses to count data
#'
#' This is the inverse of em.counts_to_responses, e.g.,
#'   em.responses_to_counts(em.counts_to_responses(counts)) == counts
#' and
#'   em.counts_to_responses(em.responses_to_counts(data)) == data.
#' @param data Vector of responses
#' keywords counts multinomial responses
#' export
em.responses_to_counts <- function(data)
{
  counts <- NULL
  for (i in 0:16)
  {
    counts <- append(counts,length(data[data == i]))
  }
  counts
}