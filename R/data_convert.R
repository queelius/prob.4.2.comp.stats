#' Convert counts to individual responses
#'
#' There are many ways to represent the same data. In problem 4.2, the data
#' is given in the form of counts, e.g., counts[j] represents the number
#' of people who responded with j risky sexual encounters.
#' The more natural form of this data is as a vector of indivual responses,
#' e.g., x[k] is the k-th persons response to the number of risky sexual
#' encounters.
#' @param counts Count data to convert
#' @keywords counts multinomial response
#' @export
em.counts_to_responses <- function(counts)
{
  # say counts = (379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1)
  # this means:
  #     379 responded 0 encounters
  #     299 responded 1 encounters
  #     222 responded 2 encounters
  #     .
  #     .
  #     .
  #     1 responded 16 encounters
  data <- NULL
  for (i in 1:length(ns))
  {
    data <- append(data,rep((i-1),ns[i]))
  }
  data
}

#' Convert individual responses to count data
#' @param data Vector of responses
#' export
em.responses_to_counts <- function(data)
{
  counts <- NULL
  for (i in 0:16)
  {
    ni <- data[data == i]
    l <-length(ni)
    counts <- append(counts,l)
  }
  counts
}
