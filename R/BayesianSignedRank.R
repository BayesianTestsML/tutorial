BayesianSignedRank <- function(diffVector, rope_min, rope_max, samples=3e4) {
  
  # Give half weight to pseudo observation to reflect prior
  weights <- c(0.5, rep(1,length(diffVector)))

  # Add difference of 0 to favour rope hypothesis via prior
  diff_vec <- c (0, diffVector)

  # Pre-compute matrix of summed difference
  sum_mat <- outer(diffVector, diffVector, "+")

  # Identify indices where sum of differences support region
  max_rope <- 2 * rope_max
  min_rope <- 2 * rope_min
  right_idx <- which(sum_mat > max_rope, arr.ind = TRUE)
  rope_idx <- which((sum_mat >= min_rope) & (sum_mat <= max_rope), arr.ind = TRUE)
  left_idx <- which(sum_mat < min_rope, arr.ind = TRUE)

  # Draw samples
  sampled_weights <- MCMCpack::rdirichlet(samples, weights)

  # Preallocate posterior probability
  win_left <- vector(length = samples, mode = "double")
  win_rope <- vector(length = samples, mode = "double")
  win_right <- vector(length = samples, mode = "double")

  for (rep in 1:samples){
    weights_rep <- sampled_weights[rep,]
    prod_mat <- outer(weights_rep, weights_rep, "*")
    win_right[rep] <- win_right[rep] + sum(prod_mat[right_idx])
    win_rope[rep] <- win_rope[rep] + sum(prod_mat[rope_idx])
    win_left[rep] <- win_left[rep] + sum(prod_mat[left_idx])
    max_wins <- max(win_right[rep], win_rope[rep], win_left[rep])
    is_right <- (win_right[rep] == max_wins) * 1
    is_rope <- (win_rope[rep] == max_wins) * 1
    is_left <- (win_left[rep] == max_wins) * 1
    n_winners <- is_right + is_rope + is_left
    win_right[rep] <- is_right / n_winners
    win_rope[rep] <- is_rope / n_winners
    win_left[rep] <- is_left / n_winners
  }
  
  results <- list(
    "winLeft" = mean(win_left),
    "winRope" = mean(win_rope),
    "winRight" = mean(win_right)
  )
  return (results)
}
