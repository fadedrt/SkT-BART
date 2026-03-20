#' @title BART Tree Prediction and MCMC Utilities
#' @description Core functions for tree traversal, prediction, and MCMC sampling.
#' @author Your Name
#' @license MIT

# --- Prediction Logic ---------------------------------------------------------

#' Fill tree details using C++ backend
#' @param curr_tree A list containing the tree_matrix.
#' @param X Feature matrix.
fill_tree_details <- function(curr_tree, X) {
  # Direct call to the compiled C++ traversal engine
  res <- fill_tree_details_cpp(curr_tree$tree_matrix, X)
  
  return(list(
    tree_matrix = res$tree_matrix,
    node_indices = res$node_indices
  ))
}

#' Get predictions from a single tree or an ensemble
#' @param trees A list of tree objects or a single tree.
#' @param X Input feature matrix.
#' @param single_tree Logical; whether to process as a single tree.
#' @export
get_predictions <- function(trees, X, single_tree = FALSE) {
  
  # Standardize input: prevent nesting issues in list of lists
  if (is.null(names(trees)) && (length(trees) == 1)) {
    trees <- trees[[1]]
  }
  
  if (single_tree) {
    # Case 1: Tree is just a root node (no splits)
    if (nrow(trees$tree_matrix) == 1) {
      predictions <- rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Case 2: Multi-node tree - Loop through node indices to get predictions
      predictions <- rep(NA, nrow(X))
      
      # Determine which node each observation in X falls into
      curr_X_node_indices <- fill_tree_details(trees, X)$node_indices
      unique_node_indices <- unique(trees$node_indices)
      
      # Explicit loop to fill in mu values per leaf node
      for (i in seq_along(unique_node_indices)) {
        target_node <- unique_node_indices[i]
        predictions[curr_X_node_indices == target_node] <- trees$tree_matrix[target_node, 'mu']
      }
    }
  } else {
    # Case 3: Ensemble of trees (Additive recursion)
    partial_trees <- trees
    partial_trees[[1]] <- NULL # Remove the first tree to recurse
    
    predictions <- get_predictions(trees[[1]], X, single_tree = TRUE) +
      get_predictions(partial_trees, X, single_tree = (length(partial_trees) == 1))
  }
  
  return(predictions)
}

# --- Sampling & Ancestry ------------------------------------------------------

#' Decide which move to make in the MCMC step
#' @param curr_tree The current tree structure.
#' @param iter Current MCMC iteration.
#' @param nburn Total burn-in iterations.
sample_move <- function(curr_tree, iter, nburn) {
  # Force 'grow' if tree is a root or during the early 10% of burn-in
  if (nrow(curr_tree$tree_matrix) == 1 || iter < max(floor(0.1 * nburn), 10)) {
    return('grow')
  } else {
    return(sample(c('grow', 'prune', 'change'), 1))
  }
}

#' Get all descendants of a specific node
get_children <- function(tree_mat, parent) {
  if (as.numeric(tree_mat[parent, 'terminal']) == 1) {
    return(parent)
  } else {
    child_l <- as.numeric(tree_mat[parent, 'child_left'])
    child_r <- as.numeric(tree_mat[parent, 'child_right'])
    
    return(c(parent, 
             get_children(tree_mat, child_l), 
             get_children(tree_mat, child_r)))
  }
}
