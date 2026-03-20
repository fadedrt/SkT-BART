fill_tree_details = function(curr_tree, X) {
  # 调用刚才编译的 C++ 后端
  res <- fill_tree_details_cpp(curr_tree$tree_matrix, X)
  
  return(list(tree_matrix = res$tree_matrix,
              node_indices = res$node_indices))
}
# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE) {
  
  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]
  
  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] =
          trees$tree_matrix[unique_node_indices[i], 'mu']
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }
  
  return(predictions)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

# Resample a single element from a vector
resample <- function(x, ...) x[sample.int(length(x), size = 1), ...]

# Update Dirichlet weights for splitting probabilities
update_s_  <- function(var_count, p, alpha_s) {
  shape <- alpha_s / p + var_count
  temp  <- rgamma(length(shape), shape, rate = 1)
  temp / sum(temp)
}

# Count the number of distinct covariates used in internal nodes
get_number_distinct_cov <- function(tree) {
  # Select the rows corresponding to internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 0)
  # Get unique covariates used in splitting
  num_distinct_cov = length(unique(tree$tree_matrix[which_terminal,'split_variable']))
  return(num_distinct_cov)
}

# Decide which type of move to make in tree MCMC
sample_move = function(curr_tree, i, nburn) {
  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1 * nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'), 1)
  }
  return(type)
}

# Get ancestors for terminal nodes
get_ancestors = function(tree) {
  save_ancestor = NULL
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  
  if(nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL, ancestor = NULL)
  } else {
    for (k in seq_len(length(which_terminal))) {
      get_parent = as.numeric(as.character(tree$tree_matrix[which_terminal[k], 'parent']))
      get_split_var = as.character(tree$tree_matrix[get_parent, 'split_variable'])
      
      save_ancestor = rbind(save_ancestor,
                            cbind(terminal = which_terminal[k],
                                  ancestor = get_split_var))
      
      while (get_parent > 1) {
        get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent']))
        get_split_var = as.character(tree$tree_matrix[get_parent, 'split_variable'])
        save_ancestor = rbind(save_ancestor,
                              cbind(terminal = which_terminal[k],
                                    ancestor = get_split_var))
      }
    }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }
  
  return(save_ancestor)
}

# Get ancestors for internal nodes
get_ancestors_internal = function(tree) {
  save_ancestor = NULL
  tree = tree$tree_matrix
  which_internal = which(tree[,'terminal'] == 0)
  
  if(nrow(tree) == 1) {
    save_ancestor = cbind(internal = NULL, ancestor = NULL)
  } else {
    for (k in length(which_internal):1) {
      internal_node = which_internal[k]
      parent = tree[internal_node, 'parent']
      get_split_var = tree[internal_node, 'split_variable']
      
      save_ancestor = rbind(save_ancestor,
                            cbind(internal = internal_node,
                                  parent   = parent,
                                  split_var = get_split_var))
      
      while (!is.na(parent) && parent > 0) {
        get_split_var = tree[parent, 'split_variable']
        parent = tree[parent, 'parent']
        save_ancestor = rbind(save_ancestor,
                              cbind(internal = internal_node,
                                    parent   = parent,
                                    split_var = get_split_var))
      }
    }
  }
  return(save_ancestor[,,drop = FALSE])
}
