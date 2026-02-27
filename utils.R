## TODO: make these into functions

{ #Generate S matrix for Example 1
  # initialize
  S = matrix(0, nrow = n, ncol = m)
  
  # helper function to locate bottom indices
  get_idx = function(r) {
    base = (r - 1) * 3
    c(B = base + 1, D = base + 2, G = base + 3)
  }
  
  # --- row 1: total country ---
  for (r in 1:n_r) {
    idx = get_idx(r)
    S[1, idx["B"]] =  1
    S[1, idx["D"]] = -1
    S[1, idx["G"]] =  1
  }
  
  # --- rows 2:(1+n_r): each region ---
  for (r in 1:n_r) {
    idx = get_idx(r)
    row = 1 + r
    
    S[row, idx["B"]] =  1
    S[row, idx["D"]] = -1
    S[row, idx["G"]] =  1
  }
  
  # --- bottom level identity ---
  S[(1 + n_r + 1):n, ] = diag(m)
}


{# Generate Amat constraint matrix for Example 1
  # number of inequality constraints = 2 per region (B and D)
  n_constr = 2 * n_r
  
  Amat = matrix(0, nrow = m, ncol = n_constr)
  bvec = rep(0, n_constr)
  
  col_id = 1
  
  for (r in 1:n_r) {
    idx = get_idx(r)
    
    # B_r >= 0
    Amat[idx["B"], col_id] = 1
    col_id = col_id + 1
    
    # D_r >= 0
    Amat[idx["D"], col_id] = 1
    col_id = col_id + 1
  }
}
