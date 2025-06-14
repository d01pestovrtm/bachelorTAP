from sympy.matrices import Matrix, eye, zeros, ones
from sympy import Symbol
def permutation_to_matrix(p):
    dim = p.size
    matrix = zeros(dim)
    for i in range(0, dim):
        matrix[ p(i), i] = 1
    return matrix

#TODO as transposed matrix
def inv_permutation_to_matrix(p):
    return permutation_to_matrix(p ** (-1))

def matrices_from_permutations(images):
    #tuple is better for performance, because we need inverse matrix also
    #[0] is inverse [1] is the normal
    return {i : (inv_permutation_to_matrix(images[i]), permutation_to_matrix(images[i])) for i in images.keys()}


#compute the Jacobi matrix without relation i and generator j
#images are permutations
#at the moment this works only for knots, not links
def compute_Jacobi_matrix_for_TAP(images, wirt_presentaion, relation_idx, generator_idx, size):
    l = len(images)
    matrix = zeros(size * (l - 1))
    relations_idx = [r for r in range(0, l) if r != relation_idx]
    for r in relations_idx:
        sign = wirt_presentaion[r][0]
        mu_idx = wirt_presentaion[r][1]
        tau_idx = (r + 1) % l
        #begin of row insertion
        shifted_row = (r - 1) * size if r > relation_idx else r * size
        #for mu
        if mu_idx != generator_idx:
            shifted_column = (mu_idx - 1) * size if mu_idx > generator_idx else mu_idx * size
            if sign < 0 :
                matrix[shifted_row : shifted_row  + size , shifted_column : shifted_column + size] = negative_der_mu(images[mu_idx], images[tau_idx])
            if sign > 0 : 
                matrix[shifted_row : shifted_row  + size , shifted_column : shifted_column + size] = positive_der_mu(images[r], size)
        
        #for tau
        if tau_idx != generator_idx:
            shifted_column = (tau_idx - 1) * size if tau_idx > generator_idx else tau_idx * size
            matrix[shifted_row : shifted_row  + size , shifted_column : shifted_column + size ] = der_tau(images[mu_idx], sign)

        #for sigma
        if r != generator_idx:
            shifted_column = (r - 1) * size if r > generator_idx else r * size
            matrix[shifted_row : shifted_row  + size , shifted_column : shifted_column + size ] = der_sigma(size)
    return matrix

def compute_Jacobi_matrix_for_AP(wirt_presentation, relation_idx, generator_idx):
      t = Symbol('t')
      l = len(wirt_presentation)
      matrix = zeros(l - 1)
      relations_idx = [r for r in range(0, l) if r != relation_idx]
      
      size = 1
      for r in relations_idx:
          sign = wirt_presentation[r][0]
          mu_idx = wirt_presentation[r][1]
          tau_idx = (r + 1) % l
          #begin of row insertion
          shifted_row = r - 1 if r > relation_idx else r
          #for mu
          if mu_idx != generator_idx:
              shifted_column = mu_idx - 1 if mu_idx > generator_idx else mu_idx
              if sign < 0 :
                  matrix[shifted_row , shifted_column ] = -t**(-1) + 1# negative_der_mu(images[mu_idx], images[tau_idx])
              if sign > 0 : 
                  matrix[shifted_row , shifted_column ] = 1 - t #positive_der_mu(images[r], size)
        
          #for tau
          if tau_idx != generator_idx:
              shifted_column = tau_idx - 1 if tau_idx > generator_idx else tau_idx 
              #matrix[shifted_row : shifted_row  + size , shifted_column : shifted_column + size ] = der_tau(images[mu_idx], sign)
              if sign < 0 :
                  matrix[shifted_row , shifted_column ] = t**(-1)# negative_tau
              if sign > 0 : 
                  matrix[shifted_row , shifted_column ] = t #positive_tau

          #for sigma
          if r != generator_idx:
              shifted_column = r - 1 if r > generator_idx else r
              matrix[shifted_row, shifted_column] = -1#der_sigma(size)
      return matrix


def compute_denominator_matrix(perm, size):
    t = Symbol('t')
    return ( permutation_to_matrix(perm)) * t - eye(size)

#TODO pass symbol earlier
#gen index check before
def positive_der_mu(perm, size):
    t = Symbol('t')
    return eye(size) - (permutation_to_matrix(perm)) * t


def negative_der_mu(perm_mu, perm_tau):
    t = Symbol('t')
    # M(-x_mu^-1 + x_mu^-1*x_tau) = -M(x_mu^-1)*t^-1 + M(x_mu^-1*x_tau)
    return -permutation_to_matrix(perm_mu ** (-1))*(t ** (-1)) + permutation_to_matrix(perm_tau * perm_mu ** (-1)) 


def der_tau(perm, sign):
    t = Symbol('t')
    return permutation_to_matrix(perm ** (sign)) * (t ** sign)

    
def der_sigma(size):
    return -eye(size) 

