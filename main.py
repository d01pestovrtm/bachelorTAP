import SVGParser as parser
from detect_crossings import wirt_presentation
from sympy.combinatorics import Permutation
from find_generators import find_min_generating_set
import matrix as m
from find_homomorphism import perform_combinatoric
from sympy import SparseMatrix
import configparser
from sympy.matrices import Matrix, eye, zeros, ones
from sympy import Symbol

config = configparser.ConfigParser()
config.read('config.ini')
config_group_size = int(config['Settings']['group_size'])
config_group_type = config['Settings']['group_type']
config_filename = config['Settings']['input_file']
strands = parser.parse_svg(config_filename)
wirt_pres = wirt_presentation(strands)
indices = find_min_generating_set(wirt_pres)
print(indices)
print(wirt_pres)
print(len(wirt_pres))

t = Symbol('t')

# # def print_Jacobi_matrix(wirt_presentation):
# #     l = len(wirt_presentation)
# #     matrix = zeros(l)

# #     for r in range(l):
# #         sign = wirt_presentation[r][0]
# #         mu_idx = wirt_presentation[r][1]
# #         tau_idx = (r + 1) % l
# #         if sign < 0 :
# #             matrix[r ,  mu_idx] = -t**(-1) + 1# negative_der_mu(images[mu_idx], images[tau_idx])
# #             matrix[r ,  tau_idx ] = t**(-1)# negative_tau

# #         if sign > 0 : 
# #             matrix[r , mu_idx ] = 1 - t #positive_der_mu(images[r], size)
# #             matrix[r , tau_idx ] = t #positive_tau

# #         matrix[r, r] = -1
# #     print(matrix)


# # print_Jacobi_matrix(wirt_pres)


# # This one is only for AP
# # TODO move it to separate file
# def compute_Jacobi_matrix(wirt_presentation, relation_idx, generator_idx):
#     l = len(wirt_presentation)
#     matrix = zeros(l - 1)
#     relations_idx = [r for r in range(0, l) if r != relation_idx]
#     #print(matrix)
#     size = 1
#     for r in relations_idx:
#         sign = wirt_presentation[r][0]
#         mu_idx = wirt_presentation[r][1]
#         tau_idx = (r + 1) % l
#         #begin of row insertion
#         shifted_row = r - 1 if r > relation_idx else r
#         #for mu
#         if mu_idx != generator_idx:
#             shifted_column = mu_idx - 1 if mu_idx > generator_idx else mu_idx
#             if sign < 0 :
#                 matrix[shifted_row , shifted_column ] = -t**(-1) + 1# negative_der_mu(images[mu_idx], images[tau_idx])
#             if sign > 0 : 
#                 matrix[shifted_row , shifted_column ] = 1 - t #positive_der_mu(images[r], size)
        
#         #for tau
#         if tau_idx != generator_idx:
#             shifted_column = tau_idx - 1 if tau_idx > generator_idx else tau_idx 
#             #matrix[shifted_row : shifted_row  + size , shifted_column : shifted_column + size ] = der_tau(images[mu_idx], sign)
#             if sign < 0 :
#                 matrix[shifted_row , shifted_column ] = t**(-1)# negative_tau
#             if sign > 0 : 
#                 matrix[shifted_row , shifted_column ] = t #positive_tau

#         #for sigma
#         if r != generator_idx:
#             shifted_column = r - 1 if r > generator_idx else r
#             matrix[shifted_row, shifted_column] = -1#der_sigma(size)
#     return matrix
    
# m = compute_Jacobi_matrix(wirt_pres, 0, 0)
# print(m)
# print(m.det())


imgs = perform_combinatoric(config_group_size, config_group_type, wirt_pres, indices)
print(imgs)
for k in imgs.keys():
    print (m.permutation_to_matrix(imgs[k]))


l = len(imgs)
for i in range(l):
    for j in range(l):
        mat = m.compute_Jacobi_matrix(imgs, wirt_pres, i, j, config_group_size)
        #print(mat)
        det_upper = SparseMatrix(mat).det()
        det_lower = m.compute_denominator_matrix(imgs[0], config_group_size).det()
        print(f'for deleted {i}th row and {j}th column:')
        print((det_upper / det_lower) * t**(6))

