from sympy.combinatorics import Permutation, Cycle
from sympy.combinatorics.perm_groups import PermutationGroup

from sympy.combinatorics.generators import symmetric
from sympy.combinatorics.generators import alternating
from sympy.combinatorics.generators import dihedral

def split_by_order(group):
        res = dict()
        for g in group:
                order = g.order()
                if order not in res.keys():
                        res[order] = []
                res[order].append(g)
        #remove identity
        res.pop(1, None)
        return res

def get_group_from_config(config_group_size, config_group_type):
    if config_group_type == "symmetric":
        return symmetric(config_group_size)
    if config_group_type == "alternating":
        return alternating(config_group_size)
    if config_group_type == "dihedral":
        return dihedral(config_group_size)
    
def is_epimorphism(images, config_group_size, config_group_type):
    #p = [Cycle(p) if not p.is_Identity else Cycle(config_group_size) for p in images.values()]
    p = []
    #TODO refactor
    for v in images.values():
        if not v.is_Identity:
            p.append(Cycle(v))
    if len(p) == 0:
        return False
    Group = PermutationGroup(*p)
    #TODO switch case 
    if config_group_type == "symmetric":
        first, second = Sn_generators(config_group_size)
    
    if config_group_type == "alternating":
        first, second = An_generators(config_group_size)

    if config_group_type == "dihedral":
        first, second = Dn_generators(config_group_size)
    
    return Group.contains(first) and Group.contains(second)

def Sn_generators(n):
    if n == 1:
        return None, None
    #transposition (1 2)
    first = [i for i in range(0, n)]
    first[0], first [1] = 1, 0

    #cycle (1 2 3 ... n)
    second = [(i + 1) % n for i in range(0, n)]
    return Permutation(first), Permutation(second)

def An_generators(n):
    if n <= 2:
        return Sn_generators(n)
    
    #cycle (1 2 3)
    first = [i for i in range(0, n)]
    first[0], first[1], first[2] = 1, 2, 0

    #cycle (1 2 ... n) for odd n or (2 3 ... n) for even n
    second = [(i + 1) % n for i in range(0, n)]
    if n % 2 == 0:        
        second[0], second[-1] = 0, 1
    return Permutation(first), Permutation(second)

def Dn_generators(n):
    if n <= 2:
        return None, None
    
    #cycle (1 2 3 ... n)
    first = [(i + 1) % n for i in range(0, n)]

    #mirror imaging
    second = [i  for i in range(0, n)]
    if n % 2 == 0:
        for i in range(1, n  // 2 ):
            second[i], second[-i] = second[-i], second[i]
    else:
        for i in range(1, (n + 1) // 2):
            second[i], second[-i] = second[-i], second[i]
    return Permutation(first), Permutation(second)
