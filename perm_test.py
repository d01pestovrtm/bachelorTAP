import permutation_utils as pu

S7 = pu.get_group_from_config(7, "dihedral")
S7 = pu.split_by_order(S7)
print(S7)

