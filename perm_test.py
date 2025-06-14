import permutation_utils as pu

S7 = pu.get_group_from_config(7, "dihedral")
S7 = pu.split_by_order(S7)
print(S7)
def run_trefoil():
    config_group_size = 3
    config_group_type = "symmetric"
    config_filename = "../knots/trefoil_polyline"
    strands = parse_svg(config_filename)
    wirt_pres = wirt_presentation(strands)
    indices = find_min_generating_set(wirt_pres)
    imgs = perform_combinatoric(config_group_size, config_group_type, wirt_pres, indices)
    print(imgs)

def run_f8():
    config_group_size = 5
    config_group_type = "dihedral"
    config_filename = "../knots/figure-eight_polyline"
    strands = parse_svg(config_filename)
    wirt_pres = wirt_presentation(strands)
    indices = find_min_generating_set(wirt_pres)
    imgs = perform_combinatoric(config_group_size, config_group_type, wirt_pres, indices)
    print(imgs)    

def run_Conway():
    config_group_size = 5
    config_group_type = "alternating"
    config_filename = "../knots/Conway_polyline"
    strands = parse_svg(config_filename)
    wirt_pres = wirt_presentation(strands)
    indices = find_min_generating_set(wirt_pres)
    imgs = perform_combinatoric(config_group_size, config_group_type, wirt_pres, indices)
    print(imgs)




def run_KT():
    config_group_size = 5
    config_group_type = "alternating"
    config_filename = "../knots/KT_polyline"
    strands = parse_svg(config_filename)
    wirt_pres = wirt_presentation(strands)
    indices = find_min_generating_set(wirt_pres)
    imgs = perform_combinatoric(config_group_size, config_group_type, wirt_pres, indices)
    print(imgs) 
