from SVGParser import parse_svg
from detect_crossings import wirt_presentation
from find_generators import find_min_generating_set
from itertools import product
import permutation_utils


#TODO test this func
#TODO need to check only generating set ???
def verify_homomorphism(images, wirt_pres):
    """Verifies that the computed images under \phi: Wirt_Pres -> Permutation_group are group homomorphism.

    At the current state this function only checks if \phi(ab) = \phi(a) * \phi(b)
    Some smarter checks will be added later

    Keyword arguments:
    images -- dictionary showing to which permutation each strand/generator was sent
    wirt_pres -- dictionary representing the Wirtinger presentation
    """
    for i in range(0, len(wirt_pres)):
        sign = wirt_pres[i][0]
        mu_index = wirt_pres[i][1]
        mu = images[mu_index]
        sigma = images[i]
        tau = images[(i + 1) % len(wirt_pres)]
        Id = sigma ** (-1) * mu ** ((-1) * sign) * tau * mu ** (sign)

        if not Id.is_Identity:
            return False
    return True


def compute_hom(images, wirt_pres, order):

    def traversal(images, wirt_pres):
        for i in wirt_pres.keys():
            compute_presentation(i, images, wirt_pres)

    #TODO refactor
    def compute_presentation(i, images, wirt_pres):
        mu = images[ wirt_pres[i][1] ]
        sigma = images[i]
        tau = images[(i + 1) % len(wirt_pres)]
        sign = wirt_pres[i][0]
        if mu is not None:

            if (sigma == tau == None or 
                sigma is not None and tau is not None):
                return
            
            if sigma is not None:
                #relation for tau
                images[(i + 1) % len(wirt_pres)] = mu ** (sign) * sigma * mu ** ((-1) * sign)
            else:
                #relation for sigma
                images[i] = mu ** ((-1) * sign) * tau * mu ** (sign)

    while None in images.values():
        traversal(images, wirt_pres) 
    
    return images


def perform_combinatoric(config_group_size, config_group_type, wirtinger_pres, generating_set):
    Group = permutation_utils.get_group_from_config(config_group_size, config_group_type)
    group = list(Group)
    slices = permutation_utils.split_by_order(group)
    
    
    for k in slices.keys():
        slice_gr = slices[k]

        for comb in product(slice_gr, repeat = len(generating_set)):
            images = {i : None for i in range(0, len(wirtinger_pres))}
            for key, value in zip(generating_set, comb):
                images[key] = value

            #TODO maybe return None if permutations have different order
            images = compute_hom(images, wirtinger_pres, k)
    
            if (permutation_utils.is_epimorphism(images, config_group_size, config_group_type) and
                verify_homomorphism(images, wirtinger_pres)):
                #print(images)
                return images
    print("no epimorhism found")

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