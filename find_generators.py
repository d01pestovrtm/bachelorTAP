from itertools import combinations


"""
This function finds minimal set of generators, such that we can compute other generators from them
"""
def find_min_generating_set(wirt_pres):
    for num_generators in range(2, len(wirt_pres)):
        res = find_solution(num_generators, wirt_pres)
        if res is not None:
            return res


"""
This function checks if there exists a generating set of size k and returns it
otherwise None
"""
def find_solution(num_generators, wirt_pres):
    def traversal(computed_generators, wirt_pres):
        updated = False
        for i in wirt_pres.keys():
            n = (i + 1) % len(wirt_pres) 
            if wirt_pres[i][1] in computed_generators:
                if i in computed_generators and n not in computed_generators:
                    updated = True
                    computed_generators.add(n)
                elif i not in computed_generators and n in computed_generators:
                    updated = True
                    computed_generators.add(i)
        return updated

    for generating_set in list(combinations(range(0, len(wirt_pres)), num_generators)):
        computed_generators = set(generating_set)
        while traversal(computed_generators, wirt_pres):
            pass
        if len(computed_generators) == len(wirt_pres):
            return generating_set
    return None
