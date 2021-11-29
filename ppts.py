from math import gcd, sqrt


def is_int(f):
    return f - int(f) == 0


def generate_ppts(limit):
    """Create all PPTs with 2 <= a < b <= limit"""
    ppts = []
    for a in range(2, limit + 1):
        for b in range(a + 1, limit + 1):
            c = sqrt(a ** 2 + b ** 2)
            if is_int(c):
                c = int(c)
                if gcd(a, b, c) == 1:
                    ppts += [(a, b, c)]
    return ppts


def all_in_ppts(limit):
    """List all elements in PPTs with 2 <= a < b <= limit"""
    found = [False] * (limit + 1)
    for ppt in generate_ppts(limit):
        for n in ppt:
            if n > limit:
                continue
            found[n] = True
    return found


def find_pythagoreans_by_c(s_limit):
    """List all PPTs from t < s <= s_limit"""
    triples = dict()
    for s in range(3, s_limit + 1, 2):
        for t in range(1, s, 2):
            a = s * t
            b = (s ** 2 - t ** 2) // 2
            c = (s ** 2 + t ** 2) // 2
            if c in triples:
                triples[c] += [(a, b)]
            else:
                triples[c] = [(a, b)]
    return triples
