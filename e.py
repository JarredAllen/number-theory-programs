def is_e_prime(e):
    """Returns if e is a prime in E (2Z)"""
    return e % 4 == 2


def e_primes(e):
    """Generate the E-primes up to e"""
    return [i for i in range(e) if is_e_prime(i)]


def e_factor(e):
    """Returns a list of all prime factorizations of e in E (2Z)"""
    if is_e_prime(e):
        return [[e]]
    factors = []
    for e_prime in e_primes(e):
        if e_prime * e_prime > e:
            break
        if e % e_prime == 0:
            subfactors = e_factor(e // e_prime)
            for factor_list in subfactors:
                factors += [[e_prime] + factor_list]
    return [f for f in factors if f == sorted(f)]


def num_factorizations(e):
    """Returns the number of distinct factorizations of e in E (2Z)"""
    return len(e_factor(e))


def n_distinct_e_factors(n):
    """Find the smallest number which has at least `n` distinct
    factorizations in the E-primes, and print both the number and the
    factorizations."""
    i = 2
    while True:
        i_factors = e_factor(i)
        if len(i_factors) >= n:
            print(i)
            print(i_factors)
            break
        i += 2
