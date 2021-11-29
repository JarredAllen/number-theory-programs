def m_seive_of_eratosthenes(n):
    """Performs the seive of eratosthenes on all M-numbers up to n
    to check for primes, returning a list of primes found"""
    prime_list = [False] + [True] * ((n - 1) // 4)
    for i in range(len(prime_list)):
        if not prime_list[i]:
            continue
        p = 4 * i + 1
        for j in range((p * p - 1) // 4, len(prime_list), p):
            prime_list[j] = False
    return [4 * i + 1 for i in range(len(prime_list)) if prime_list[i]]


MAX_M_PRIMES = 10 ** 6
m_primes = m_seive_of_eratosthenes(MAX_M_PRIMES)
m_primes_set = set(m_primes)


def is_m_prime(m):
    """Returns if m is a prime in M"""
    if m < MAX_M_PRIMES:
        return m in m_primes_set


def iter_m_primes(m):
    """Generate the E-primes up to e"""

    def first_larger_index(lst, n):
        for i in range(len(lst)):
            if lst[i] >= n:
                return i
        return len(lst)

    return m_primes[: first_larger_index(m_primes, m)]


def m_factor(m):
    """Returns a list of all prime factorizations of e in E (2Z)"""
    if is_m_prime(m):
        return [[m]]
    factors = []
    for m_prime in iter_m_primes(m):
        if m_prime * m_prime > m:
            break
        if m % m_prime == 0:
            subfactors = m_factor(m // m_prime)
            for factor_list in subfactors:
                factors += [[m_prime] + factor_list]
    return [f for f in factors if f == sorted(f)]


def num_factorizations(m):
    """Returns the number of distinct factorizations of m in M"""
    return len(m_factor(m))


def n_distinct_m_factors(n):
    """Find the smallest number which has at least `n` distinct
    factorizations in the M-primes, and print both the number and the
    factorizations."""
    i = 2
    while True:
        i_factors = m_factor(i)
        if len(i_factors) >= n:
            print(i)
            print(i_factors)
            break
        i += 2
