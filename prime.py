import math


def eratosthenes(n):
    """Perform the Sieve of Eratosthenes to identify all prime numbers
    up to n.

    Returns all such primes as a list.
    """
    sieve = [False, False] + [True] * (n - 1)
    for i in range(math.isqrt(n)):
        if not sieve[i]:
            continue
        for j in range(i * i, n + 1, i):
            sieve[j] = False
    return [n for (n, p) in enumerate(sieve) if p]


def is_prime(n):
    """Returns true iff n is prime"""
    n = abs(n)
    if n == 0 or n == 1:
        return False
    for p in eratosthenes(math.isqrt(n)):
        if n % p == 0:
            return False
    return True


def prime_factor(n):
    """Returns the prime factorization of n"""
    factors = []
    for p in eratosthenes(int(n ** 0.5 + 1)):
        while n % p == 0:
            n //= p
            factors += [p]
    if n > 1:
        factors += [n]
    return factors


def rabin_miller(n, iters=100):
    """Perform `iter` iterations on the Rabin Miller test, returning
    True iff n can still be a prime"""
    if iters > (n + 3) // 4:
        iters = (n + 3) // 4
    q = n - 1
    k = 0
    while q % 2 == 0:
        q //= 2
        k += 1
    for i in range(1, iters + 1):
        if i % n == 0:
            continue
        if not rabin_miller_test_case(n, q, k, i):
            return False
    return True


def rabin_miller_test_case(n, q, k, a):
    """Check the base a in the Rabin-Miller primality test, returning
    True iff a does not prevent $n = 1+q2^k$ from being prime"""
    base = pow(a, q, n)
    if base == 1 or base == n - 1:
        return True
    for _ in range(k):
        base = pow(base, 2, n)
        if base == n - 1:
            return True
        if base == 1:
            return False
    return False


def korselts(n, display=False):
    """Return true iff n passes Korselt's Criterion"""

    def output(s, *args, **kwargs):
        if display:
            print(s.format(*args, **kwargs))

    if n % 2 == 0:
        output("{} is even, so it is not Carmichael.", n)
        return False
    factors = prime_factor(n)
    if len(factors) == 1:
        output("{} is prime, so it is not Carmichael.", n)
        return False
    output("{} has prime factors: ${}$\n", n, factors)
    for factor in factors:
        if n % (factor ** 2) == 0:
            output(
                "${p}$ is a prime factor of ${n}$, and ${p}^2|n$, "
                "so $n$ cannot be Carmichael.",
                p=factor,
                n=n,
            )
            return False
        if (n - 1) % (factor - 1) != 0:
            output(
                "${p}$ is a prime factor of ${n}$, and ${p}-1\\nmid{n}-1$, "
                "so $n$ cannot be Carmichael.",
                p=factor,
                n=n,
            )
            return False
        output(
            "${p}$ is a prime factor of ${n}$, and both ${p}\\parallel{n}$ "
            "and ${p}-1\\mid{n}-1$.\n",
            p=factor,
            n=n,
        )
    output("All prime factors pass, so ${n}$ is Carmichael.", n=n)
    return True


def write_korselts_list(nums):
    """Given a list of numbers, output the LaTeX for finding whether or
    not they're Carmichael"""
    for n in nums:
        print("\\item ", end="")
        korselts(n, True)


def pollard(n, output=False):
    """Use Pollard's p-1 method to try to factor n (only for n=p*q)"""

    def display(*args, **kwargs):
        if output:
            print(*args, **kwargs)

    L = 2
    multiplier = 2
    while True:
        display(f"Let's try $L = {multiplier}! = {L}$: ", end="")
        probe = (pow(2, L, n) + n - 1) % n
        g = math.gcd(probe, n)
        display(
            f"We find $2^L - 1 \\equiv {probe} \\pmod{{{n}}}$, "
            f"and thus that $\\gcd({n}, {probe}) = {g}",
            end="",
        )
        if g != 1:
            display(f" \\neq 1$, so we have found the smooth prime $p = {g}$. ", end="")
            display(f"We can divide into $n$ to get $q = {n // g}$.\n")
            display(f"Thus, we can factor ${n} = {g} \\cdot {n // g}$.")
            return g, n // g
        display("$. This didn't work, so let's try the next value for $L$.\n")
        multiplier += 1
        L *= multiplier


def n_plus_b_squared(n, output=False):
    """Factor N by looking at N+b^2 values until a perfect square is reached"""

    def display(*args, **kwargs):
        if output:
            print(*args, **kwargs)

    b = 1
    while True:
        display(f"Let's try $b = {b}$: ", end="")
        nb2 = n + b ** 2
        display(f"We find $N+b^2 = {nb2}", end="")
        if math.isqrt(nb2) ** 2 == nb2:
            a = math.isqrt(nb2)
            display(f" = {a}^2$. Thus we find $N = (a+b)(a-b) = {a+b} \\cdot {a-b}$.\n")
            return a + b, a - b
        else:
            display("$, which is not a perfect square, so this didn't work.\n")
        b += 1


def k_n_plus_b_squared(k, n, output=False):
    """Factor N by looking at k*N+b^2 values until a perfect square is reached"""

    def display(*args, **kwargs):
        if output:
            print(*args, **kwargs)

    b = 1
    while True:
        display(f"Let's try $b = {b}$: ", end="")
        knb2 = k * n + b ** 2
        display(f"We find $k\\cdot N+b^2 = {knb2}", end="")
        if math.isqrt(knb2) ** 2 == knb2:
            a = math.isqrt(knb2)
            display(
                f" = {a}^2$. Thus we find $k \\cdot N = (a+b)(a-b) = {a+b}"
                f" \\cdot {a-b}$, ",
                end="",
            )
            f1 = math.gcd(a + b, n)
            f2 = math.gcd(a - b, n)
            display(f"so $N = {f1} \\cdot {f2}$.\n")
            return (f1, f2)
        else:
            display("$, which is not a perfect square, so this didn't work.\n")
        b += 1
