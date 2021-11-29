import itertools
from math import gcd


def iterate_powers(g, p):
    """Iterate over the powers of g mod p"""
    i = 1
    while True:
        yield i
        i = (i * g) % p


def find_generator(p):
    """Find a generator for (Z/pZ)*, for any prime p"""
    for a in range(2, p):
        if len(set(itertools.islice(iterate_powers(a, p), p - 1))) == p - 1:
            return a
    return None


def find_all_generators(p):
    """Find a list of all generators for (Z/pZ)*, for any prime p"""
    g = find_generator(p)
    return [pow(g, a, p) for a in range(1, p - 1) if gcd(a, p - 1) == 1]


def discrete_log(g, p, a):
    """Solve g^h = a (mod p) for h."""
    return next(i for (i, pow) in enumerate(iterate_powers(g, p)) if pow == a)


def make_table(g, p, latex_output=False):
    """Output the LaTeX for a table of discrete logs"""
    targets = list(range(1, p))
    logs = [discrete_log(g, p, i) for i in targets]
    if latex_output:
        print("\\begin{tabular}{ c|c }")
        print("$a$ & $\\log_g(a)$ \\\\\n\\hline")
        for target, log in zip(targets, logs):
            print("{} & {} \\\\".format(target, log))
        print("\\end{tabular}")
    return dict(zip(targets, logs))


def check_diffs(g, p):
    """Check the differences of I(a) and I(b) in Z/pZ using g as a
    primitive root. Outputs all values which the difference cannot
    be."""
    table = make_table(g, p, latex_output=False)
    seen = [False for i in range(p - 1)]
    for a in range(1, p):
        b = (p + 1 - a) % p
        if b == 0:
            continue
        diff = (table[a] - table[b] + p - 1) % (p - 1)
        seen[diff] = True
    return [i for i, s in enumerate(seen) if not s]


def check_diffs_all_generators(p):
    """Returns the output of check_diffs for all generators"""
    return {g: check_diffs(g, p) for g in find_all_generators(p)}


def check_sums(g, p):
    """Check the sums of I(a) and I(b) in Z/pZ using g as a
    primitive root. Outputs all values which the difference cannot
    be."""
    table = make_table(g, p, latex_output=False)
    seen = [False for i in range(p - 1)]
    for a in range(1, p):
        b = (p + 1 - a) % p
        if b == 0:
            continue
        sum = (table[a] + table[b]) % (p - 1)
        seen[sum] = True
    return [i for i, s in enumerate(seen) if not s]


def check_sums_all_generators(p):
    """Returns the output of check_diffs for all generators"""
    return {g: check_sums(g, p) for g in find_all_generators(p)}
