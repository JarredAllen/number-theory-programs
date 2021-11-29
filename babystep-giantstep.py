import math

from euclid import extended_euclidean


def babystep_giantstep(p, g, h):
    """Solves the DLP: find $x$ such that $g^x \equiv h \pmod{p}$

    Uses Baby Step-Giant Step algorithm from Shanks"""
    n = 1 + math.isqrt(p)
    # Have a list and a set so we can quickly check containment
    # and we can also find what i value produced it
    baby_steps = [pow(g, i, p) for i in range(n + 1)]
    baby_steps_set = set(baby_steps)
    u = extended_euclidean(baby_steps[n], p)[0]
    probe = h
    for j in range(n + 1):
        if probe in baby_steps_set:
            i = baby_steps.index(probe)
            return i + n * j
        probe = (probe * u) % p
