from math import gcd


def find_order(a, m):
    """Compute the order of a mod m"""
    if gcd(a, m) != 1:
        return None
    i = 1
    power = a
    while power != 1:
        i += 1
        power = (power * a) % m
    return i
