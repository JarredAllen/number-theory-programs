def jacobi(a, b):
    """Evaluate the Jacobi symbol (a|b)"""
    a = a % b
    if a == 0:
        return 0
    elif a == 1:
        return 1
    elif a == -1:
        return 1 if b % 4 == 1 else -1
    elif a == 2:
        return 1 if b % 8 in [1, 7] else -1
    elif b == 2:
        return a % 2
    elif b == 3:
        return a == 1
    prod = 1
    j2b = jacobi(2, b)
    while a % 2 == 0:
        a //= 2
        prod *= j2b
        if a == 2:
            return prod * j2b
    flip_factor = -1 if a % 4 == 3 and b % 4 == 3 else 1
    return prod * jacobi(b, a) * flip_factor
