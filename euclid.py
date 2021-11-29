def print_euclidean(a, b):
    while b != 0:
        q, r = divmod(a, b)
        print("{} &= {}\\times{} + {} \\\\".format(a, q, b, r), end="\n")
        a, b = b, r


def extended_euclidean(a, b, output=False):
    def display(*args, **kwargs):
        if output:
            print(*args, **kwargs)

    coeffs = [(0, 1), (1, 0)]
    og_a = a
    og_b = b
    while b != 0:
        q, r = divmod(a, b)
        display("{} &= {}\\times{} + {}".format(a, q, b, r), end="")
        display(" & " + "\\;" * 4 + " & ", end="")
        x = coeffs[-2][1] - q * coeffs[-1][1]
        y = coeffs[-2][0] - q * coeffs[-1][0]
        display("{} &= {}\\times{} &+ {}\\times{} \\\\".format(r, x, og_a, y, og_b))
        a, b = b, r
        coeffs += [(y, x)]
    x, y = coeffs[-2][::-1]
    if a * x + b * y == -1:
        x *= -1
        y *= -1
    return x, y


def print_lightningbolt(a, b):
    coeffs = []
    while b != 0:
        q, r = divmod(a, b)
        print("{} &= {}\\times{} + {} \\\\".format(a, q, b, r), end="\n")
        a, b = b, r
        coeffs.append(q)
    row1 = [0, 1]
    row2 = [1, 0]
    for coeff in coeffs:
        row1.append(row1[-2] + row1[-1] * coeff)
        row2.append(row2[-2] + row2[-1] * coeff)
    print("\\begin{tabular}{ c c|%s }" % " ".join(["c"] * len(coeffs)))
    print("\\hline")
    print("&&" + " & ".join(map(str, coeffs)) + "\\\\")
    print(" & ".join(map(str, row1)) + "\\\\")
    print(" & ".join(map(str, row2)))
    print("\\end{tabular}")
