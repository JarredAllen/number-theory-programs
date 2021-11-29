def list_quadratic(n):
    hit = set()
    for a in range(n):
        hit.add(pow(a, 2, n))
    return sorted(hit)


def list_cubic(p, display_output=True):
    def output(*args, **kwargs):
        if display_output:
            print(*args, **kwargs)

    hit = set()
    output("\\begin{tabular}{ c|c }\n$a$ & $a^3$ \\\\\n\\hline")
    for a in range(1, p):
        hit.add(pow(a, 3, p))
        output("{} & {} \\\\".format(a, pow(a, 3, p)))
    output("\\end{tabular}")
    ans = sorted(hit)
    output(
        "Thus, we observe the cubic residues modulo {} are {}".format(
            p, ", ".join(map(str, ans))
        )
    )
    return ans
