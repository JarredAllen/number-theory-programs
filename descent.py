def descent(A, B, p, num_iters=None, print_output=True):
    """Perform Fermat's Descent, from two numbers A, B such that
    A**2 + B**2 = Mp for some multiple 1 < M < p.

    Repeat for the given number of times, or until you find A, B such that
    A**2 + B**2 = p.
    """

    def output(*args, **kwargs):
        if print_output:
            print(*args, **kwargs)

    if (A ** 2 + B ** 2) % p != 0:
        return
    iter = 0
    while A ** 2 + B ** 2 != p:
        if iter == num_iters:
            return A, B
        m = (A ** 2 + B ** 2) // p
        u = A % m
        if u >= m / 2:
            u -= m
        v = B % m
        if v >= m / 2:
            v -= m
        iter += 1
        oldA, oldB = A, B
        A, B = abs((u * A + v * B) // m), abs((v * A - u * B) // m)
        output(
            f"From ${oldA}^2 + {oldB}^2 = {p}M$, we solve and find that $M = {m}$. Reducing modulo $M$, we find $A \\equiv {u}, B \\equiv {v}$. Thus, we can take $A = {A}, B = {B}$ as a new solution equal to a smaller multiple of {p}.",
            end="\n\n",
        )
    output(f"Now, ${A}^2 + {B}^2 = {p}$, so we have found the solution.")
    return A, B


def find_sum_of_squares(n, print_output=True):
    """Find n as a sum of two squares, if such a solution exists"""

    def output(*args, **kwargs):
        if print_output:
            print(*args, **kwargs)

    x = None
    for i in range(n):
        if pow(i, 2, n) == n - 1:
            x = i
            output(
                "We find that the congruence $x^2+1 \\equiv 0 \\pmod{%d}$ has a solution at $x = %d$.\n"
                % (n, x)
            )
            break
    if x is None:
        output(
            "There is no solution to the congruence $x^2 + 1 \\equiv 0 \\pmod{%d}$, so no such solution exists."
            % n
        )
        return
    return descent(x, 1, n, print_output=print_output)
