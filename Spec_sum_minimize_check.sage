import numpy as np

def lambda1_plus_lambda2_symmetric(G):
    """
    lambda_1 + lambda_2 for adjacency matrix of G using symmetric eigensolver.
    Guaranteed real up to float rounding.
    """
    A = np.array(G.adjacency_matrix(), dtype=float)  # symmetric
    vals = np.linalg.eigvalsh(A)                     # ascending, real
    return float(vals[-1] + vals[-2])

def path_sum_closed_form(n):
    # For P_n, adjacency eigenvalues are 2cos(k*pi/(n+1)).
    return float(2*cos(pi/(n+1)) + 2*cos(2*pi/(n+1)))

def Y_tree(l1, l2, l3):
    """
    One degree-3 center; three arms of edge-lengths l1,l2,l3 >= 1.
    Total vertices = 1 + l1 + l2 + l3.
    """
    G = Graph()
    center = 0
    G.add_vertex(center)
    nxt = 1
    for L in (l1, l2, l3):
        prev = center
        for _ in range(L):
            v = nxt; nxt += 1
            G.add_edge(prev, v)
            prev = v
    return G

def double_star_leaf_subdiv(a, b, c, d):
    """
    Two degree-3 vertices u=0 and v=1 joined by an unsubdivided edge (0,1).
    Two arms off 0 of lengths a,b >=1 and two arms off 1 of lengths c,d >=1.
    Total vertices = 2 + a + b + c + d.
    """
    G = Graph()
    u, v = 0, 1
    G.add_vertices([u, v])
    G.add_edge(u, v)   # core edge NOT subdivided
    nxt = 2

    def attach_arm(root, L, nxt):
        prev = root
        for _ in range(L):
            w = nxt; nxt += 1
            G.add_edge(prev, w)
            prev = w
        return nxt

    nxt = attach_arm(u, a, nxt)
    nxt = attach_arm(u, b, nxt)
    nxt = attach_arm(v, c, nxt)
    nxt = attach_arm(v, d, nxt)
    return G

def min_over_Y_family(n):
    """
    Min of lambda1+lambda2 over Y-trees on n vertices.
    Parameterization: l1<=l2<=l3, l1+l2+l3 = n-1.
    """
    best = +infinity
    arg = None
    for l1 in range(1, n-2):
        for l2 in range(l1, n-1-l1):
            l3 = (n-1) - l1 - l2
            if l3 < l2 or l3 < 1:
                continue
            s = lambda1_plus_lambda2_symmetric(Y_tree(l1, l2, l3))
            if s < best:
                best, arg = s, (l1, l2, l3)
    return float(best), arg

def min_over_double_star_family(n):
    """
    Min of lambda1+lambda2 over the double-star-with-unsubdivided-core family.
    Parameterization: a,b,c,d >=1 with n = 2+a+b+c+d.
    Reduce duplicates by:
      - sorting arms on each side: a<=b, c<=d
      - ordering sides: (a,b) <= (c,d) lexicographically
    """
    S = n - 2
    best = +infinity
    arg = None

    for a in range(1, S-2):
        for b in range(a, S-1-a):
            R = S - a - b
            if R < 2:   # need c,d >=1
                continue
            for c in range(1, R):
                d = R - c
                if d < c or d < 1:
                    continue

                # canonicalize swapping the two degree-3 centers
                if (a, b) > (c, d):
                    continue

                s = lambda1_plus_lambda2_symmetric(double_star_leaf_subdiv(a, b, c, d))
                if s < best:
                    best, arg = s, (a, b, c, d)

    return float(best), arg

def check_n(n, tol=1e-9):
    path = path_sum_closed_form(n)

    bestY, argY = min_over_Y_family(n) if n >= 4 else (+infinity, None)
    bestD, argD = min_over_double_star_family(n) if n >= 6 else (+infinity, None)

    if bestY <= bestD:
        best, fam, arg = bestY, "Y", argY
    else:
        best, fam, arg = bestD, "DoubleStar", argD

    ok = (best + tol >= path)
    return ok, path, best, fam, arg, (best - path)

def run_upto(N=60, tol=1e-9):
    for n in range(15, N+1):
        ok, path, best, fam, arg, gap = check_n(n, tol=tol)
        print(f"n={n:2d} ok={ok}  path={path:.12f}  best={best:.12f}  "
              f"gap={gap:+.3e}  family={fam}  arg={arg}")
        if not ok:
            print("COUNTEREXAMPLE FOUND at n =", n)
            break

run_upto(60)
