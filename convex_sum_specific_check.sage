#Convex combination codes: Instructions, copy the approapriate code and run separately.

#Thm: 1.5: max \alpha = 0.5
n=18 # This can be run on personal computer upto 18.

G=graphs.trees(n)

M=0

for g in G:

    X=g.spectrum()

    if X[0]+X[1]>=M:

        M=X[0]+X[1]

        g.show()

        print(M)




# Minimization: n\le 18: Thm 1.8
n=16
G=graphs.trees(n)
M=4
for g in G:

    X=g.spectrum()

    if X[0]+X[1]<=M:

        M=X[0]+X[1]

        g.show()

        print(M)


# Minimization n = 18, max degree\ge 4. Thm 1.8
def maxdegree(G):

    D=G.degree_sequence()

    d=max(D)

    return d

n=18
G = graphs.trees(n)
for g in G: 

    if maxdegree(g)>=4:

        X=g.spectrum()

        if X[0]+X[1]<4:

            print(X[0]+X[1])

            g.show()
