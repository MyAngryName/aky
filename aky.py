#!/bin/python

from si import chinese_remainder as chrSatz
from si import extended_gcd as euklid
from si import EC, ECelem
from si import *
from ordnung import *
from math import ceil, sqrt

"""
for pow(n,d) change d, else keep it like it is
"""
def phi(n, d = 0):
    if (d == 0):
        primCheck = True
        ggtCheck = False
        for i in range(2,n):
            if((n % i) == 0):
                p = i
                q = (n/i)
                primCheck = False
                e = euklid(p, q, verbose=0)
                if ((e[0] == 1)):
                    ggtCheck = True
                    break
                else:
                    continue

        if (primCheck):
            return (n-1)
        elif (ggtCheck):
            return (p-1)*(q-1)
        else:
            primeFac = prime_factors(n,n)
            checker = []
            ph = n
            for i in range(len(primeFac)):
                if primeFac[i] in checker:
                    continue
                else:
                    ph *= (1 - (1/primeFac[i]))
                    checker.append(primeFac[i])
            return ph



    else:
        return (pow(n,(d-1))*(n-1))
"""
Checks if something is a Prime Number, if not prints out if it is a Carmichael Number
"""
def fermat_test(p):
    primes = list_primes(0,p)

    for i in primes:
        if(pow(i,p,p) == i):
            continue
        else:
            return False
    if (miller_rabin(p)):
        print("Fermat Test: True | Millar Rabin (40R): True")
        return True
    else:
        print("Fermat Test: True | Millar Rabin (40R): False = Carmichaelzahl")
        return True
"""
Bestimmte Anzahl an Runden
No x for 40 Rounds Test
x < 0 to Test how many rounds are needed
x > 0 for costum amount of rounds
"""
def x_miller(p, x = 0):
    if (x == 0):
        e = miller_rabin(p)
        return e
    elif(x < 0):
        e = True
        i = 1
        while(e):
            e = miller_rabin(p, max_rounds = i)
            i += 1
            return (i)
    else:
        e = miller_rabin(p, max_rounds = x)
        return e

"""

if 0 > x < 1 Give it a possiblity and get the amount of reounds needed
If x > 1 it returns the possiblity in Percent
"""
def miller_rounds_for_possibility(x):
    checker = 0
    i = 0
    if (x > 0 and x < 1):
        while True:
            checker = pow(0.25,i)
            if (checker <= x):
                return i
            else:
                i += 1
    elif (x > 1):
        return pow(0.25, x)
    else:
        return ("ERROR x < 0")
            

"""
Miller Rabin zur Basis a
n = prim
"""
def miller_rabin_basis(a, n, max_rounds = 40):
    """ use the Miller-Rabin primality test to check
    if n is a pseudoprime."""
    # factor n-1 into 2**s * k with k odd
    k = n-1
    s = 0
    while( k%2 == 0 ):
        k //= 2
        s += 1
    # test for at most max_rounds rounds
    for _ in range( max_rounds ):
        witness = True
        # choose random basis a
        a = pow( a, k, n )
        if a==1:
            # a**k == 1
            continue
        for r in range(s):
            if a==n-1:
                # a**(2**r)*k == -1
                # not a witness
                witness = False
                continue
            a = pow( a, 2, n )
        if witness:
            # witness found
            return False
    # max_rounds rounds survived
    perc = miller_rounds_for_possibility(max_rounds)
    return (True, perc) 


"""
g = Basis
h = Logarithmus von h
p = Primzahl Zp
"""
def bsgs(g, h, p):
    '''
    Solve for x in h = g^x mod p given a prime p.
    If p is not prime, you shouldn't use BSGS anyway.
    '''
    N = ceil(sqrt(p - 1))  # phi(p) is p-1 if p is prime

    # Store hashmap of g^{1...m} (mod p). Baby step.
    tbl = {pow(g, i, p): i for i in range(N)}

    # Precompute via Fermat's Little Theorem
    c = pow(g, N * (p - 2), p)

    # Search for an equivalence in the table. Giant step.
    for j in range(N):
        y = (h * pow(c, j, p)) % p
        if y in tbl:
            return j * N + tbl[y]

    # Solution not found
    return None
"""
g = Basis (x, y)
h = Logarithmus von (a, b)
e = EC()  - Elliptische Kurve im EC Format
ordnung = ordnung...
"""
def bsgs_ec(g ,h, e, ordnung):
    N = ceil(sqrt(ordnung))  # phi(p) is p-1 if p is prime
    q = ECelem(e, g)
    hEC = ECelem(e, h)
    list_baby = []
    for i in range(N):
        tmp = q.mult(i)
        list_baby.append(str(tmp))
    qInv = q.inverse()
    giantBasis = qInv.mult(N)
    
    for ii in range(N):
        tmp1 = giantBasis.mult(ii)
        tmp2 = tmp1.add(hEC)
        if (str(tmp2) in list_baby):
            for i in range(N):
                if(str(tmp2) == list_baby[i]):
                    return N * ii + i
    return None

"""
p = Prim
(r,s) = Signatur
g = g
A = Verifikationsschlüssel
h = hash Wert der Message
attack = True, wenn eine Verifikation probiert, wobei die erste Bedinung nicht geprüft werden soll.
"""
def elgamal_signatur_pruefung( r, s, h, g, A, p, attack = False):
    if not attack and not(1 <= r and r <= (p-1)):
        return False
    tmp1 = (pow(A,r,p) * pow(r,s,p)) % p
    tmp2 = pow(g, h, p)

    if not(tmp1 == tmp2):
        return False
    return True
"""
h1 = Bereits erhaltene Signatur
h2 = Hashwert für den die Signatur zu fälschen ist
"""
def elgamal_attack_signatur_pruefung(r,s,h1,h2,p):
    
    iMod = inverse_mod(h1, (p-1))
    if (iMod == None):
        return False
    u = h2 * iMod
    s1 = (s * u) % (p-1)

    r1 = chinese_remainder((p-1,p),(r*u,r))
    return (r1, s1)

def elgamal_attack_parameter(r, s1, s2, h1, h2, p):
    s = s1 - s2
    if (inverse_mod(s, p-1) == None):
        return False
    else:
        k = (inverse_mod(s, p-1) * (h1 - h2)) % (p-1)
    a = ((h1 - k * s1) % (p - 1))

    return a

def lineare_diophantische_gleichung(x, y, z):
    if (x > y):
        e = extended_gcd(x, y)
        xG = e[1]
        yG = e[2]
    elif(x < y):
        e = extended_gcd(y, x)    
        yG = e[1]
        xG = e[2]
    else:
        return None
    x1 = xG * z
    y1 = yG * z

    return ((x1 / e[0]), (y1 / e[0]), e[0])

def linear_diopahntisch_postiv(x,y,z):
    e = lineare_diophantische_gleichung(x,y,z)
    k = 0
    xk = e[0]
    yk = e[1]
    d = e[2]
    x /= d
    y /= d
    z /= d
    x1 = (xk - (y * k))
    while (x1 <= 0):
        k += 1
        x1 = (xk + (y * k))

    y1 = yk - (x * k)

    return (x1, y1)
