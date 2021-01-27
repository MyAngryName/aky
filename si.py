"""
list_primes( a, b ):
    return a list of all primes p, st. a <= p < b
prime_factors(n,m): 
    compute all prime factors of n that are smaller than m (with multiplicity)
unique_prime_factors(n,m): 
    compute all prime factors of n that are smaller than m (without multiplicity)
equalmod(a, b, n):
    returns true, iff a = b mod n 
extended_gcd(a, b, verbose=0): 
    computes the extended GCD, i.e. gcd, x, y, st. gcd = ax + by
    shows all steps, when verbose>0
inverse_mod(a, m, verbose=0):
    computes the modular inverse of a mod m or None, if it does not exist
continued_fraction( a, b ):
    computes a continued fraction for a/b.
    returns a list
cf_approx( a ):
    returns a list of approximations for a continued fraction a
    approximations are pairs representing fractions
miller_rabin( n, max_rounds = 40 ):
    use the Miller-Rabin primality test to check, if n is a pseudoprime
EC(a,b,p):
    construct an elliptic curve y^2=x^3+ax+b mod p
    methods and attributes:
        order(): number of points on the curve (computes all points!)
        points(): list of all elements of the curve
        zero(): returns the neutral element of the group
ECelem(c,(x,y)):
    construct point (x,y) on the elliptic curve c
    tests, whether (x,y) is on c, but only warns, if it is not
    methods and attributes:
        p == q: checks, if two points are equal
        p.inverse(): compute additive inverse of point p
        p.add(q): add the points p and q (checks, if curves are identical)
        p.double(): double the point p
        p.mult(n): compute n*p for any integer n
        p.order(): compute the order of the point p
        
See examples at the end of the file.
"""

from functools import reduce
from math import gcd as math_gcd
from random import randint
from sympy import sieve

def list_primes( a, b ):
    """ return a list of all primes p, st. a <= p < b """
    return list( sieve.primerange( a, b ) )

def prime_factors( n, m ):
    """compute all prime factors of n that are smaller than m (with multiplicity)"""
    factors = []
    primes = list_primes(1,m)
    index = 0
    while n>1:
        if n%(primes[index]) == 0:
            factors.append(primes[index]) 
            n /= primes[index]
        else: 
            index += 1
            if index >= len(primes):
                break 
    return factors    

def unique_prime_factors( n, m ):
    """compute all prime factors of n that are smaller than m (without multiplicity)"""
    factors = []
    primes = list_primes(1,m)
    index = 0
    found = False
    while n>1:
        if n%(primes[index]) == 0:
            if not(found):
                factors.append(primes[index]) 
                found = True 
            n /= primes[index]
        else: 
            index += 1
            found = False
            if index >= len(primes):
                break 
    return factors
                
def equalmod( a, b, n ):
    """ returns true, iff a = b mod n """
    return (a-b)%n == 0

def extended_gcd( a, b, verbose=0 ):
    """ computes the extended GCD, i.e.
        gcd, x, y, st. gcd = ax + by """
    if a<b:
        return extended_gcd( b, a, verbose=verbose )
    max_width = max(len(str(a)),len(str(b)))+3
    max_width_plus_margin = 2*max_width-5
    x,y, u,v = 0,1, 1,0
    if verbose>0:
        print (max_width*" "+"{:{width_plus}d}{:{width_plus}d}".format(a, b, width_plus=max_width_plus_margin))
    while a != 0:
        q, r = b//a, b%a
        m, n = x-u*q, y-v*q
        b,a, x,y, u,v = a,r, u,v, m,n
        if verbose>0:
            if q==0:
                print ("{:{width}d}{:{width_plus}d}{:{width_plus}d}".format(b, x, y, width=max_width, width_plus=max_width_plus_margin))
            else:
                print ("{:{width}d}{:{width_plus}d}{:{width_plus}d}{:{width}d}".format(b, x, y, q, width=max_width, width_plus=max_width_plus_margin))
    gcd = b
    return gcd, x, y
    
def inverse_mod( a, m, verbose=0 ):
    """ computes the modular inverse of a mod m
        or None, if it does not exist """
    a = a % m
    gcd, x, y = extended_gcd(m, a, verbose)
    if gcd != 1:
        return None  # modular inverse does not exist
    else:
        return y % m

def continued_fraction( a, b ):
    """computes a continued fraction for a/b.
       returns a list
    """
    cf = []
    while b!=0:
       cf.append(a//b) 
       a, b = b, a%b
    return cf

def cf_approx( a ):
    """returns a list of approximations for a continued fraction a
       approximations are pairs representing fractions
    """
    p = [a[0],a[0]*a[1]+1]
    q = [1,a[1]]
    for n in range(2,len(a)):
        p.append(a[n]*p[n-1]+p[n-2])
        q.append(a[n]*q[n-1]+q[n-2])
    return list(zip(p,q))

def chinese_remainder( n, a ):
    """ compute the result of the Chinese Remainder Theorem
        for a list of moduli n and a list of remainders a
    """
    prod = 1
    for ni in n:
        if math_gcd( ni, prod ) > 1:
            raise ZeroDivisionError("Moduli have a common factor.")
        prod *= ni        
    sum = 0
    for n_i, a_i in zip(n, a):
        p = prod // n_i
        sum += a_i * inverse_mod(p, n_i) * p
    return sum % prod

def miller_rabin( n, max_rounds = 40 ):
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
        a = randint( 2, n-1 )
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
    return True
    
def recursive_order( pt, o ):            
        plist = unique_prime_factors(o,o+1)
        for p in plist:
            opt = pt.mult(o//p)
            if opt.iszero:
                return recursive_order(pt,o//p)
        return o 
        
class EC:
    """ class to represent an eliptic curve
        holds parameters prime, a and b """
    def __init__( self, a, b, prime ):
        self.a = a
        self.b = b
        self.prime = prime
    def __str__(self):
        return "y**2 = x**3 + "+str(self.a)+"*x + "+str(self.b)+" mod "+str(self.prime)
        
    def __repr__(self):
        return "EC("+str(self.a)+","+str(self.b)+","+str(self.prime)+")"
    
    def zero(self):
        return ECelem(self,None)
            
    def points(self):
        """compute a list of all points on the curve"""
        prime = self.prime
        #y2list = [(y**2)%prime for y in range(prime)]
        return [ECelem(self,(x,y)) for x in range(prime) for y in range(prime) if equalmod(y**2,x**3+self.a*x+self.b,prime)] + [ECelem(self,None)]
    def order(self):
        """compute the order of the group"""
        return len(self.points()) 

class ECelem:
    """ class for points on an EC
        ec: elliptic curve parameters as an object of class EC
        coord: a pair of integers or None (for the group identity) """
    def __init__( self, ec, coord ):
        if coord == None:
            self.iszero = True
            self.x = None
            self.y = None
        else: 
            self.iszero = False
            self.x, self.y = coord
            if not( equalmod( self.y**2, self.x**3+ec.a*self.x+ec.b, ec.prime ) ):
                print( "Warning: Point is not on the curve.")
                print( "Continuing ..." )
        self.ec = ec

    def __str__(self):
        if self.iszero:
#            return "zero on "+str(self.ec)
            return "infty"
        else: 
#            return "("+str(self.x)+","+str(self.y)+") on "+str(self.ec)
            return "("+str(self.x)+","+str(self.y)+")"
    
    def __repr__(self):
        if self.iszero:
            return "ECelem("+self.ec.__repr__()+", None)"
        else:
            return "ECelem("+self.ec.__repr__()+","+self.__str__()+")"

    def inverse(self):
        """ compute the inverse of a point """
        if self.iszero:
            return self
        else:
            return ECelem( self.ec, ( self.x, (-self.y)%self.ec.prime ) )
            
    def __neg__(self):
        """overload unary minus operator"""
        return self.inverse()
            
    def add( self, q, verbose=0 ):
        """ add two arbitrary points (may be the same point twice) """
        if verbose>0:
            print ('__________________________')
            print (str(self)+" + "+str(q)+":")
        if self.ec != q.ec:
            raise ("Cannot add points of two different elliptic curves")
        if self.iszero:
            if verbose>0:
                print (q)
                print ('==========================')
            return q
        if q.iszero:
            if verbose>0:
                print(self)
                print ('==========================')
            return self
        p = self.ec.prime
        x1, y1 = self.x, self.y
        x2, y2 = q.x, q.y
        if equalmod(x1,x2,p):
            if equalmod(y1,y2,p):
                return self.double(verbose)
            else: 
                return self.ec.zero()
                if verbose>0:
                    print (self.ec.zero())
                    print ('==========================')
        else: 
            k = ( (y1-y2)*inverse_mod(x1-x2,p,verbose-1) ) % p
            if k == None:
                raise("modular division by zero")
            x3 = (k**2 - x1 - x2) % p
            y3 = (-y1+k*(x1-x3)) % p
            if verbose>0:
                print ("     k = ("+str(y1)+"-"+str(y2)+") / ("+str(x1)+"-"+str(x2)+")")
                print ("       = "+str((y1-y2)%p)+" / "+str((x1-x2)%p) )
                print ("       = "+str((y1-y2)%p)+" * "+str(inverse_mod(x1-x2,p)))
                print ("       = "+str(k))
                print ("    x3 = "+str(k)+"**2 - "+str(x1)+" - "+str(x2))
                print ("       = "+str(x3))
                print ("    y3 = -"+str(y1)+" + "+str(k)+"*("+str(x1)+"-"+str(x3)+")")
                print ("       = -"+str(y1)+" + "+str(k)+"*"+str((x1-x3)%p))
                print ("       = "+str(y3))
                print (ECelem( self.ec, (x3, y3) ))
                print ('==========================')
            return ECelem( self.ec, (x3, y3) )

    def __add__( self, q ):
        """overload plus operator"""
        return self.add(q)

    def __sub__( self, q ):
        """overload binary minus operator"""
        return self.add(q.inverse())

    def double( self, verbose=0 ):
        """ double a point """
        if verbose>0:
            print ('__________________________')
            print (str(self)+" + "+str(self)+":")
        if self.iszero:
            if verbose>0:
                print(self)
                print ('==========================')
            return self
        p = self.ec.prime
        a = self.ec.a
        x1, y1 = self.x, self.y
        if equalmod(y1,0,p):
            if verbose>0:
                print (self.ec.zero())
                print ('==========================')
            return self.ec.zero()
        k = ( (3*x1**2+a)*inverse_mod(2*y1,p,verbose-1) ) % p
        x3 = (k**2 - 2*x1) % p
        y3 = (-y1+k*(x1-x3)) % p 
        if verbose>0:
            print ("     k = (3*"+str(x1)+"^2+"+str(a)+") / (2*"+str(y1)+")")
            print ("       = "+str((3*x1**2+a)%p)+" / "+str((2*y1)%p) )
            print ("       = "+str((3*x1**2+a)%p)+" * "+str(inverse_mod(2*y1,p)))
            print ("       = "+str(k))
            print ("    x3 = "+str(k)+"**2 - 2*"+str(x1))
            print ("       = "+str(x3))
            print ("    y3 = -"+str(y1)+" + "+str(k)+"*("+str(x1)+"-"+str(x3)+")")
            print ("       = -"+str(y1)+" + "+str(k)+"*"+str((x1-x3)%p))
            print ("       = "+str(y3))
        if verbose>0:
            print (ECelem( self.ec, (x3, y3) ))
            print ('==========================')
        return ECelem( self.ec, (x3, y3) )
        
    def mult( self, k, verbose=0 ):
        """ computes the result of adding the point k times with itself 
            using the binary method"""
        if verbose>0:
            print ("computing ",k,"*",self)
        if k<0:
            if verbose>0:
                print ("computing ",-k,"* (-",self,")")
            return self.inverse().mult(-k,verbose)
        s = self.ec.zero()
        base = self.add(s)
        while k!=0:
            if k%2==1:
                if verbose>0:
                    print ("adding ... ",end="")
                if verbose>1:
                    print ()
                s = s.add(base,verbose-1)
            if k>1:
                if verbose>0:
                    print ("doubling base ... ",end="")
                if verbose>1:
                    print ()
                base = base.double(verbose-1)
            k = k//2
        if verbose==1:
            print ("")
        return s

    def __rmul__( self, k ):
        """Overload the '*' operator"""
        return self.mult(k)

    def order(self):
        """compute the order of the point"""
        return recursive_order(self,self.ec.order())
        
    def __eq__( self, other ):
        return (self.x, self.y) == (other.x, other.y)        


if __name__ == "__main__":
    # ==============
    # == Examples ==
    # ==============
        
    # number theory
    print( list_primes(7,19) )
    print( prime_factors( 57480192000, 7 ) )
    print( unique_prime_factors( 57480192000, 100 ) )
    print( equalmod(199,99,10) )
    print( extended_gcd(1035,336) )
    print( extended_gcd(1035,336, verbose=1 ) )
    print( inverse_mod(17,33,verbose=1) )
    
    #define the curve y**2 = x**3 + 2x + 3 mod 541                    
    curve = EC(2,3,541)
    
    # compute all curve points
    cp = curve.points()
    print( cp[100] )
    # compute the order of the curve
    print ("curve oder is ",curve.order())
    
    # define points (3,6) and (7,3) on the curve
    q = ECelem(curve,(4,285))
    r = ECelem(curve, (93,151))
    
    print( -q )
    # define the zero of the curve
    z = curve.zero()
    # z = ECelem(curve,None) is equivalent
    
    # double q
    q2 = q.double()
    print( q+q )
    # verbose output
    q.double(1)
    # more verbose output
    q.double(2)
    
    # add q and r
    s = q.add(r)
    print( q+r )
    s = q+r
    # verbose
    s = q.add(r,2)
    
    # compute orders of points
    for pt in cp[18:24]:
        print (pt, " has order ",pt.order())
    
    # compute multiples of points
    # 9*r
    print( 9*r )
    r.mult(9,1)
    q.mult(3,2)
