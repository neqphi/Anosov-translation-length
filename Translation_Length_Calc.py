import math
import pdb

"""
This program calculates a translation length of a Anosov element on the Curve Graph of Torus(Farey graph)
Moreover, it specifies a segment of geodesic axis in Farey graph stabilized by given Anosov element,
in a form of extended rationals to represent a vertex of Farey graph.
#Programmed by
 Sanghoon Kwak, Mathematical Sciences Dept., KAIST (k_science_h@kaist.ac.kr)
#Based on the paper
 'On Translation Length of Curve Graph of Torus', by Hyungryul Baik, Changsub Kim, Sanghoon Kwak, Hyunshik Shin.
LAST UPDATE : Mar 19, 2019.
"""

class ExtRational:
    """
    Implements Extended Rational Number; i.e.
    1. Rational Numbers : p/q (p,q are relatively prime and q is positive)
    2. Infinity : 1/0(Displayed as "INF")

    *Property: numerator, denominator.
    *Method : 2 Getters for above properties.
    """
    def __init__(self, p, q):
        #Create a new ExtRational Instance
        #p = numerator, q = denominator
        if q == 0:
            if p > 0:
                p = 1; q = 0
            elif p < 0:
                p = -1; q = 0
            else:
                raise InputError # "Error: 0/0 is NOT a form of ExtRational!"
            
        else:
            g = gcd(p,q)
            p = p / g
            q = q / g
            if q < 0:
                p = -p; q = -q
        self.numerator = p
        self.denominator = q
        #sets value
        if self.denominator == 0:
            if self.numerator > 0:
                self.value = float("inf")
            elif self.numerator < 0:
                self.value = -float("inf")

        else:
            self.value = float(self.numerator)/self.denominator

    def __str__(self):
        #Specifies a string representation of ExtRational instance.
        if self.denominator == 0:
            if self.numerator > 0:
                return "INF"
            else:
                return "-INF"
        else:
            return str(self.numerator) + "/" + str(self.denominator)

    def __eq__(self, other):
        #Overrides the default equality.
        #Returns True if two ExtRationals represent the same value.
        if isinstance(other, ExtRational):
            return self.value == other.value
        return False

    def __ne__(self, other):
        #Overrides the default !=.
        #Functions the opposite of __eq__.
        return not self.__eq__(other)

    def getNumer(self):
        #Returns Numerator
        return self.numerator

    def getDenom(self):
        #Returns Denominator
        return self.denominator

    def getValue(self):
        #Returns its value
        return self.value

    def isIn(self, pq, rs):
        #Returns if 'self' is in between (as value) 'pq' and 'rs'.
        #Exclusive as endpoints.
        a = pq.getValue()
        b = rs.getValue()
        sV = self.getValue()
        return (a < sV < b) or (b < sV < a)

class Error(Exception):
    """Base class for exceptions in this module."""
    pass
    
class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = "Invalid input."
        self.message = "Invalid input"
    
def BezoutCoeff(a, b):
    #Input: (Int, Int), Output: (Int, Int)
    #Returns Bezout coefficient (x,y) of a,b; A solution of ax+by=1.
    #Ref: http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode

    s = 0; old_s = 1
    t = 1; old_t = 0
    r = b; old_r = a

    while r != 0:
        q = old_r / r
        (old_r, r) = (r, old_r - q * r)
        (old_s, s) = (s, old_s - q * s)
        (old_t, t) = (t, old_t - q * t)
    return (old_s, old_t)

def gcd(a,b):
    #Input: (Int, Int), Output: Positive Integer
    #Returns the greatest common divisor of a,b.
    #Euclidean Algorithm
    if a * b == 0:
        return abs(a + b)
    else:
        r = a % b
        return gcd(b,r)
    
def isNbd(pq,rs):
    p = pq.getNumer()
    q = pq.getDenom()
    r = rs.getNumer()
    s = rs.getDenom()
    return (abs(p * s - q * r) == 1)

def FareyCompare(pq,rs):
    # pq = rs is impossible; as they are not Farey nbd.
    if not isNbd(pq,rs):
        raise InputError # "Error: %s and %s are NOT nbd." % (pq, rs)

    p = abs(pq.getNumer()) # abs: To manange negative numbers. bigger absVal -> bigger tier.
    q = pq.getDenom()
    r = abs(rs.getNumer())
    s = rs.getDenom()
    
    if (q > s) or (q == s and p > r):
        return 1
    else:
        return -1

def samePair(p1,p2):
    #Returns True if two pairs p1, p2 are the same as sets
    return ((p1[0] == p2[0] and p1[1] == p2[1])) or ((p1[0] == p2[1]) and (p1[1] == p2[0]))

def FareySort(pq,rs):
    #Decreasing Order
    if FareyCompare(pq,rs) == 1:
        return pq, rs
    else:
        return rs, pq
    
def FareySum(pq,rs):
    if not isNbd(pq,rs):
        raise InputError # "Error: %s and %s are NOT nbd." % (pq, rs)

    p = pq.getNumer()
    q = pq.getDenom()
    r = rs.getNumer()
    s = rs.getDenom()
    result = ExtRational(p+r,q+s)
    return result

def FareySub(pq,rs):
    if not isNbd(pq,rs):
        raise InputError # "Error: %s and %s are NOT nbd." % (pq, rs)
    pq, rs = FareySort(pq, rs)
    p = pq.getNumer()
    q = pq.getDenom()
    r = rs.getNumer()
    s = rs.getDenom()
    return ExtRational(p-r,q-s)
   
def FareyNext(pq,rs,direction):
    if direction == 1:
        return FareySum(pq,rs)
    elif direction == -1:
        return FareySub(pq,rs)
    else:
        raise InputError
    
def ancestorOf(pq):
    #Input: ExtRational, Output: ExtRational
    #Returns the First Parent "\alpha" of pq. (The Map defined in Section 4 of [1].)

    q = pq.getDenom()
    p = pq.getNumer()

    if q == 0:
        return None #INF has no parents

    nbds =[] #Farey Neighborhoods of pq which have a smaller denominator than pq.
    for s in range(q):
        if (- 1 + p*s) % q == 0: #ps - qr =  1
            r = (- 1 + p*s) // q
            nbds.append((r,s))
        if (1 + p*s) % q == 0: #ps - qr = -1
            r = (1 + p*s) // q
            nbds.append((r,s))
    # Find the one with minimal denominator.
    # If two ExtRationals have the same denominator, then choose one with smaller numerator.
    nbds.sort(key=lambda element: (element[1], element[0]))
    result = ExtRational(nbds[0][0],nbds[0][1])
    return result

def ancestorPath2Inf(pq):
    #Input: ExtRational, Output: List of ExtRationals
    #Returns a geodesic path from pq to INF.
    path = [pq]
    while ancestorOf(pq) != None:
        pq = ancestorOf(pq)
        path.append(pq)
    return path

def mobius(a,b,c,d,pq):
    #Input: Int,Int,Int,Int,ExtRational, Output: ExtRational
    #Returns the image of Mobius Transformation (a,b,c,d) of given ExtRational pq :
    p = pq.getNumer()
    q = pq.getDenom()
    return ExtRational(a*p + b*q, c*p + d*q)

def ancestorPath(pq,rs):
    #Input: ExtRational, Output: List of ExtRationals
    #Returns a geodesic path from pq to rs.
    r = rs.getNumer()
    s = rs.getDenom()
    if s == 0:
        return ancestorPath2Inf(pq)

    #Find a Mobius Transformation (a,b,c,d) mapping rs to INF.
    (c, d) = (-s, r)
    (a, b) = BezoutCoeff(r, s)

    #Transform pq by the above transformation. Say it (pq)'
    pq = mobius(a, b, c, d, pq)

    #Take the inverse Mobius transformation of (a,b,c,d) on the path from (pq)' to INF.
    result = [mobius(d, -b, -c, a, xy) for xy in ancestorPath2Inf(pq)]
    return result

"""
LADDER GENERATING MODULE
"""

def genLadder(pq,rs,xy,zw):

    if samePair((pq,rs),(xy,zw)):
        raise InputError

    pivotList = []
    typeList = []
    prev = None
    mn = FareySum(xy,zw) #This helps to figure out the relative position of (xy,zw) to (pq,rs), without considering the xy == pq or zw == rs cases.
    while True:
        ###1. New Point 
        if mn.isIn(pq,rs):
            tu = FareySum(pq,rs)
        else:
            tu = FareySub(pq,rs)

        ###2. Rung Choice
        (d1, d2) = FareySort(pq,rs) #d1 =q

        if mn.isIn(d1,tu): 
            choice, nonchoice = d1, d2
        else:
            choice, nonchoice = d2, d1

        ###3. Update Lists
        if choice == prev:
            typeList[-1] += 1
        else:
            if prev == None:
                pivotList.append(nonchoice)
            pivotList.append(choice)
            typeList.append(1)
            prev = choice

        ###4. Check Point
        if samePair((choice,tu),(xy,zw)):
            pivotList.append(tu)
            return pivotList, typeList

        ###5. Redefine pq, rs
        (pq,rs) = (choice, tu)
    
def calibrateLadder(pivotList, ladderType, mn):
    #We calibrate ladder odd length ladder > 1.
    if len(ladderType) % 2 == 0:
        return pivotList, ladderType

    intTranslate = ladderType[0]
    ladderType = ladderType[1:] #well defined; ladder length >= 3
    lastElt = pivotList.pop()
    lastPivot = pivotList[-1] #well defined; pivot list >= 5.
    
    for i in range(intTranslate):
        if mn.isIn(lastPivot,lastElt): #well defined ladder length >= 5.
            lastElt = FareySum(lastPivot,lastElt)
        else:
            lastElt = FareySub(lastPivot,lastElt)

    pivotList = pivotList[1:]
    pivotList.append(lastElt)
    ladderType[-1] += intTranslate

    return pivotList, ladderType

"""
EFFICIENT GEODESIC FINDING MODULE
"""
def ignoreLastDigit(ladderType):
    #To ignore last digit(Make it Infinite)
    #Force the last move to be "t"
    return ladderType[:-1] + [99]

def alterLadder(ladderType):
    #start from top; cyclically translate initial digit to the last digit
    return ladderType[1:len(ladderType)] + [ladderType[0]]

def concatLadder(ladderType,num):
    result = ladderType
    for i in range(num-1):
        result = result + ladderType
    return result

def getGeodCode(ladderType):
    n = len(ladderType)
    ladder1 = ignoreLastDigit(ladderType) #path starting from bottom
    ladder2 = ignoreLastDigit(alterLadder(ladderType)) #path starting from top
    ladderSet = [ladder1, ladder2]
    strSet = []
    for i in range(2):
        ladder = ladderSet[i]
        strGeodesic = ""
        idx = 0
        while idx < n:
            if ladder[idx] == 1:
                strGeodesic += "p"
                idx += 2
            else:
                strGeodesic += "t"
                idx += 1
        strSet.append(strGeodesic)
    if len(strSet[0]) <= len(strSet[1]):
        return strSet[0] + "O" #Original Path Tag (Starting from first digit)
    else:
        return strSet[1] + "A" #Alternative Path Tag (Starting from second digit)

def getGeodesic(pivotList, geodCode):
    result = []
    if geodCode[-1] == "O":
        idx = 0
    elif geodCode[-1] == "A":
        idx = 1

    result.append(pivotList[idx])
                  
    for move in range(len(geodCode)-1): #-1 to exclude Terminal Path Tag
        if geodCode[move] == "t":
            idx += 1
        elif geodCode[move] == "p":
            idx += 2
        result.append(pivotList[idx])

    return result
    
"""
MAIN
"""

def getPivotsWithType(a,b,c,d):
    """
    Returns the efficient geodesic associated to given PSL(2,Z) element.
    """
    Det = a * d - b * c
    Tr = a + d
    if Det != 1:
        print "Error: (%d,%d,%d,%d) is not a PSL(2,Z) element: Determinant = %d" % (a,b,c,d,Det)
        return
    else:
        if Tr ** 2 <= 4:
            print "Error: (%d,%d,%d,%d) is not an Anosov element: Trace = %d" % (a,b,c,d,Tr)
            return
                        
    #Fixed Points
    fp1 = (a - d - math.sqrt((a+d)**2 - 4)) / (2 * c)
    fp2 = (a - d + math.sqrt((a+d)**2 - 4)) / (2 * c)
    mid = ExtRational((a - d), (2 * c))                              

    if fp1 * fp2 < 0:
        ladderPoint = ExtRational(0,1)
        ladderPoint2 = ExtRational(1,0)
    else:
        while True:
            if not (fp1 < ancestorOf(mid).getValue() < fp2):
                ladderPoint = mid
                ladderPoint2 = ancestorOf(mid)
                break
            else:                
                mid = ancestorOf(mid)

    temp = genLadder(ladderPoint,ladderPoint2, mobius(a,b,c,d,ladderPoint), mobius(a,b,c,d,ladderPoint2))
    mn = FareySum(mobius(a,b,c,d, temp[0][1]),mobius(a,b,c,d, temp[0][2]))
    
    temp = calibrateLadder(temp[0],temp[1],mn)
    (pivotList,ladderType) = (temp[0],temp[1])
    return (pivotList, ladderType)

def getEffGeod(pivotList, ladderType):
    effGeodesic = getGeodesic(pivotList, getGeodCode(ladderType))
    return effGeodesic

def getTrlen(a,b,c,d):
    """
    Returns the translation length of given PSL(2,Z) element.
    """
    (pivotList, ladderType) = getPivotsWithType(a,b,c,d)
    
    effGeodesic = getEffGeod(pivotList, ladderType)
    return len(effGeodesic)-1

def main():
    temp = raw_input("\nType the 4 entries of an Anosov element \n \t\t\t|a  b|\n \t\t\t|c  d|\n in PSL(2,Z) in the following form: a b c d: \t")

    tempList = [int(i) for i in temp.split(" ")]
    (a,b,c,d) = (tempList[0],tempList[1],tempList[2],tempList[3])

    pivotList, ladderType = getPivotsWithType(a,b,c,d)
    effGeodesic = getEffGeod(pivotList, ladderType)
    trlen = len(effGeodesic)-1
    
    imgGeod = [mobius(a,b,c,d,ext) for ext in effGeodesic]
    invimgGeod = [mobius(d,-b,-c,a,ext) for ext in effGeodesic]
    
    print "\n <####################RESULT#####################>\n"
    print "0. INPUT"
    print "\t\t\t|%4d %4d|" % (a,b)
    print "\t\t\t|%4d %4d|" % (c,d)
    print ""
    print "1. TRANSLATION LENGTH: \t",
    print trlen
    print ""
    print "2. INVARIANT LADDER"
    print " 1) Pivot List:\t",
    for ext in pivotList:
        print ext,
    print ""
    print " 2) Ladder Type: ",
    for tp in ladderType:
        print tp,
    print "\n"

    print "3. INVARIANT GEODESICS(EFFICIENT GEODESICS)"
    print " 1) Invariant Geodesic:\t\t",
    for ext in effGeodesic:
        print ext,
    print ""
    
    print " 2) Image:\t \t \t \t",
    for ext in imgGeod:
        print ext,
    print ""
    
    print " 3) Inverse Image:   ",
    for ext in invimgGeod:
        print ext,
    print ""
    
    print " 4) Concatenation[3) + 1) + 2)]: \t",
    print "~ ~ ~",
    for ext in invimgGeod[:-1] + effGeodesic[:-1] + imgGeod:
        print ext,
    print "~ ~ ~"
    print ""
    
    print "<#####################DONE######################>"
    

while(True):
    main()
    prompt = raw_input("Try another inputs?(y/n)")
    if (prompt == "n"):
        break
