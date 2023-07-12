"""
On Length Spectrum of Translation Length(as a subproject of Translation Length of a Curve Graph of a Torus)
Researcher : Hyungryul Baik, Changsub Kim, Sanghoon Kwak.
Coded by Kwak on 03/14/2019.
Last Updated : 04/24/2019
"""

import math
import pdb
import time
import datetime

from Translation_Length_Calc import getTrlen

M = 1000 #Trace

def main(M):
    """
    1. Collect All Positive Matrices whose trace is less than or equal to M

    """

    time_start = time.clock()

    print "###Step 1 : Collect All Positive Matrices whose trace is less than or equal to %d" % (M)

    result = []
    for a in range(M):
        for d in range(M-a+1):
            if abs(a + d) <= 2:
                continue
            else:
                bc = a * d - 1
                for b in range(1,int(math.sqrt(abs(bc)))+1):
                    if bc % b == 0:
                        c = bc / b
                        result.append((a,b,c,d))
                        if b != c:
                            result.append((a,c,b,d))

    total = len(result)
    count = 0

    time_step1 = time.clock() - time_start
    print "###Step 1 Over; time elapsed : %f(s)" % (time_step1)
    
    """
    2. Calculate the translation length of each matrix, and classify them by the length.
    """

    print "###Step 2 : Calculate the translation length of each matrix, and classify them by the length."

    dicTrlen = dict()

    for mat in result:        
        trlen = getTrlen(mat[0],mat[1],mat[2],mat[3])
        try:
            dicTrlen[trlen].append(mat)
        except KeyError:
            dicTrlen[trlen] = [mat]
            count += 1
        if count % 50000 == 0:
            print "%s%%Done; %s out of %s matrices" % (int(float(count)/total * 100), count, total)

    time_step2 = time.clock() - time_start
    print "###Step 2 Over; time spent for step 2 : %f(s)" % (time_step2)

    """
    3. Print out the result.
    """

    print "*********************RESULT********************"
    print "Among %d positive matrices with the trace bounded above %d," % (total, M)
    print "1. The Length Spectrum : " + str(dicTrlen.keys())
    print "2. The Distribution of matrices by translation length"
    for trlen in dicTrlen.keys():
        print "The Number of matrices with Translation Length %d : %d" %(trlen, len(dicTrlen[trlen]))

    """
    4. Save it as a file.
    """
    
    today = datetime.date.today()
    filename = "trace_result(%d)_%s.txt" % (M, today.strftime("%B-%d-%Y"))
    file = open(filename, 'w')
    file.write("Time Elapsed : Step1 - %f, Step 2 - %f, Total - %f \n" % (time_step1, time_step2, time_step1+ time_step2))
    file.write("\n#############Summary############\n")
    file.write("Among %d matrices with a norm bounded above %d," % (total, M))
    file.write("\n1. The Length Spectrum : " + str(dicTrlen.keys()))
    file.write("\n2. The Distribution of matrices by translation length\n")
    for trlen in dicTrlen.keys():
        file.write("The Number of matrices with Translation Length %d : %d\n" %(trlen, len(dicTrlen[trlen])))

    file.write("\n############Matrices###########\n")
    for i in range(len(dicTrlen)):
        file.write("Length : %s \n" % str(i+1))
        for matrix in dicTrlen[i+1]:
            file.write(str(matrix)+",\n")

    file.close()
    """
    for l in dicTrlen.keys():
        dicTr = dict()
        for mat in dicTrlen[l]:
            tr = mat[0]+mat[3]
            try:
                dicTr[tr].append(mat)
            except KeyError:
                dicTr[tr] = [mat]
        for k in dicTr.keys():
            print "The Num of matrices with Trace %d and Translation Length %d : %d" %(k, l, len(dicTr[k]))
    """

main(50)
main(100)
main(150)
main(200)
main(300)
main(400)
main(500)
main(10000)
    
#print "3. The Equidistribution Ratio"
