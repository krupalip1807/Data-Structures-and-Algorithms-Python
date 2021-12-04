import sys
sys.path.append('/u1/h0/jkinne/public_html/cs458-f2020/code//GENERIC')
from all_subsets import *

# ./p1.py --x TTCGAAG --y TCAG --alg lcs_naive
# ./p1.py --x GGAAATCCGC --y CATGAACC --alg lcs_naive
def lcs_naive(X, Y, alg):
    
    Lx = first_subset()
    longest = ''
    longest_Lx = []
    longest_Ly = []
    verbose = alg.verbose
    while len(Lx) < len(X):
        Lx = next_subset(Lx, len(X))

        Ly = first_subset()
        while len(Ly) < len(Y):
            Ly = next_subset(Ly, len(Y))

            # Lx and Ly
            s_x = string_subset(Lx, X)
            s_y = string_subset(Ly, Y)
            if s_x == s_y and len(s_x) > len(longest):
                longest = s_x
                longest_Lx = Lx
                longest_Ly = Ly

    #print(longest, longest_Lx, longest_Ly)
        
    if verbose>=0:
        print(len(longest))
    if verbose>=1:
        max = 0
        seqX = ''
        seqY = ''
        if len(X)>len(Y):
            max = len(X)
        else:
            max = len(Y)
        #print X
        clear = 0
        for i in range(0, max):
            if(len(X)>i):
                if(i==longest_Lx[clear]):
                    if(longest_Lx[clear]>=longest_Ly[clear]):
                        seqX += X[i]
                    else:
                        dif = longest_Ly[clear] - longest_Lx[clear]
                        sp = ''
                        for b in range(0,dif):
                            sp += ' '
                        seqX += (sp + X[i])                                    
                    clear +=1
                else:
                    seqX += X[i]
        print(seqX)
        #print |
        #print Y
        clear = 0
        for i in range(0, max):
            if(len(Y)>i):
                if(i==longest_Ly[clear]):
                    if(longest_Ly[clear]>=longest_Lx[clear]):
                        seqY += Y[i]
                    else:
                        dif = longest_Lx[clear] - longest_Ly[clear]
                        sp = ''
                        for b in range(0,dif):
                            sp += ' '
                        seqY += (sp + Y[i])                                    
                    clear +=1
                else:
                    seqY += Y[i]
        print(seqY)
        
    #if verbose>=2:
         
        

    
# ./p1.py --x_file /u1/junk/kinne/genomes/gene_transcripts/GRCh38_latest_rna.fna.100k --y "GACCCCAAAATCAGCGAAAT" --alg lcs_dp    
def lcs_dp(X, Y, alg):
    verbose = alg.verbose
    m = len(X)
    n = len(Y)
    c = [[0]*(n+1) for i in range(0,m+1)]
    b = [['']*(n+1) for i in range(0,m+1)]

    for i in range(1, m+1):
        for j in range(1, n+1):
            if X[i-1] == Y[j-1]:
                c[i][j] = c[i-1][j-1] + 1
                b[i][j] = '\\'
            elif c[i-1][j] >= c[i][j-1]:
                c[i][j] = c[i-1][j]
                b[i][j] = '|'
            else:
                c[i][j] = c[i][j-1]
                b[i][j] = '-'
    result = {'c': c, 'b': b}

    print(X)
    print(Y)
    if verbose>=0:
        print(result['c'][len(X)][len(Y)])
       # print(result)
    if verbose>=1:
        i=1
        k=0
        sequence=''
        matchX = [None]*len(X)
        matchY = [None]*len(Y)
        while i<=len(X):
             j = 1
             while j<=len(Y):
                  if b[i][j]=='\\':
                      #print(i, j)
                      sequence += X[i-1]
                      matchX[k] = i-1
                      matchY[k] = j-1
                  j = j+1
             i = i+1
        cnt = 0
        j = 0
        fX = ''
        fY = ''

        for s in range(len(matchX)):
             while j<len(X):
                   sp = ''
                   if j == matchX[s]:
                        if(matchX[s]<matchY[s]):
                              for b in range(0,(matchY[s]-matchX[s])):
                                    sp = sp + ' '
                   fX = fX + sp + X[j]
                   j = j+1

        j = 0
        for s in range(len(matchY)):
             while j<len(Y):
                   sp = ''
                   if j == matchY[s]:
                        if(matchY[s]<matchX[s]):
                              for b in range(0,(matchX[s]-matchY[s])):
                                    sp = sp + ' '
                   fY = fY + sp + Y[j]
                   j = j+1

        #print(sequence)
        print(fX)
        print(fY)

    if verbose>=2:
        print('', list(Y))
        for r in range(len(c)):
            if r==0:
                print('',c[r])
            else:
                print(X[r-1], c[r])
            #..print(row)
        print(result) 


def edit_dp(X, Y, verbose):
    #print('algorithm not yet implemented')
    subcost = 0
    sizeX = len(X)+1
    sizeY = len(Y)+1
    mat = [[0]*(sizeX) for i in range(0,sizeY)]     #declare int d[0..m, 0..n]  set each element in d to zero
    for i in range(sizeX):
        mat[i][0] = i                               #for i from 1 to m: d[i, 0] := i
    for j in range(sizeY):
        mat[0][j] = j                               #for j from 1 to n: d[0, j] := j

    for j in range(sizeY):                          #for j from 1 to n:
        for i in range(sizeX):                      #for i from 1 to m:
            if X[i] == Y[j]:                        #if s[i] = t[j]:
                subcost = 0                         #substitutionCost := 0
            else:
                subcost = 1                         #substitutionCost := 1
            
            mat[i][j] = min(mat[i-1][j]+1,          #d[i, j] := minimum(d[i-1, j] + 1,                   // deletion
                            mat[i][j-1]+1,          #                   d[i, j-1] + 1,                   // insertion
                            mat[i-1][j-1]+subcost)  #                   d[i-1, j-1] + substitutionCost)  // substitution
    print(mat)                                      #return d[m, n]        

def match_naive(X, Y, verbose):
    txt = X
    pat = Y
    M = len(txt) 
    N = len(pat) 
    pos = -1
    found = 0
    brkchain = 0
  
    print("1")
    for i in range(0, (M-N)):
        print("2")
        for j in range (0, N):
            print("3")
            if X[i] == Y[j]:
                print("4")
                if pos<0:
                    pos = i
                found = found+1
                break
            elif found>0:
                print("5")
                brkchain = 1
                break
            j += 1
        if(brkchain):
            print("6")
            break
  
    if found==N:
        print("Found at index " + str(pos))
    else:
        print("Not found")



def match_RK(X, Y, verbose):
    A = len(X)
    B = len(Y)
    t = 0    #hash of X
    p = 0    #hash of Y
    h = 1
    d = 256
    q = 101
 
    for i in range(B-1):
        h = (h*d)%q

    for i in range(B):
        t = (d*t + ord(X[i]))%q
        p = (d*p + ord(Y[i]))%q
 
    for i in range(A-B+1):
        if p==t:
            for j in range(B):
                if X[i+j] != Y[j]:
                    break
                else: j+=1
            if j==B:
                print("Found at index " + str(i))
        if i < A-B:
            t = (d*(t-ord(X[i])*h) + ord(X[i+B]))%q
            if t < 0:
                t = t+q



algorithms = {'lcs_naive': lcs_naive,
              'lcs_dp': lcs_dp,
              'edit_dp': edit_dp,
              'align_NW': align_NW,
              'align_H': align_H,
              'align_SW': align_SW,
              'blast': blast,
              'match_naive': match_naive,
              'match_RK': match_RK}

help = '''lcs_naive: computes LCS using brute force (from class).
lcs_dp: uses dynamic programming (algorithm from the text, code from class).  
edit_dp: computes the edit distance using dynamic programming (as described in wikipedia, with unit-cost operations insert, delete, substitute, transpose).  
align_NW: aligns strings using the Needleman–Wunsch algorithm as described on wikipedia.  
align_H: aligns strings using the Hirschberg algorithm as described on wikipedia.
align_SW: aligns strings using the Smith–Waterman algorithms as described on wikipedia.  
blast: aligns strings using the BLAST algorithm/process, TBD precisely which one we will do.
match_naive: string matching using naive algorithm from text/classe (code from class).
match_RK: string matching using RK from the text.'''
