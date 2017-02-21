import collections
import itertools
import numpy as np

def sigmoid(score):
    return 1.0/(1.0+np.exp(score))

def ldn(a, b):
    """
    Leventsthein distance normalized
    :param a: word
    :type a: str
    :param b: word
    :type b: str
    :return: distance score
    :rtype: float
    """
    m = [];
    la = len(a) + 1;
    lb = len(b) + 1
    for i in range(0, la):
        m.append([])
        for j in range(0, lb): m[i].append(0)
        m[i][0] = i
    for i in range(0, lb): m[0][i] = i
    for i in range(1, la):
        for j in range(1, lb):
            s = m[i - 1][j - 1]
            if (a[i - 1] != b[j - 1]): s = s + 1
            m[i][j] = min(m[i][j - 1] + 1, m[i - 1][j] + 1, s)
    la = la - 1;
    lb = lb - 1
    return float(m[la][lb])# / float(max(la, lb))
    
def LD(x,y,lodict=None):
    """
    x is target, y is source
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    n,m = len(x),len(y)
    dp = np.zeros((n+1,m+1))
    pointers = np.zeros((n+1,m+1),np.int32)
    for i in range(1,n+1):
        if not lodict:
            dp[i,0] = dp[i-1,0]+1
        else:
            dp[i,0] = dp[i-1,0]+lodict['-',x[i-1]]
        pointers[i,0]=1
    for j in range(1,m+1):
        if not lodict:
            dp[0,j] = dp[0,j-1]+1
        else:
            dp[0,j] = dp[0,j-1]+lodict[y[j-1],'-']
        pointers[0,j]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            if not lodict:
                match = dp[i-1,j-1]
                if x[i-1] != y[j-1]:
                    match = match+1                
                insert = dp[i-1,j]+1
                delet = dp[i,j-1]+1
            else:
                match = dp[i-1,j-1]+lodict[x[i-1],y[j-1]]
                insert = dp[i-1,j]+lodict['-',x[i-1]]
                delet = dp[i,j-1]+lodict[y[j-1],'-']
            min_score = min([match,insert,delet])
            dp[i,j] = min_score
            pointers[i,j] = [match,insert,delet].index(min_score)
    alg = []
    i,j = n,m
    while(i>0 or j>0):
        pt = pointers[i,j]
        if pt==0:
            i-=1
            j-=1
            alg = [[x[i],y[j]]]+alg
        if pt==1:
            i-=1
            alg = [[x[i],'-']]+alg
        if pt==2:
            j-=1
            alg = [['-',y[j]]]+alg
    return dp[-1,-1], alg


def nw(x,y,lodict=None,gp1=-2.5,gp2=-1.75):
    """
    Needleman-Wunsch algorithm for pairwise string alignment
    with affine gap penalties.
    'lodict' must be a dictionary with all symbol pairs as keys
    and match scores as values.
    gp1 and gp2 are gap penalties for opening/extending a gap.
    Returns the alignment score and one optimal alignment.
    """
    n,m = len(x),len(y)
    dp = np.zeros((n+1,m+1))
    pointers = np.zeros((n+1,m+1),np.int32)
    for i in range(1,n+1):
        dp[i,0] = dp[i-1,0]+(gp2 if i>1 else gp1)
        pointers[i,0]=1
    for j in range(1,m+1):
        dp[0,j] = dp[0,j-1]+(gp2 if j>1 else gp1)
        pointers[0,j]=2
    for i in range(1,n+1):
        for j in range(1,m+1):
            if not lodict:
                if x[i-1] == y[j-1]:
                    match = dp[i-1,j-1]+1
                else:
                    match = dp[i-1,j-1]-1
            else:
                match = dp[i-1,j-1]+lodict[x[i-1],y[j-1]]
            insert = dp[i-1,j]+(gp2 if pointers[i-1,j]==1 else gp1)
            delet = dp[i,j-1]+(gp2 if pointers[i,j-1]==2 else gp1)
            max_score = max([match,insert,delet])
            dp[i,j] = max_score
            pointers[i,j] = [match,insert,delet].index(max_score)
    alg = []
    i,j = n,m
    while(i>0 or j>0):
        pt = pointers[i,j]
        if pt==0:
            i-=1
            j-=1
            alg = [[x[i],y[j]]]+alg
        if pt==1:
            i-=1
            alg = [[x[i],'-']]+alg
        if pt==2:
            j-=1
            alg = [['-',y[j]]]+alg
    return dp[-1,-1], alg
