from __future__ import division
import numpy as np

int = np.int32

def weighted_high_median(a, wts):
    N = len(a)
    wtotal = 0
    wdiscardedlow = 0

    for i in range(N):
        wtotal += wts[i]

    nn = N
    while True:
        assert (nn > 0 and len(a) == nn)

        trial = sorted(a)[int(nn/2)]
        # Count up the weight to the left of and at the trial point.
        # Weight to the right of it isn't needed
        wleft = wtrial = 0
        for i in range(nn):
            if a[i] < trial:
                wleft += wts[i]
            elif a[i] == trial:
                wtrial += wts[i]

        if 2*(wdiscardedlow + wleft) > wtotal:
            # Trial value is too high
            ncandidates = 0
            #for i = 1:nn
            for i in range(nn):
                if a[i] < trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates

        elif 2*(wdiscardedlow + wleft + wtrial) > wtotal:
            # Trial value is just right
            return trial

        else:
            # Trial value is too low
            ncandidates = 0
            #for i = 1:nn
            for i in range(nn):
                if a[i] > trial:
                    a[ncandidates] = a[i]
                    wts[ncandidates] = wts[i]
                    ncandidates += 1
            nn = ncandidates
            wdiscardedlow += wleft+wtrial

        a=a[:nn]
        wts=wts[:nn]

def qn(data):
    # sort data
    data = np.sort(data)

    n = len(data)
    h = int(n/2) + 1
    k = int(h*(h-1)/2)

    left = np.arange(n+1,1,-1)
    right = np.full(n,n, dtype= np.int32)

    work = np.zeros(n, dtype=np.float64)
    weight = np.zeros(n, np.int32)
    P = np.zeros(n, np.int32)
    Q = np.zeros(n, np.int32)

    jhelp = int((n*(n+1))/2)
    knew = k+jhelp
    nL = jhelp
    nR = n*n
    found = False
    Qn = 0*data[0]

    while (nR-nL) > n:
        j = 1
        for i in range(1,n,1):
            if left[i] <= right[i]:
                weight[j-1] = right[i] - left[i] + 1
                jhelp = left[i] + int(weight[j-1]/2)
                work[j-1] =data[i] - data[n+1-jhelp-1]
                j += 1

        trial = weighted_high_median(work[:j-1], weight[:j-1])

        j=0
        for i in range(n-1, -1, -1):
            while (j < n) and (data[i]-data[n-j-1] < trial):
                j += 1
            P[i] = j

        j = n+1
        for i in range(n):
            while data[i]-data[n-j+2-1] > trial: # 55
                j -= 1
            Q[i] = j

        sumP = sum(P)
        sumQ = sum(Q)-n # 60

        if knew <= sumP:
            right[:] = P[:]
            nR = sumP
        elif knew > sumQ:
            left[:] = Q[:]
            nL = sumQ
        else:
            Qn = trial
            found = True
            break

    if found == False:
        j=1
        for i in range(1,n,1):
            if left[i] <= right[i]:
                for jj in range(left[i], right[i]+1, 1):
                    work[j-1] = data[i]-data[n-jj]
                    j += 1

        Qn = sorted(work[:j])[knew-nL-1]

    if n<10:
        nscale = [0, .399, .994, .512, .844, .611, .857, .669, .872][n-1]
    elif n%2 == 1:
        nscale = n/(n+1.4)
    else:
        nscale = n/(n+3.8)

    Qn = Qn*2.2219*nscale

    return Qn

def get_cluster_size_distribution(clusterid_to_taxa):
    clusterlen_to_frequency = {}
    for clusterid, taxa in clusterid_to_taxa.items():
        try:
            clusterlen_to_frequency[len(taxa)] += 1
        except:
            clusterlen_to_frequency[len(taxa)] = 1

    return [i for j in [[clusterlen] * frequency for clusterlen, frequency in clusterlen_to_frequency.items()] for i in j]
