import sklearn.metrics as metrics

def ΜΙ(truth, algorithm):
    return metrics.mutual_info_score(truth, algorithm)

def NMI(truth, algorithm):
    return metrics.normalized_mutual_info_score(truth, algorithm)

def ARI(truth, algorithm):
    return metrics.adjusted_rand_score(truth, algorithm)

def RI(truth, algorithm):
    return metrics.rand_score(truth, algorithm)

# Read grouth-true communities and sort the nodes in list
def readPredictedCommunities(path):
    communities = {}
    cc = []
    with open(path) as file:
        for line in file:
            temp = line.strip("\n").split()
            c = temp.pop(0).strip()
            cc.append(int(c))
            for n in temp:
                communities[int(n)] = int(c)

    sort = sorted(communities.items())
    index = 0
    ret = [0]*len(sort)
    for n,c in sort:
        ret[index] = c
        index += 1

    return ret

def evaluate(algorithm,path):
    truth = readPredictedCommunities(path)
    print("============\n","EVALUATION","\n============")
    print("ΜΙ Score:",ΜΙ(truth,algorithm))
    print("NMI Score:",NMI(truth,algorithm))
    print("RI Score:",RI(truth,algorithm))
    print("ARI Score:",ARI(truth,algorithm))
