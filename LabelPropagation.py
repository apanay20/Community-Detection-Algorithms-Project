# Label Propagation(LPA) Community Detection (Overlap Propagation Algorithm)
#
# Near linear time algorithm to detect community structures in large-scale networks
# (Usha Nandini Raghavan, Reka Albert and Soundar Kumara)
#
# Undirected unweighted graph
#
# Algorithm:
#   1) Every node is initialized with a unique community label (an identifier).
#          These labels propagate through the network.
#   2) At every iteration of propagation, each node updates its label to the one that the maximum numbers of its neighbours belongs to.
#          LPA reaches convergence when each node has the majority label of its neighbours.
#   3) LPA stops if either convergence, or the user-defined maximum number of iterations is achieved.
#
# Authors: Andreas Panayiotou, Theodoros Kyriakou

import time
import random
import networkx as nx
from collections import defaultdict
from Utilities.customPlot import showCommunities
from Utilities.evaluate import evaluate

ITERATIONS = 200

def readGraph(path):
    print("Reading Graph")
    graph = nx.Graph()
    tempEdges = []
    # Start getting edges
    with open(path) as file:
        for line in file:
            temp = line.strip("\n").split()
            if temp[0] != temp[1]:
                tempEdges.append( (int(temp[0]), int(temp[1]), {'weight':1.0}) )
    # Create a list with nodes, take the unique values and sort them
    tempNodes = []
    for (n1,n2,w) in tempEdges:
        tempNodes.append(n1)
        tempNodes.append(n2)
    nodes = sorted(list(set(tempNodes)))
    # Create a networkx graph from nodelist and edgelist
    graph.add_nodes_from(nodes)
    graph.add_edges_from(tempEdges)
    # Relabel nodes in case they are not sequential
    count = 0
    mapping = {}
    for n in nodes:
        mapping[n] = count
        count += 1

    return nx.relabel_nodes(graph, mapping)

def checkAssignment(graph):
	# For every node in graph
    for node in list(graph.nodes(data=False)):
        occurrence = {}
        # Iterate throught all neigbors
        for neighbor in graph.neighbors(node):
        	# Get label of neigbor
            neighborLabel = graph.nodes[neighbor]['label']
            # Increase occurence of current label if exists or initialize with zero
            occurrence[neighborLabel] = occurrence.setdefault(neighborLabel, 0) + 1
        # Sort neighbors labels based on label occurence 
        srt = sorted(occurrence.items(), key=lambda item: item[1], reverse=True)
        
        # Check if current node has as label the max label of its neigbors labels. 
        if graph.nodes[node]['label'] != srt[0][0]:
            return False
    return True

def LPA(graph):
    # Initialy assign every node to its own community 
    for node, data in list(graph.nodes(data=True)):
        data['label'] = node

    counter = 1
    # Iterative label propagation process until max iterations or
    # LPA reaches convergence(every ) 
    while counter <= ITERATIONS:
    	# For every node in graph
        for node in list(graph.nodes()):
            occurrence = {}
            # Iterate throught all neigbors 
            for neighbor in graph.neighbors(node):
            	# Get label of neigbor
                neighborLabel = graph.nodes[neighbor]['label']
                # Increase occurence of current label if exists or initialize with zero.
                # If exists increase by 1, else set to zero and increase by 1
                occurrence[neighborLabel] = occurrence.setdefault(neighborLabel, 0) + 1

            # Find the label with the largest count
            maxLabel = sorted(occurrence.items(), key=lambda item: item[1], reverse=True)[0][1]
            # Insert to labels array all labels with equal max count
            labels = [k for k,v in occurrence.items() if v == maxLabel]
            # Pick a random label and assign to current node
            graph.nodes[node]['label'] = random.choice(labels)
        # If every node has a label that the maximum number of their neighbors have, then stop.
        if checkAssignment(graph):
            return
        counter += 1

# Extract and return communities from networkx graph
def getCommunities(graph):
    communities = defaultdict(dict)
    for v in graph:   
        communities[graph.nodes[v]['label']] = []
    for v in graph:   
        communities[graph.nodes[v]['label']].append(v)
    return communities

def printResults(communities):
    count = 1
    for c in communities.keys(): 
        print("Community:", count, "Nodes:", communities[c])
        count += 1    

# Convert partition format to plot the graph
def convertPartitionFormat(nodes, communities):
    temp = {}
    count = 0
    for c in communities.keys():         
        for n in communities[c]:
            temp[n] = count
        count += 1

    export = {}
    i = 0
    for n in nodes:
        export[i] = temp[n]
        i += 1
    return export

# Compute modularity of final partition
def modularity(G, partition):
    if len(partition) == 1:
        return 0.0

    m = 0
    for e in G.edges():
        m += 1

    nodes = G.nodes()
    communities = [0] * len(nodes)

    for comm,nodes in partition.items():
        for n in nodes:
            communities[n] = comm

    Q = 0
    for ni in nodes:
        for nj in nodes:
            if ni != nj:
                ki = G.degree(ni)
                kj = G.degree(nj) 
                if G.has_edge(ni, nj):
                    Aij = 1
                else:
                    Aij = 0
                # d is 1 if both nodes of current edge in same community
                if communities[ni] == communities[nj]:
                    d = 1
                else:
                    d = 0
                Q += ( Aij - ((ki*kj)/(2*m)) )*d
    Q = Q / (2*m)
    return Q

# Change the format of partition to compute evaluation 
def evaluateFormat(communities, nodesCount):
    count = 0
    newComm = {}
    for c in communities.keys():
        newComm[count] = communities[c]
        count += 1        

    ret = [0]*nodesCount
    for c in newComm.keys():
        for n in newComm[c]:
            ret[n] = int(c)
    return ret

if __name__ == "__main__":
    path = "Datasets/karate.txt"
    G = readGraph(path)
    print('Running Label Propagation on',path,"\n")

    start_time = time.time()
    LPA(G)
    ExeTime = time.time() - start_time
    
    communities = getCommunities(G)

    print("Execution time:","{:.2f}".format(ExeTime), "seconds")
    printResults(communities)
    print("Number of communities:",len(communities.keys()))
    print("Modularity:", modularity(G, communities))
    evaluateFormat = evaluateFormat(communities,len(G.nodes()))
    evaluate(evaluateFormat, "Datasets/Truth/karateTrue.txt")

    finalPartition = convertPartitionFormat(G.nodes(), communities)
    showCommunities(G, finalPartition)