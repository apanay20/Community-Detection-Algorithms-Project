# Girvan & Newman Community Detection (node-removal approache)
#
# Community structure in social and biological networks
# (M. Girvan and M. E. J. Newman)
#
# Undirected weighted/unweighted graph
#
# Algorithm:
# 	1) The betweenness of all existing edges in the network is calculated first.
# 	2) The edge with the highest betweenness is removed.
# 	3) The betweenness of all edges affected by the removal is recalculated.
# 	4) Steps 2 and 3 are repeated until no edges remain.
#
#
# Authors: Andreas Panayiotou, Theodoros Kyriakou

import time
import networkx as nx
import math
import random
import operator
from Utilities.customPlot import showCommunities
from Utilities.evaluate import evaluate

# If k is not None use k node samples to estimate betweenness.
# The value of k <= n where n is the number of nodes in the graph.
# Higher values give better approximation, but more running time.
centralityValue = None;

# Remove edges with max betweenness centrality from Graph until one of the connected components splits into two
def runStep(G):
    # Compute betweenness centrality of graph
    # Betweenness centrality of an edge e is the sum of the fraction of all-pairs shortest paths that pass through e.
    # Centrality of edge e is: SUM( c(s,t|e) / c(s,t) ), where c(s,t) is the number of shortest path, and
    # c(s,t|e) is the number of those paths passing  throught edge e. 
    centrality = nx.edge_betweenness_centrality(G, k=centralityValue, weight='weight')
    # Get number of connected components
    initialComponentsCount = nx.number_connected_components(G)
    newComponentsCount = initialComponentsCount

    # Run the proccess until one or more components are detached
    while newComponentsCount == initialComponentsCount:
        centrality = nx.edge_betweenness_centrality(G, k=centralityValue, weight='weight')
        # Find the edge which has the max centrality and then remove it
        sortedCentrality = sorted(centrality.items(), key=operator.itemgetter(1), reverse=True)
        # Find all edges which may have same max centrality
        temp = []
        for e,v in sortedCentrality:
            if v == sortedCentrality[0][1]:
                temp.append(e)
            else:
                break
        # Pick a random edge from edges with max modularity
        maxC = random.choice(temp)
        # Remove the edge
        G.remove_edge(maxC[0],maxC[1])

        # Get connected components after the removal
        newComponents = list(nx.connected_components(G))
        newComponentsCount = len(newComponents)
        # Compute betweenness centrality of nodes ONLY affected by the edge removal
    return newComponents

# Compute modularity of graph in current partition
def computeModularity(G, m, partition):
    if len(partition) == 1:
        return 0
    
    nodes = G.nodes()
    communities = [0] * len(nodes)
    for n in nodes:
        for p in partition:
            if n in p:
                communities[n] = partition.index(p)

    Q = 0
    for ni in nodes:
        for nj in nodes:
            if ni != nj:
                try:
                    ki = G.degree(ni,weight="weight")
                    kj = G.degree(nj,weight="weight")
                    # Aij is 1 if both nodes of current edge in same component 
                    if G.has_edge(ni, nj):
                        Aij = G.get_edge_data(ni,nj)['weight']
                    else:
                        Aij = 0
                except:                 
                    ki = 1
                    kj = 1
                    # Aij is 1 if both nodes of current edge in same component 
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

# Runs GirvanNewman algorithm and by maximizing modularity find the best graph partition
def GirvanNewman(G):
    edges = G.edges(data=True)
    m = 0
    try:
        for n1,n2,data in G.edges(data=True):
            m += float(data['weight'])
    except:
        for n1,n2,data in G.edges(data=True):
            m += 1

    Qbest, Qcurrent = 0.0, 0.0
    # Run until all edges have been removed from graph
    while G.number_of_edges() > 0:
        # Run one step of algorithm
        components = runStep(G)
        # Compute modularity of the current version of graph
        Qcurrent = computeModularity(G, m, components)
        # If new modularity is greater than the previous best then update list of communities
        if Qcurrent > Qbest:
            Qbest = Qcurrent
            communities = components
    return communities, Qbest

def readGraph(path):
    print("Reading Graph")
    graph = nx.Graph()
    tempEdges = []
    # Start getting edges
    with open(path) as file:
        for line in file:
            temp = line.strip("\n").split()
            if len(temp) == 2:
                if temp[0] != temp[1]:
                    tempEdges.append( (int(temp[0]), int(temp[1]), {'weight':1.0}) )
            else:
                if temp[0] != temp[1]:
                    tempEdges.append( (int(temp[0]), int(temp[1]), {'weight':float(temp[2])}) )
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

def printResults(communities, modularity, G):
    count = 1
    for c in communities:         
        print("Community:", count, "Nodes:", c)
        for n in c:
            G.nodes[n]['c'] = count
        count += 1
    print("Number of communities:",count-1)

# Convert partition format to plot the graph
def convertPartitionFormat(nodes, communities):
    temp = {}
    count = 0
    for c in communities:         
        for n in c:
        	temp[n] = count
        count += 1

    export = {}
    i = 0
    for n in nodes:
    	export[i] = temp[n]
    	i += 1
    return export

# Change the format of partition to compute evaluation 
def evaluateFormat(communities, nodesCount):
    count = 0
    ret = [0]*nodesCount

    for p in communities:
        for n in p:
            ret[n] = count
        count += 1

    return ret

if __name__ == "__main__":
    path = "Datasets/karate.txt"
    G = readGraph(path)
    tempG = G.copy()
    print('Running Girvan Newman on',path,"\n")

    start_time = time.time()
    communities, modularity = GirvanNewman(G)
    ExeTime = time.time() - start_time

    print("Execution time:","{:.2f}".format(ExeTime), "seconds")
    printResults(communities, modularity, tempG)
    print("Modularity:", modularity)
    evaluateFormat = evaluateFormat(communities,len(G.nodes()))
    evaluate(evaluateFormat, "Datasets/Truth/karateTrue.txt")

    finalPartition = convertPartitionFormat(tempG.nodes(), communities)
    showCommunities(tempG, finalPartition)