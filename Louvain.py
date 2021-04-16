# Louvain Community Detection (Modularity-Based, Overlapping)
#
# Fast unfolding of communities in large networks
# (Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte, Etienne Lefebvre)
#
# # Undirected weighted/unweighted graph
#
# Algorithm:
#   1) Each node start in its own community
#   2) Repeat until convergence
#       for each node:
#           for each neighbor:
#               If adding node to its neighbor community increase modularity, do it
#   3) When converged, create an induced network
#       Each community becomes a node
#       Edge weight is the sum of weights of edges between them
#
#
# Authors: Andreas Panayiotou, Theodoros Kyriakou

import time
import numpy as np
import random
import networkx as nx
from Utilities.customPlot import showCommunities
from Utilities.evaluate import evaluate

# If yes, set random traversal visit sequence in phase1
randomBool = True
random.seed(130)

# Compute modularity Gain when moving node to specific community
def modularityGain(A, m, communities, node, nodeDegree, cluster):
    m2 = 2*m
    ki = nodeDegree
    # Get the nodes that are in current cluster
    nodesInCluster = np.where(communities==cluster)
    # Sum of degree of every node in given cluster
    Stot = 0
    # Sum of the weights of edges with both ends in current cluster
    Sin = 0
    # Sum of the weights of edges between node and other nodes in the community
    kin = 0
    for n1 in nodesInCluster[0]:
        Stot += np.sum(A[n1])
        for n2 in nodesInCluster[0]:
            Sin += 2*A[n1][n2]
            if n1 == node and n2 != node:
                kin += A[n1][n2]

    #communities[node] = storeOriginalComm
    DQ = (((Sin+(2*kin))/m2) - ((Stot+ki)/m2)**2) - (Sin/m2 - (Stot/m2)**2 - (ki/m2)**2)
    return DQ
 
# Compute Graph modularity 
def modularityGraph(G, partition):
    uniqueComm = list(set(partition))
    if len(uniqueComm) == 1:
        return 0

    m = 0
    try:
    	for n1,n2,data in G.edges(data=True):
            m += 2*float(data['weight'])
    except:
    	for n1,n2,data in G.edges(data=True):
            m += 2

    nodes = G.nodes()

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
                if partition[ni] == partition[nj]:
                    d = 1
                else:
                    d = 0
                Q += ( Aij - ((ki*kj)/(2*m)) )*d
    Q = Q / (2*m)

    return Q

# Return all unique communmities of node's neighbours
def getNeighCommunities(node, neighbours, partition):
    temp = []
    for n in neighbours:
        part = partition[n]
        if not part in temp:
       		temp.append(part)
    return temp;

# Run phase 1 of Louvain algorithm
def phase1(A, partition):
    # Get number of nodes
    countNodes = len(A)
    # If needed, set random traversal visit sequence
    visit = np.arange(countNodes)
    visit = random.sample(list(visit), countNodes) if randomBool else visit

    # Get sum of all degrees
    m = np.sum(A)/2

    # List to store partition of every run
    partitions = []
    bestGain = -1
    gains = []
    while True:
        for node in visit:
      	    # Get community of current node
            nodeComm = partition[node]
            # Get node degree
            nodeDegree = np.sum(A[node])
            # Find community and store of current node
            nodeComm = partition[node]
            # Evaluate gain in modularity in its current community
            Qcurrent = modularityGain(A, m, partition, node, nodeDegree, nodeComm)
            # Get all neighbours of node
            neighbours = np.argwhere(A[node]>0).flatten()
            # Get all unique communmities of node's neighbours
            neighCommunities = getNeighCommunities(node, neighbours, partition)

            # Evaluate gain in modularity by removing the node from its current
            # community and adding it to its neighbours' communities
            Qtotal = []
            for newComm in neighCommunities:
            	# Assign new community to node
                partition[node] = newComm
                # Calculate to medolarity in new community
                Qadd = modularityGain(A, m, partition, node, nodeDegree, newComm)
                # Recover the original cokmmuntiy of node
                partition[node] = nodeComm
                # Calculate the gain in modularity
                Qgain = Qadd - Qcurrent
                Qtotal.append(Qgain)
            
            # Get max modularity
            gain = max(Qtotal)
            # Assign the node to best community only if max gain is positive
            # else, staty in its current community
            if gain > 0.0:
    	        # Get neighbour's community with max modularity
                bestComm = neighCommunities[np.argmax(Qtotal)]
    	        # Set best community to the current node
                partition[node] = bestComm
            # Append current gain to gains array
            gains.append(gain)

        if max(gains) <= bestGain:
            break
        bestGain = max(gains)
        partitions.append(partition.copy())

    # Return best partition and gain
    return partitions[len(partitions)-2]

# Get community of every node(Rename community ids to start from 0)
def createLabels(countNodes, partition):
    # Get the number of unique partitions
    countComm = len(np.unique(partition))
    # Rename communities' id to starting from 0
    temp = dict(enumerate(np.unique(partition)))
    uniqueComm = dict((value,key) for key,value in temp.items())
    # Assign the communities to the nodes
    labels = []
    for n in partition:
    	labels.append(uniqueComm[n])

    return labels

# Run phase 2 of Louvain algorithm
def phase2(A, labels):
    uniqueLabels = list(set(labels))
    # Create new adj. matrix filled with zeros
    newAjx = np.zeros((len(uniqueLabels), len(uniqueLabels)))
    # Iterrate all edges in Adj. matrix
    # If two nodes in same community add one to community super node self-loop
    # If two nodes in different communitites add 1 to the edge connecting two communities 
    for i in range (0,len(A)):
     	for j in range (0,len(A)):
            if A[i][j] > 0:
                if labels[i] == labels[j]:
                    comm = labels[i]
                    indexComm = uniqueLabels.index(comm)
                    newAjx[indexComm][indexComm] += A[i][j]
                else:
                    comm1 = labels[i]
                    comm2 = labels[j]
                    indexComm1 = uniqueLabels.index(comm1)
                    indexComm2 = uniqueLabels.index(comm2)
                    newAjx[indexComm1][indexComm2] += A[i][j] 

    return newAjx

# Run louvain algorithm
def louvain(A):
	# Run phase 1 of Louvain algorithm and pass a partition where every node
    # located in its own community
    partition = phase1(A, np.arange(len(A)))
   	# Run phase 2 of Louvain algorithm until 1 node left in graph
    levels = []
    while True:
        labels = createLabels(len(A), partition)    
        A = phase2(A, labels)
        levels.append(labels.copy())
        # If one node left then stop the phase 2
        if len(A) == 1:
            break
        # Fill adj. matrix with zeros in diagonal to run phase1 again
        for i in range(0,len(A)):
        	for j in range(0,len(A)):
        		if i == j:
        			A[i][j] = 0
        # Run phase 1 of Louvain algorithm and pass a partition where super node
        # located in its own community
        partition = phase1(A, np.arange(len(A)))

    return levels

# Convert partition representation for plotting
def convertPartitionFormat(communities):
    temp = {}
    count = 0
    for n in communities:
    	temp[count] = n
    	count += 1
    return temp;

# Extract partition of every level
def buildPartitions(G,levels):
    partition = []
    # Add first level whuch contains all graph nodes
    partition.append(levels[0])
    # Unpack partitions from levels based on the level above
    for i in range(0,len(levels)-1):
        tempPartition = []
        label = levels[i+1] 
        for p in partition[i]:           
            tempPartition.append(label[p])
        partition.append(tempPartition)
    # Compute graph modularity for partition ov every level
    modularity = []

    for p in partition:
        modularity.append(modularityGraph(G,p))

    # Find max modularity
    bestMod = max(modularity)
    # Find index of max modularity
    bestModIndex = modularity.index(bestMod) 
    # Return the partition with max modularity and the value of modularity
    return partition[bestModIndex], max(modularity)

def printResults(communities,nodes):
    count = 1
    commUnique = set(communities)
    
    for c in commUnique:
        tempComm = []
        nodeCount = 0
        for n in communities:
            if n == c:
                tempComm.append(nodes[nodeCount])
            nodeCount += 1         
        print("Community:", count, "Nodes:", tempComm)
        count += 1

    print("Number of communities:",len(commUnique))

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

if __name__ == "__main__":
    path = "Datasets/karate.txt" 
    G = readGraph(path)
    A = nx.to_numpy_array(G.copy())

    print('Running Louvain on',path,"\n")
    start_time = time.time()
    levels = louvain(A)
    bestPartition, modularity = buildPartitions(G, levels)
    ExeTime = time.time() - start_time

    print("Execution time:","{:.2f}".format(ExeTime), "seconds")
    printResults(bestPartition, list(G.nodes()))
    print("Modularity:", modularity)
    evaluate(bestPartition, "Datasets/Truth/karateTrue.txt")

    finalPartition = convertPartitionFormat(bestPartition)
    showCommunities(G, finalPartition)