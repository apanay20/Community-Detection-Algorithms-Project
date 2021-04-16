import networkx as nx
import matplotlib.pyplot as plt

def showCommunities(G, partition):
    # Compute graph layout and plot the graph
    pos = nx.spring_layout(G)  
    plt.axis('off')
    nx.draw(G, pos, node_size=150, cmap=plt.cm.seismic, node_color=list(partition.values()))
    nx.draw_networkx_edges(G, pos, alpha=0.5)
    plt.show()
