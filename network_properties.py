import networkx as nx
from web import Network

if __name__ == '__main__':
    network_file = '../networks_lucie/graph1.edgelist'
    graph = nx.read_edgelist(network_file)
    net = Network(graph)
    
    print 'connectance = ', net.connectance()
    
    tls = net.get_trophic_levels()
    
    # top, top_preds = net.top_predators()
    # basal, basal_sps = net.basal()
    # for u,v in net.edges():
    #     if u in basal_sps and v in top_preds and tls[v] == 3:
    #         net.remove_edge(u,v)
    #
    # print 'new connectance = ', net.connectance()
    # print tls
    
    # print net.find_trophic_positions()
    
    print sorted(nx.algorithms.cycles.simple_cycles(net))