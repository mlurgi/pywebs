

import csv
import numpy as np
import matplotlib.pyplot as plt

from utils import DatabaseReader
from web import Network
from network_processing import NetworkCreator
from networkx.linalg.graphmatrix import adjacency_matrix

from plotting import NetworkPlotter

if __name__ == '__main__':
    dr = DatabaseReader()
    nc = NetworkCreator(dr)
    
    years = [2000, 2001, 2002, 2003, 2004]
    
    grids = dr.read_species_per_grid()
    
    #print grids
    
    sp_locations = dict()
    
    for y in grids.keys():
        if not y in years:
            continue
        for sp in grids[y].keys():
            
            for c in grids[y][sp]:
                if not c in sp_locations.keys():
                    sp_locations[c] = set([sp]) 
                else:
                    sp_locations[c].add(sp)
      
##### this is to have it per species, but we need it per grid cell      
#             if not sp in sp_locations.keys():
#                 sp_locations[sp] = grids[y][sp]
#             else:
#                 sp_locations[sp] = sp_locations[sp].union(grids[y][sp])
#
#      quad_set = sorted(sp_locations.values()[0].union(*sp_locations.values()))
    
    
    for gr in sp_locations.keys():
        print gr
        nc.create_networks_from_cell_and_species(gr, sp_locations[gr])
        
    
    nets = nc.get_networks()  
    
#     netpl = NetworkPlotter()
#     net = nets['DH00']['all']
#     netpl.plot_network(net, node_diameter='degree') 
#     plt.show()
    
    repetitions = 50
     
    for i in range(1,(repetitions+1)):      
        grid_sequence = np.random.permutation(sp_locations.keys())
         
        #merged = Network()
        species_agg = set()
        islands = 1
        for gr in grid_sequence:
            current_net = nets[gr]['all']
             
#             merged.add_nodes_from(current_net.nodes(data=True))
#             merged.add_edges_from(current_net.edges(data=True))
#              
             
            species_agg = species_agg.union(current_net.nodes())
#             
#             print species_agg
            
            a,b,c,d = nc.create_nets_from_seasonality(species_agg)
            
            z = d.get_copy_removed_nodes(d.basal()[1])
            
            np.savetxt('/home/miguel/Documents/networks-pyrennees/network-'+str(i)+'-'+str(islands)+'.txt', adjacency_matrix(z), fmt='%d')
            islands += 1
        
    #print d.basal()    
        #print 'connectance = ', str(d.connectance())
        
#     
#     for key, value in grids.items():
#         writer = csv.writer(open('dict-'+ str(key) +'.csv', 'wb'))
#         writer.writerow(['species', 'locations'])
#         
#         for k2, v2 in value.items():
#             writer.writerow([k2, list(v2)])