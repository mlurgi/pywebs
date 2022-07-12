
#import csv
import matplotlib.pyplot as plt

from config import SAVE_NETWORK, INTERVAL, START_YEAR, END_YEAR, JACCARD 
from config import HIGHLIGHTED, NODE_DIAMETER, REMOVED_NODES
from config import REMOVAL, REMOVED_SPECIES, REMOVAL_BETA, REMOVAL_EXT_THRESHOLD
from config import RGRAPH_SEED, RGRAPH_ITERATION_F, RGRAPH_COOLING_F
from network_processing import NetworkCreator
from plotting import DataStatsPlotter, NetworkPlotter
from utils import *

from web import Network

import networkx as nx

if __name__ == '__main__':
    dr = DatabaseReader()
    
    #dr.differences_between_datasets()
    
    observations, years_sps = dr.get_species_observations()
    #years = observations.keys()
    
    nc = NetworkCreator(dr)
    
    net = nc.create_whole_network(years_sps.keys())
    
    print net.order(), net.size(), net.connectance()
    
    nx.write_edgelist(net, '../output/net-adj-list')
    
#    clustering = True
#    nc.create_networks_from_data_and_habitats(observations, clustering)
    
#    mismatches = nc.create_networks_from_data_and_habitats(observations, clustering)
#    
#    print 'mismatches...'
#    
#    habitats = dr.read_habitats_names()
#    for species in sorted(mismatches.keys()):
#        habitats_names = set()
#        for h in mismatches[species]['habitats']:
#            habitats_names.add(habitats[h])
#            
#        print species, '; ', sorted(mismatches[species]['grids']),'; ', sorted(habitats_names)
#    
#    #nc.create_networks_from_data(observations)
    
#    nets_dict = nc.get_networks()
#
#    #the following code writes the output tables according to above conditions
#    plots_data = create_output_csv(nets_dict)
    
    #calculate_migrant_species(nets_dict)
    
#    masses = dr.read_species_mass()
#    habitats = dr.read_habitats()
#    sps_in_post = degree_per_periods(nets_dict, masses, habitats)
    
    if JACCARD:
        calculate_jaccard_matrix(nets_dict)
    
#    output = csv.DictReader(open('../output/output.csv', 'rb'), delimiter=',')
#    years = []
#    plotting_data = dict()
#    for row in output:
#        year = row['year']
#        years.append(year)
#        
#        habitat = row['habitat']
    
#    plots_data_new = dict()
#    for year in sorted(plots_data.keys()):
#        
#        if not plots_data_new.has_key('years'):
#            plots_data_new['years'] = [year]
#        else:
#            plots_data_new['years'].append(year)
#        
#        for habitat in plots_data[year].keys():
#            
#            if not plots_data_new.has_key(habitat):
#                plots_data_new[habitat] = dict()
#            
#            for elevation in plots_data[year][habitat].keys():
#                
#                species_number = plots_data[year][habitat][elevation]['species_number']
#                
#                min_biomass = plots_data[year][habitat][elevation]['min_biomass']
#                if min_biomass == None:
#                    min_biomass = 0.0
#                     
#                max_biomass = plots_data[year][habitat][elevation]['max_biomass']
#                if max_biomass == None:
#                    max_biomass = 0.0
#                
#                mean_biomass = plots_data[year][habitat][elevation]['mean_biomass']
#                if mean_biomass == None:
#                    mean_biomass = 0.0
#                    
#                linkage_density = plots_data[year][habitat][elevation]['linkage_density']
#                if linkage_density == None:
#                    linkage_density = 0.0
#                    
#                connectance = plots_data[year][habitat][elevation]['connectance']
#                if connectance == None:
#                    connectance = 0.0
#                    
#                mean_fcl = plots_data[year][habitat][elevation]['mean_fcl']
#                if mean_fcl == None:
#                    mean_fcl = 0.0
#                
#                if not plots_data_new[habitat].has_key(elevation):
#                    plots_data_new[habitat][elevation] = dict()
#                    plots_data_new[habitat][elevation]['species_number'] = [species_number]
#                    plots_data_new[habitat][elevation]['min_biomass'] = [min_biomass]
#                    plots_data_new[habitat][elevation]['max_biomass'] = [max_biomass]
#                    plots_data_new[habitat][elevation]['mean_biomass'] = [mean_biomass]
#                    plots_data_new[habitat][elevation]['linkage_density'] = [linkage_density]
#                    plots_data_new[habitat][elevation]['connectance'] = [connectance]
#                    plots_data_new[habitat][elevation]['mean_fcl'] = [mean_fcl]
#    
#                else:
#                    plots_data_new[habitat][elevation]['species_number'].append(species_number)
#                    plots_data_new[habitat][elevation]['min_biomass'].append(min_biomass)
#                    plots_data_new[habitat][elevation]['max_biomass'].append(max_biomass)
#                    plots_data_new[habitat][elevation]['mean_biomass'].append(mean_biomass)
#                    plots_data_new[habitat][elevation]['linkage_density'].append(linkage_density)
#                    plots_data_new[habitat][elevation]['connectance'].append(connectance)
#                    plots_data_new[habitat][elevation]['mean_fcl'].append(mean_fcl)
#    
#    
#    print plots_data_new
#    
#    figpl = DataStatsPlotter()
#    
#    years = plots_data_new['years']
#    y1 = plots_data_new[5][0]['species_number']
#    y2 = plots_data_new[5][1]['species_number']
#    y3 = plots_data_new[5][2]['species_number']
#    x_label = 'Years'
#    y_label = 'Species number'
#    title = 'Year vs Number of species for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
#    
#    y1 = plots_data_new[5][0]['min_biomass']
#    y2 = plots_data_new[5][1]['min_biomass']
#    y3 = plots_data_new[5][2]['min_biomass']
#    x_label = 'Years'
#    y_label = 'Min biomass'
#    title = 'Year vs Min Biomass for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
#    
#    y1 = plots_data_new[5][0]['max_biomass']
#    y2 = plots_data_new[5][1]['max_biomass']
#    y3 = plots_data_new[5][2]['max_biomass']
#    x_label = 'Years'
#    y_label = 'Max biomass'
#    title = 'Year vs Max Biomass for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
#    
#    y1 = plots_data_new[5][0]['mean_biomass']
#    y2 = plots_data_new[5][1]['mean_biomass']
#    y3 = plots_data_new[5][2]['mean_biomass']
#    x_label = 'Years'
#    y_label = 'Mean biomass'
#    title = 'Year vs Mean Biomass for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
#    
#    y1 = plots_data_new[5][0]['linkage_density']
#    y2 = plots_data_new[5][1]['linkage_density']
#    y3 = plots_data_new[5][2]['linkage_density']
#    x_label = 'Years'
#    y_label = 'Linkage density'
#    title = 'Year vs Linkage Density for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
#    
#    y1 = plots_data_new[5][0]['connectance']
#    y2 = plots_data_new[5][1]['connectance']
#    y3 = plots_data_new[5][2]['connectance']
#    x_label = 'Years'
#    y_label = 'Connectance'
#    title = 'Year vs Connectance for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
#    
#    y1 = plots_data_new[5][0]['mean_fcl']
#    y2 = plots_data_new[5][1]['mean_fcl']
#    y3 = plots_data_new[5][2]['mean_fcl']
#    x_label = 'Years'
#    y_label = 'Mean Food Chain Length'
#    title = 'Year vs MFCL for 3 altitude ranges in habitat 5'
#    figpl.plot_three_series_points(years, years, years, years, y1, y2, y3, x_label, y_label, title, style_lines=True)
    
    
#    #here we store in a file the network indicated in SAVE_NETWORK
#    if clustering:
#        try:
#            y, h, e, s = SAVE_NETWORK.split('-')
#            year = int(y)
#            cluster = int(h)
#            elevation = int(e)
#            net_original = nets_dict[year][cluster][elevation][s]
#            
#            net = net_original.get_copy_removed_nodes(REMOVED_NODES)
#            
#            filepath = '../output/'+SAVE_NETWORK+'.edgelist'
#            nx.write_pajek(net, filepath)
#        except:
#            print 'There are no networks for the specified parameters in SAVE_NETWORK'
#    else:
#        try:
#            y, e, s = SAVE_NETWORK.split('-')
#            year = int(y)
#            elevation = int(e)
#            net = nets_dict[year][elevation][s]
#            filepath = '../output/'+SAVE_NETWORK+'.graphml'
#            nx.write_graphml(net, filepath)
#        except:
#            print 'There are no networks for the specified parameters in SAVE_NETWORK'
    
    
    
    #the networks structure is organised in the following way:
    #first index = year
    #second = cluster number (1-6)
    #third = altitude interval code (depends upon INTERVAL_SLOTS values from 0 to INTERVAL_SLOTS-1)
    #fourth = season of the network 'ss' 'aw' 'ay'
    
#    netpl = NetworkPlotter()
#    net = nets_dict[2000][1][2]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    
#    for n in net.nodes():
#        if n in sps_in_post[0]:
#            net.node[n]['in_post'] = True
#        else:
#            net.node[n]['in_post'] = False
#            
#    #net.write_network_files('../output/2001-0')
#    
#    netpl.plot_network(net, node_diameter='biomass')
    

    
#    netpliii = NetworkPlotter()
#    net = nets_dict[2001][1][2]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    
#    for n in net.nodes():
#        if n in sps_in_post[2]:
#            net.node[n]['in_post'] = True
#        else:
#            net.node[n]['in_post'] = False
#            
#    net.write_network_files('../output/2001-2')
#    
#    netpliii.plot_network(net)
#    
#    netpl_ii = NetworkPlotter()
#    net = nets_dict[1990][1][1]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    netpl_ii.plot_network(net)
#    
#    netpl_iii = NetworkPlotter()
#    net = nets_dict[1990][1][2]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    netpl_iii.plot_network(net)
#    
#    netpl = NetworkPlotter()
#    net = nets_dict[2000][1][0]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    netpl.plot_network(net)
#    
#    netpl_ii = NetworkPlotter()
#    net = nets_dict[2000][1][1]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    netpl_ii.plot_network(net)
#    
#    netpl_iii = NetworkPlotter()
#    net = nets_dict[2000][1][2]['all'].get_copy_removed_nodes(REMOVED_NODES)
#    netpl_iii.plot_network(net)
#    
    
    
    
#    plt.show()
    
#    
    netpl_iv = NetworkPlotter()
    #net = nets_dict[1980][1][0]['ss']
    
    network_file = '../output/pyrenees.graphml'
    nx.write_graphml(net, network_file)
    
    mod, n_mods = net.modularity_rgraph(RGRAPH_SEED, RGRAPH_ITERATION_F, RGRAPH_COOLING_F, 0)
    
    print 'modularity = ', mod
    print 'modules = ', n_mods
    
#    graph = nx.read_graphml(network_file)
#    net = Network(graph)
#    
#    #print nets_dict
#    
    for n in net.nodes():
        if n == 'Gyps fulvus':
            print 'gyps is here!!!'
        net.node[n]['in_post'] = False
        net.node[n]['mut_prod'] = False
        net.node[n]['mut'] = False
    
    netpl_iv.plot_network(net, highlight='Gyps fulvus', node_diameter=None)
    
    
    
#    
#    
#    netpl_v = NetworkPlotter()
#    net = nets_dict[1995][1][2]['aw']
#    netpl_v.plot_network(net)
#    
#    netpl_vi = NetworkPlotter()
#    net = nets_dict[1995][3][2]['aw']
#    netpl_vi.plot_network(net)
#    
#    netpl_vii = NetworkPlotter()
#    net = nets_dict[1995][4][2]['aw']
#    netpl_vii.plot_network(net)
#    
#    netpl_viii = NetworkPlotter()
#    net = nets_dict[1995][5][2]['aw']
#    netpl_viii.plot_network(net)
#    
#    
#    
#    netpl_ix = NetworkPlotter()
#    net = nets_dict[1995][3][0]['aw']
#    netpl_ix.plot_network(net)
#    
#    netpl_x = NetworkPlotter()
#    net = nets_dict[1995][3][1]['aw']
#    netpl_x.plot_network(net)
#    
#    netpl_xi = NetworkPlotter()
#    net = nets_dict[1995][3][2]['aw']
#    netpl_xi.plot_network(net)
    
    
#    if REMOVAL:
#        sps_to_remove = dr.get_species_to_remove()
#        if len(sps_to_remove) == 0:
#            sps_to_remove = set(REMOVED_SPECIES)
#        
#        print sps_to_remove
#        if len(sps_to_remove & set(net.nodes())) != 0:
#            net_pl_removal_i = NetworkPlotter()
#            net_removed, r = net.robustness_to_removal(sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
#            net_pl_removal_i.plot_network(net_removed)
#            
    #print net.longest_path_length(True)
 
#    netpl_iv = NetworkPlotter()
#    net = nets_dict[1999][5][2]['aw']
#    netpl_iv.plot_network(net, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    if REMOVAL:
#        if len(sps_to_remove & set(net.nodes())) != 0:
#            net_pl_removal_ii = NetworkPlotter()
#            net_removed, r = net.robustness_to_removal(sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
#            net_pl_removal_ii.plot_network(net_removed, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    
#    netpl_v = NetworkPlotter()
#    net = nets_dict[1990][3][1]['ss']
#    netpl_v.plot_network(net, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    if REMOVAL:
#        if len(sps_to_remove & set(net.nodes())) != 0:
#            net_pl_removal_ii = NetworkPlotter()
#            net_removed, r = net.robustness_to_removal(sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
#            net_pl_removal_ii.plot_network(net_removed, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#
#    netpl_vi = NetworkPlotter()
#    net = nets_dict[1990][3][1]['aw']
#    netpl_vi.plot_network(net, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    if REMOVAL:
#        if len(sps_to_remove & set(net.nodes())) != 0:
#            net_pl_removal_ii = NetworkPlotter()
#            net_removed, r = net.robustness_to_removal(sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
#            net_pl_removal_ii.plot_network(net_removed, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    
#    netpl_vii = NetworkPlotter()
#    net = nets_dict[1990][3][2]['ss']
#    netpl_vii.plot_network(net, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    if REMOVAL:
#        if len(sps_to_remove & set(net.nodes())) != 0:
#            net_pl_removal_ii = NetworkPlotter()
#            net_removed, r = net.robustness_to_removal(sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
#            net_pl_removal_ii.plot_network(net_removed, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#
#    netpl_viii = NetworkPlotter()
#    net = nets_dict[1990][3][2]['aw']
#    netpl_viii.plot_network(net, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)
#    if REMOVAL:
#        if len(sps_to_remove & set(net.nodes())) != 0:
#            net_pl_removal_ii = NetworkPlotter()
#            net_removed, r = net.robustness_to_removal(sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
#            net_pl_removal_ii.plot_network(net_removed, highlight=HIGHLIGHTED, node_diameter=NODE_DIAMETER)



    #print net.longest_path_length(True)
#    
#    netpl_v = NetworkPlotter()
#    net = nets_dict[1984][5][1]['ss']
#    netpl_v.plot_network(net)
#    #print net.longest_path_length(True)
#    
#    netpl_vi = NetworkPlotter()
#    net = nets_dict[1984][5][2]['ss']
#    netpl_vi.plot_network(net)
#    #print net.longest_path_length(True)
    
  
#    links, distro = net.get_cummulative_degree_dist()
#    figpl.plot_single_series(links, distro, 'Degree', 'Number of nodes', 'Degree distribution for the 2000 summer-spring network', 'ro', True, fitting='all')
#    
#    net = nc.get_networks()[2001][1][1]['ss']
#    print net.longest_path_length(True)
#    
#    figpl.plot_network(net)
#    links, distro = net.get_cummulative_degree_dist()
#    figpl.plot_single_series(links, distro, 'Degree', 'Number of nodes', 'Degree distribution for the 2000 summer-spring network', 'ro', True, fitting='powerlaw')
    
#    #plot years vs species numbers
#    years = nc.get_years_data()
#    y1, y2 = nc.get_nsp_data()
#    x_label = 'Years'
#    y_label = 'S'
#    title = 'Year vs Number of Species'
#    figpl.plot_series_points(years, years, years, y1, y2, x_label, y_label, title, style_lines=True)
#    
#    #plot species vs links numbers
#    x1, x2 = nc.get_nsp_data()
#    y1, y2 = nc.get_nlinks_data()
#    x_label = 'S'
#    y_label = 'L'
#    title = 'Number of species vs Number of links, for '+str(INTERVAL)+'-year aggregated networks between '+str(START_YEAR)+' and '+str(END_YEAR)
#    figpl.plot_series_points(years, x1, x2, y1, y2, x_label, y_label, title, point_labels=True)
#
    #plot connectance vs number of species
#    x1, x2 = nc.get_nsp_data()
#    y1, y2 = nc.get_cons_data()
#    y_label = 'Connectance'
#    x_label = 'S'
#    title = 'Number of species vs C = L/S^2, for '+str(INTERVAL)+'-year aggregated networks between '+str(START_YEAR)+' and '+str(END_YEAR)
#    figpl.plot_series_points(years, x1, x2, y1, y2, x_label, y_label, title, point_labels=True)
#        
#    #plot connectance vs number of observations
#    years = nc.get_years_data()
#    x1, x2 = nc.get_obs_data()
#    y1, y2 = nc.get_cons_data()
#    x_label = 'Number of observations'
#    y_label = 'Connectance'
#    title = 'Number of observations vs C = L/S^2, for '+str(INTERVAL)+'-year aggregated networks between '+str(START_YEAR)+' and '+str(END_YEAR)
#    figpl.plot_series_points(years, x1, x2, y1, y2, x_label, y_label, title, point_labels=True)
#    
#    x,y,z = process_data_observations(years_sps, observations)
#    
#    years = observations.keys()
#    figpl.plot_observations_per_year(years, x, y, z)
    
    plt.show()
    

    print '\n\n\nUsing this application you can compare two given networks of the ones just obtained based on their species composition.\n'
    correct_input = False
    while not correct_input:
        nets_to_compare = raw_input("Enter the names of 2 networks to compare (format: year-habitat-elevation-season), separated by space: ")
        
        try:
            net1, net2 = nets_to_compare.split()
        except:
            print 'Input error: You must specify two (and only two) network names following the format specified\n'
            continue
        
        #we obtain the information for network 1
        correct_1 = False
        try:
            y1, h1, e1, s1 = net1.split('-')
            try:
                net1_y = nets_dict[int(y1)]
                try:
                    net1_y_h = net1_y[int(h1)]
                    try:
                        net1_y_h_e = net1_y_h[int(e1)]
                        try:
                            net1_final = net1_y_h_e[s1] 
                            correct_1 = True
                        except:
                            print 'Input error: There are no networks for season ', s1, ' in the elevation ', e1, ' for the habitat ', h1, ' in ',y1
                    except:
                        print 'Input error: There are no networks in the elevation ', e1, ' for the habitat ', h1, ' in ',y1
                except:
                    print 'Input error: There are no networks for the habitat ', h1, ' in ', y1
            except:
                print 'Input error: There are no networks for the year: ' + y1
        except:
            print 'Input error: Name for network 1 has incorrect format'
        
        #we obtain the info for network 2
        correct_2 = False
        try:
            y2, h2, e2, s2 = net2.split('-')
            try:
                net2_y = nets_dict[int(y2)]
                try:
                    net2_y_h = net2_y[int(h2)]
                    try:
                        net2_y_h_e = net2_y_h[int(e2)]
                        try:
                            net2_final = net2_y_h_e[s2] 
                            correct_2 = True
                        except:
                            print 'Input error: There are no networks for season ', s2, ' in the elevation ', e2, ' for the habitat ', h2, ' in ', y2
                    except:
                        print 'Input error: There are no networks in the elevation ', e2, ' for the habitat ', h2, ' in ',y2
                except:
                    print 'Input error: There are no networks for the habitat ', h2, ' in ', y2
            except:
                print 'Input error: There are no networks for the year: ' + y2
        except:
            print 'Input error: Name for network 2 has incorrect format'
        
        
        print net1,net2
        
        correct_input = correct_1 and correct_2
        if not correct_input:
            continue
        
        print 'Network comparison: \n'
        
        sps_1 = set(net1_final.nodes())
        sps_2 = set(net2_final.nodes())
        both_species = sorted(sps_1 & sps_2)
        
        if len(both_species) == 0:
            print 'The networks have no species in common.'
        else:
            print 'Species that are in both networks: \n'
            for i in both_species:
                print i
        
        print'\n\n'
        
        exclusive_1 = sorted(sps_1 - sps_2)
        if len(exclusive_1) == 0:
            print 'All of the species in network 1 are also present in network 2.'
        else: 
            print 'Species in network 1 that are not in network 2: \n'
            for i in exclusive_1:
                print i
        
        print'\n\n'
        
        exclusive_2 = sorted(sps_2 - sps_1)
        if len(exclusive_2) == 0:
            print 'All of the species in network 2 are also present in network 1.'
        else:
            print 'Species in network 2 that are not in network 1: \n'
            for i in exclusive_2:
                print i
        
        print'\n\n'
        
        print 'Writing the comparison values to the output file ' + net1 + '_vs_' + net2 + '.csv...'
        
        create_compare_output(net1, net1_final, net2, net2_final)
        
        quit = ''
        while quit != 'yes' and quit != 'no':
            quit = raw_input("Do you want to quit the application? (yes/no) ")
            if quit == 'yes':
                print 'Thank you for choosing Pyrenees networks for your analyses...'
            if quit == 'no':
                correct_input = False
        

    