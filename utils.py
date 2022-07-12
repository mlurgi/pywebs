
from files_data_reader import *
from database_reader import *
from numpy import log, exp, arange

import threading
import scipy
import networkx as nx

import matplotlib.pyplot as plt

from web import NotInvadableNetwork

from config import REMOVED_NODES, ROBUSTNESS_INTERVAL, CALCULATE_SEASONS
from config import ROBUSTNESS_MAX_REMOVED, ROBUSTNESS_CUMULATIVE, ROBUSTNESS_BETA, ROBUST_EXT_THRESHOLD, CALCULATE_ROBUSTNESS 
from config import INVASION, INTRODUCED_GROUP, BIOMASS_INTRODUCED, INTRODUCED_SPECIES
from config import FRACTION_OF_GENERALISM, FRACTION_OF_PREDATORS, INVASION_EXT_THRESHOLD
from config import REMOVAL, REMOVAL_BETA, REMOVAL_EXT_THRESHOLD, REMOVED_SPECIES
from config import RGRAPH_MOD, RGRAPH_SEED, RGRAPH_ITERATION_F, RGRAPH_COOLING_F, RGRAPH_RANDOMIZATIONS
from config import INVASION_SCENARIOS, INTRODUCED_GROUPS, PARAMETERS_OFFSET

    
def plfit_lsq(x,y):
    """
    Returns A and B in y=Ax^B
    http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
    """
    n = len(x)
    btop = n * (log(x)*log(y)).sum() - (log(x)).sum()*(log(y)).sum()
    bbottom = n*(log(x)**2).sum() - (log(x).sum())**2
    b = btop / bbottom
    a = ( log(y).sum() - b * log(x).sum() ) / n

    A = exp(a)
    return A,b


def create_output_csv(nets_dict):
    
    print 'number of networks = ', str(len(nets_dict.keys()) )
    
    if CALCULATE_ROBUSTNESS:
        header_names = ['year','habitat','elevation','season','S', 'L', 'L/S','C', 'T', 'B', 'I', 'Ca', 'Loop', 'NCycles', 'O', 'T-B', 'T-I', 'I-I', 'I-B', 'GenSD', 'VulSD', 'MxSim', 'MaxChainLength', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'MeanShortestPathLength', 'complexity', 'components', 'cc', 'cc_rnd', 'diameter', 'diameter_rnd', 'radius', 'radius_rnd', 'compartmentalisation', 'modularity_rgraph', 'modules_rgraph', 'modularity_rgraph_randoms', 'sd_modularity_rgraph_randoms', 'chi_square', 'p_value', 'mean_mass_ratio', 'min_biomass', 'max_biomass', 'mean_biomass', 'mean_ratio', 'variance_ratio', 'skew_ratio', 'kurt_ratio', 'R_removed', 'removed_species', 'R_CA_10', 'R_CA_20', 'R_CA_30', 'R_CA_40', 'R_CA_50', 'R_CA_60', 'R_CD_10', 'R_CD_20', 'R_CD_30', 'R_CD_40', 'R_CD_50', 'R_CD_60', 'R_MA_10', 'R_MA_20', 'R_MA_30', 'R_MA_40', 'R_MA_50', 'R_MA_60', 'R_MD_10', 'R_MD_20', 'R_MD_30', 'R_MD_40', 'R_MD_50', 'R_MD_60','R_TPA_10', 'R_TPA_20', 'R_TPA_30', 'R_TPA_40', 'R_TPA_50', 'R_TPA_60', 'R_TPD_10', 'R_TPD_20', 'R_TPD_30', 'R_TPD_40', 'R_TPD_50', 'R_TPD_60', 'R_R_10', 'R_R_20', 'R_R_30', 'R_R_40', 'R_R_50', 'R_R_60', 'Rdyn_CA_10', 'Rdyn_CA_20', 'Rdyn_CA_30', 'Rdyn_CA_40', 'Rdyn_CA_50', 'Rdyn_CA_60', 'Rdyn_CD_10', 'Rdyn_CD_20', 'Rdyn_CD_30', 'Rdyn_CD_40', 'Rdyn_CD_50', 'Rdyn_CD_60', 'Rdyn_MA_10', 'Rdyn_MA_20', 'Rdyn_MA_30', 'Rdyn_MA_40', 'Rdyn_MA_50', 'Rdyn_MA_60', 'Rdyn_MD_10', 'Rdyn_MD_20', 'Rdyn_MD_30', 'Rdyn_MD_40', 'Rdyn_MD_50', 'Rdyn_MD_60', 'Rdyn_TPA_10', 'Rdyn_TPA_20', 'Rdyn_TPA_30', 'Rdyn_TPA_40', 'Rdyn_TPA_50', 'Rdyn_TPA_60', 'Rdyn_TPD_10', 'Rdyn_TPD_20', 'Rdyn_TPD_30', 'Rdyn_TPD_40', 'Rdyn_TPD_50', 'Rdyn_TPD_60', 'Rdyn_R_10', 'Rdyn_R_20', 'Rdyn_R_30', 'Rdyn_R_40', 'Rdyn_R_50', 'Rdyn_R_60', 'invaded']
    else:
        header_names = ['year','habitat','elevation','season','S', 'L', 'L/S','C', 'T', 'B', 'I', 'Ca', 'Loop', 'NCycles', 'O', 'T-B', 'T-I', 'I-I', 'I-B', 'GenSD', 'VulSD', 'MxSim', 'MaxChainLength', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'MeanShortestPathLength', 'complexity', 'components', 'cc', 'cc_rnd', 'diameter', 'diameter_rnd', 'radius', 'radius_rnd', 'compartmentalisation', 'modularity_rgraph', 'modules_rgraph', 'modularity_rgraph_randoms', 'sd_modularity_rgraph_randoms', 'chi_square', 'p_value', 'mean_mass_ratio', 'min_biomass', 'max_biomass', 'mean_biomass', 'mean_ratio', 'variance_ratio', 'skew_ratio', 'kurt_ratio', 'R_removed', 'removed_species', 'invaded']
    out = csv.DictWriter(open('../output/output.csv', 'w'), header_names, delimiter=',')
    
    #we write the headers
    headers_dict = dict()
    for n in header_names:
        headers_dict[n] = n
    
    out.writerow(headers_dict)
    
    data_for_plots = dict()
    
    invasion = INVASION
    if invasion:
        #output file for the invaded networks
        if CALCULATE_ROBUSTNESS:
            header_names = ['year','habitat','elevation','season','invader','S', 'L', 'L/S','C', 'T', 'B', 'I', 'Ca', 'Loop', 'NCycles', 'O', 'T-B', 'T-I', 'I-I', 'I-B', 'GenSD', 'VulSD', 'MxSim', 'MaxChainLength', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'MeanShortestPathLength', 'complexity', 'components', 'cc', 'cc_rnd', 'diameter', 'diameter_rnd', 'radius', 'radius_rnd', 'compartmentalisation', 'modularity_rgraph', 'modules_rgraph', 'modularity_rgraph_randoms', 'sd_modularity_rgraph_randoms', 'chi_square', 'p_value', 'mean_mass_ratio', 'min_biomass', 'max_biomass', 'mean_biomass', 'mean_ratio', 'variance_ratio', 'skew_ratio', 'kurt_ratio', 'R_removed', 'removed_species', 'R_CA_10', 'R_CA_20', 'R_CA_30', 'R_CA_40', 'R_CA_50', 'R_CA_60', 'R_CD_10', 'R_CD_20', 'R_CD_30', 'R_CD_40', 'R_CD_50', 'R_CD_60', 'R_MA_10', 'R_MA_20', 'R_MA_30', 'R_MA_40', 'R_MA_50', 'R_MA_60', 'R_MD_10', 'R_MD_20', 'R_MD_30', 'R_MD_40', 'R_MD_50', 'R_MD_60', 'R_TPA_10', 'R_TPA_20', 'R_TPA_30', 'R_TPA_40', 'R_TPA_50', 'R_TPA_60', 'R_TPD_10', 'R_TPD_20', 'R_TPD_30', 'R_TPD_40', 'R_TPD_50', 'R_TPD_60', 'R_R_10', 'R_R_20', 'R_R_30', 'R_R_40', 'R_R_50', 'R_R_60', 'Rdyn_CA_10', 'Rdyn_CA_20', 'Rdyn_CA_30', 'Rdyn_CA_40', 'Rdyn_CA_50', 'Rdyn_CA_60', 'Rdyn_CD_10', 'Rdyn_CD_20', 'Rdyn_CD_30', 'Rdyn_CD_40', 'Rdyn_CD_50', 'Rdyn_CD_60', 'Rdyn_MA_10', 'Rdyn_MA_20', 'Rdyn_MA_30', 'Rdyn_MA_40', 'Rdyn_MA_50', 'Rdyn_MA_60', 'Rdyn_MD_10', 'Rdyn_MD_20', 'Rdyn_MD_30', 'Rdyn_MD_40', 'Rdyn_MD_50', 'Rdyn_MD_60', 'Rdyn_TPA_10', 'Rdyn_TPA_20', 'Rdyn_TPA_30', 'Rdyn_TPA_40', 'Rdyn_TPA_50', 'Rdyn_TPA_60', 'Rdyn_TPD_10', 'Rdyn_TPD_20', 'Rdyn_TPD_30', 'Rdyn_TPD_40', 'Rdyn_TPD_50', 'Rdyn_TPD_60', 'Rdyn_R_10', 'Rdyn_R_20', 'Rdyn_R_30', 'Rdyn_R_40', 'Rdyn_R_50', 'Rdyn_R_60', 'invaded']
        else:
            header_names = ['year','habitat','elevation','season','invader','S', 'L', 'L/S','C', 'T', 'B', 'I', 'Ca', 'Loop', 'NCycles', 'O', 'T-B', 'T-I', 'I-I', 'I-B', 'GenSD', 'VulSD', 'MxSim', 'MaxChainLength', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'MeanShortestPathLength', 'complexity', 'components', 'cc', 'cc_rnd', 'diameter', 'diameter_rnd', 'radius', 'radius_rnd', 'compartmentalisation', 'modularity_rgraph', 'modules_rgraph', 'modularity_rgraph_randoms', 'sd_modularity_rgraph_randoms', 'chi_square', 'p_value', 'mean_mass_ratio', 'min_biomass', 'max_biomass', 'mean_biomass', 'mean_ratio', 'variance_ratio', 'skew_ratio', 'kurt_ratio', 'R_removed', 'removed_species', 'invaded']
        
        out_invasive = csv.DictWriter(open('../output/output_invasive.csv', 'w'), header_names, delimiter=',')
        
        #we write the headers
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
        
        out_invasive.writerow(headers_dict)
    
    #for the 'per node' values
    header_names = ['net_name', 'node_name', 'species_mass(g)', 'indegree', 'outdegree', 'norm_indegree', 'norm_outdegree', 'betweeness', 'deg_centrality', 'cc', 'TL', 'TP', 'NumberOfPaths', 'O', 'module_rgraph', 'role_rgraph']
    out_per_node = csv.DictWriter(open('../output/output_per_node.csv', 'w'), header_names, delimiter=',')
    
    headers_dict = dict()
    for n in header_names:
        headers_dict[n] = n
        
    out_per_node.writerow(headers_dict)
    
    #for the 'per edge' values
    header_names = ['net_name', 'predator', 'predator_mass(g)', 'prey', 'prey_mass(g)', 'strength', 'mass_ratio']
    out_per_edge = csv.DictWriter(open('../output/output_per_edge.csv', 'w'), header_names, delimiter=',')
    
    headers_dict = dict()
    for n in header_names:
        headers_dict[n] = n
        
    out_per_edge.writerow(headers_dict)
    
    if REMOVAL:
        dr = DatabaseReader()
    #species_mass = dr.read_species_mass()
    
    out_row = dict()
    out_row_node = dict()
    out_row_edge = dict()
    
    if invasion:
        out_row_invasive = dict()
    
    if INVASION_SCENARIOS:
        threads = []
    
    years = sorted(nets_dict.keys())
    for year in years:
        invasion = INVASION
        out_row['year'] = year
        
        data_for_plots[year] = dict()
        
        if invasion:
            out_row_invasive['year'] = year
        
        for habitat in nets_dict[year].keys():
            invasion = INVASION
            out_row['habitat'] = habitat
            
            data_for_plots[year][habitat] = dict()
            
            if invasion:
                out_row_invasive['habitat'] = habitat
                
            for elevation in nets_dict[year][habitat].keys():
                invasion = INVASION
                out_row['elevation'] = elevation
                
                data_for_plots[year][habitat][elevation] = dict()
                
                if invasion:
                    out_row_invasive['elevation'] = elevation
                    
                for season in nets_dict[year][habitat][elevation].keys():
                    
                    if not season in CALCULATE_SEASONS:
                        continue
                    
                    invasion = INVASION
                    out_row['season'] = season
                    
                    net_original = nets_dict[year][habitat][elevation][season]
                    net_original.longest_path_length()
                    net = net_original.get_copy_removed_nodes(REMOVED_NODES)
                    
                    #net = net_original.get_aggregated_network(0.75)
                    
                    if net.size() == 0 or net.number_of_nodes() == 0:
                        continue
                    
                    if INVASION_SCENARIOS:
                        inv_thread = ThreadInvadedNets(net, year, habitat, elevation, season)
                        threads.append(inv_thread)
                        inv_thread.start()
                        #generate_invaded_networks(net, year, habitat, elevation, season)
                    
                    if invasion:
                        try:
                            net_invasive = net.get_invaded_network(introduced_group=INTRODUCED_GROUP, invasive_name=INTRODUCED_SPECIES, biomass_introduced=BIOMASS_INTRODUCED, generalism=FRACTION_OF_GENERALISM, predators=FRACTION_OF_PREDATORS, ext_threshold=INVASION_EXT_THRESHOLD)
                        except NotInvadableNetwork:
                            invasion = False
                    
                    if invasion:
                        out_row['invaded'] = 'True'
                        out_row_invasive['season'] = season                        
                        net_invasive.longest_path_length()
                        net_inv = net_invasive.get_copy_removed_nodes(REMOVED_NODES)
                        
                        try:
                            inv_name = INTRODUCED_GROUP + ' - ' +str(net_invasive.node[INTRODUCED_SPECIES]['biomass'])+' g'
                            #inv_name = 1
                        except:
                            inv_name = 'No invasion'
                            continue
                        
                        out_row_invasive['invader'] = inv_name
                    
                    else:
                        out_row['invaded'] = 'False'
                    
                    out_row['S'] = net.number_of_nodes()
                    out_row['L'] = net.size()
                    ld = net.linkage_density()
                    out_row['L/S'] = ld
                    c = net.connectance()
                    out_row['C'] = c
                    out_row['T'], top_sps = net.top_predators()
                    out_row['B'], basal_sps = net.basal(heterotrophs=False)
                    out_row['I'], inter_sps = net.intermediate(heterotrophs=False)
                    
                    fractions = net.get_links_fractions_between_levels(basal_sps, inter_sps, top_sps)
                    out_row['T-B'] = fractions['tb']
                    out_row['T-I'] = fractions['ti']
                    out_row['I-I'] = fractions['ii']
                    out_row['I-B'] = fractions['ib']
                    
                    gen_sd, vul_sd = net.generality_vulnerability_sd()
                    out_row['GenSD'] = gen_sd
                    out_row['VulSD'] = vul_sd
                    
                    out_row['MxSim'] = net.maximum_similarity()
                    
                    out_row['Ca'] = net.cannibalism()
                    out_row['Loop'], out_row['NCycles'] = net.fraction_in_loops()
                    out_row['MaxChainLength'] = net.longest_path_length()
                    out_row['MeanShortestPathLength'] = net.mean_path_length()
                    out_row['O'], omni_sps = net_original.omnivory()
                    out_row['complexity'] = net.complexity()
                    out_row['components'] = net.components()
                    
                    if invasion:
                        out_row_invasive['S'] = net_inv.number_of_nodes()
                        out_row_invasive['L'] = net_inv.size()
                        out_row_invasive['L/S'] = net_inv.linkage_density()
                        out_row_invasive['C'] = net_inv.connectance()
                        out_row_invasive['T'], top_sps = net_inv.top_predators()
                        out_row_invasive['B'], basal_sps = net_inv.basal(heterotrophs=False)
                        out_row_invasive['I'], inter_sps = net_inv.intermediate(heterotrophs=False)
                        
                        fractions = net_inv.get_links_fractions_between_levels(basal_sps, inter_sps, top_sps)
                        out_row_invasive['T-B'] = fractions['tb']
                        out_row_invasive['T-I'] = fractions['ti']
                        out_row_invasive['I-I'] = fractions['ii']
                        out_row_invasive['I-B'] = fractions['ib']
                        
                        gen_sd, vul_sd = net_inv.generality_vulnerability_sd()
                        out_row_invasive['GenSD'] = gen_sd
                        out_row_invasive['VulSD'] = vul_sd
                        
                        out_row_invasive['MxSim'] = net_inv.maximum_similarity()
                        
                        out_row_invasive['Ca'] = net_inv.cannibalism()
                        out_row_invasive['Loop'], out_row_invasive['NCycles'] = net_inv.fraction_in_loops()
                        out_row_invasive['MaxChainLength'] = net_inv.longest_path_length()
                        out_row_invasive['MeanShortestPathLength'] = net_inv.mean_path_length()
                        out_row_invasive['O'], omni_sps = net_invasive.omnivory()
                        out_row_invasive['complexity'] = net_inv.complexity()
                        out_row_invasive['components'] = net_inv.components()
                        
                    
                    net_und = net.to_undirected()
                    #clustering coefficient for the network
                    cc = 0.0
                    cc_rnd = 0.0
                    diameter = 0
                    diameter_rnd = 0
                    radius = 0
                    radius_rnd = 0
                    if net_und.number_of_nodes() > 0:
                        net_rnd = net.generate_random_graph()
                        net_rnd_und = net_rnd.to_undirected()
                        
                        cc = nx.average_clustering(net_und)
                        cc_rnd = nx.average_clustering(net_rnd_und)
                        
                        #we calculate structural measures of diameter and radius for the original network
                        try:
                            eccs_net = nx.eccentricity(net_und)
                            diameter = max(eccs_net.values())
                            radius = min(eccs_net.values())
                        except:
                            pass
                        #we calculate structural measures of diameter and radius for the random network
                        try:
                            eccs_net_rnd = nx.eccentricity(net_rnd_und)
                            diameter_rnd = max(eccs_net_rnd.values())
                            radius_rnd = min(eccs_net_rnd.values())
                        except:
                            pass
                            
                    out_row['cc'] = cc
                    out_row['cc_rnd'] = cc_rnd
                    out_row['diameter'] = diameter
                    out_row['diameter_rnd'] = diameter_rnd
                    out_row['radius'] = radius
                    out_row['radius_rnd'] = radius_rnd
                    
                    tps, nops, mean_length = net.find_trophic_positions()
                    out_row['MeanFoodChainLength'] = mean_length
                    
                    length_var, length_sd, paths_no = net.get_path_length_feats()
                    out_row['ChnSD'] = length_sd
                    out_row['ChnNo'] = paths_no
                    
                    out_row['compartmentalisation'] = net.degree_of_compartmentalization()
                    
                    if RGRAPH_MOD:
                        if net.order() < 7:
                            randoms = 0
                        else:
                            randoms = RGRAPH_RANDOMIZATIONS 
        
                        modularity, no_modules = net.modularity_rgraph(seed=RGRAPH_SEED, iter_factor=RGRAPH_ITERATION_F, cooling_factor=RGRAPH_COOLING_F, randoms=randoms)
                        
                        print 'modularity done'
                        
                        out_row['modularity_rgraph'] = modularity
                        out_row['modules_rgraph'] = no_modules
                        out_row['modularity_rgraph_randoms'] = net.modularity_of_randomizations
                        out_row['sd_modularity_rgraph_randoms'] = net.sd_mod_of_randomizations
                    else:
                        out_row['modularity_rgraph'] = 'N/A'
                        out_row['modules_rgraph'] = 'N/A'
                        out_row['modularity_rgraph_randoms'] = 'N/A'
                        out_row['sd_modularity_rgraph_randoms'] = 'N/A'
                    
                    chi_sq, p_val = net.chi_square_deg_freq()
                    out_row['chi_square'] = chi_sq
                    out_row['p_value'] = p_val
                    
                    out_row['mean_mass_ratio'] = net.mean_body_mass_ratio()
                    
                    min_bm, max_bm, mean_bm = net.get_min_max_biomasses()
                    out_row['min_biomass'] = min_bm
                    out_row['max_biomass'] = max_bm
                    out_row['mean_biomass'] = mean_bm
                    
                    if REMOVAL:
                        sps_to_remove = dr.get_species_to_remove()
                        if len(sps_to_remove) == 0:
                            sps_to_remove = set(REMOVED_SPECIES)
                        if len(sps_to_remove) == 0:
                            out_row['R_removed'] = 'N/A'
                            out_row['removed_species'] = 'N/A'
                        else:
                            n_removed, r_removed = net.robustness_to_removal(to_remove=sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
                            out_row['R_removed'] = r_removed
                            out_row['removed_species'] = sorted(set(net.nodes()) - set(n_removed.nodes()))
                    else:
                        out_row['R_removed'] = 'N/A'
                        out_row['removed_species'] = 'N/A'
                          
                    if season == 'ss':
                        data_for_plots[year][habitat][elevation]['species_number'] = net.number_of_nodes()
                        data_for_plots[year][habitat][elevation]['min_biomass'] = min_bm
                        data_for_plots[year][habitat][elevation]['max_biomass'] = max_bm
                        data_for_plots[year][habitat][elevation]['mean_biomass'] = mean_bm
                        data_for_plots[year][habitat][elevation]['mean_fcl'] = mean_length
                        data_for_plots[year][habitat][elevation]['connectance'] = c
                        data_for_plots[year][habitat][elevation]['linkage_density'] = ld
                        
                    ##now we calculate the above metrics for the invaded network
                    if invasion:
                        net_und_inv = net_inv.to_undirected()
                        cc = 0.0
                        cc_rnd = 0.0
                        diameter = 0
                        diameter_rnd = 0
                        radius = 0
                        radius_rnd = 0
                        if net_und_inv.number_of_nodes() > 0:
                            net_rnd = net_inv.generate_random_graph()
                            net_rnd_und = net_rnd.to_undirected()
                            
                            cc = nx.average_clustering(net_und_inv)
                            cc_rnd = nx.average_clustering(net_rnd_und)
                            
                            #we calculate structural measures of diameter and radius for the original network
                            try:
                                eccs_net = nx.eccentricity(net_und_inv)
                                diameter = max(eccs_net.values())
                                radius = min(eccs_net.values())
                            except:
                                pass
                            #we calculate structural measures of diameter and radius for the random network
                            try:
                                eccs_net_rnd = nx.eccentricity(net_rnd_und)
                                diameter_rnd = max(eccs_net_rnd.values())
                                radius_rnd = min(eccs_net_rnd.values())
                            except:
                                pass
                                
                        out_row_invasive['cc'] = cc
                        out_row_invasive['cc_rnd'] = cc_rnd
                        out_row_invasive['diameter'] = diameter
                        out_row_invasive['diameter_rnd'] = diameter_rnd
                        out_row_invasive['radius'] = radius
                        out_row_invasive['radius_rnd'] = radius_rnd
                        
                        tps_inv, nops_inv, mean_length_inv = net_inv.find_trophic_positions()
                        out_row_invasive['MeanFoodChainLength'] = mean_length_inv
                        
                        length_var_inv, length_sd_inv, paths_no_inv = net_inv.get_path_length_feats()
                        out_row_invasive['ChnSD'] = length_sd_inv
                        out_row_invasive['ChnNo'] = paths_no_inv
                        
                        out_row_invasive['compartmentalisation'] = net_inv.degree_of_compartmentalization()
                        
                        if RGRAPH_MOD:
                            if net_inv.order() < 7:
                                randoms = 0
                            else:
                                randoms = RGRAPH_RANDOMIZATIONS 
        
                            modularity, no_modules = net_inv.modularity_rgraph(seed=RGRAPH_SEED, iter_factor=RGRAPH_ITERATION_F, cooling_factor=RGRAPH_COOLING_F, randoms=randoms)
                            out_row_invasive['modularity_rgraph'] = modularity
                            out_row_invasive['modules_rgraph'] = no_modules
                            out_row_invasive['modularity_rgraph_randoms'] = net_inv.modularity_of_randomizations
                            out_row_invasive['sd_modularity_rgraph_randoms'] = net_inv.sd_mod_of_randomizations
                        else:
                            out_row_invasive['modularity_rgraph'] = 'N/A'
                            out_row_invasive['modules_rgraph'] = 'N/A'
                            out_row_invasive['modularity_rgraph_randoms'] = 'N/A'
                            out_row_invasive['sd_modularity_rgraph_randoms'] = 'N/A'
                            
                        chi_sq_inv, p_val_inv = net_inv.chi_square_deg_freq()
                        out_row_invasive['chi_square'] = chi_sq_inv
                        out_row_invasive['p_value'] = p_val_inv
                        
                        out_row_invasive['mean_mass_ratio'] = net_inv.mean_body_mass_ratio()
                        
                        min_bm_inv, max_bm_inv, mean_bm_inv = net_inv.get_min_max_biomasses()
                        out_row_invasive['min_biomass'] = min_bm_inv
                        out_row_invasive['max_biomass'] = max_bm_inv
                        out_row_invasive['mean_biomass'] = mean_bm_inv
                        
                        if REMOVAL:
                            sps_to_remove = dr.get_species_to_remove()
                            if len(sps_to_remove) == 0:
                                sps_to_remove = set(REMOVED_SPECIES)
                            if len(sps_to_remove) == 0:
                                out_row_invasive['R_removed'] = 'N/A'
                                out_row_invasive['removed_species'] = 'N/A'
                            else:
                                n_removed, r_removed = net_inv.robustness_to_removal(to_remove=sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
                                out_row_invasive['R_removed'] = r_removed
                                out_row_invasive['removed_species'] = sorted(set(net_inv.nodes()) - set(n_removed.nodes()))
                        else:
                            out_row_invasive['R_removed'] = 'N/A'
                            out_row_invasive['removed_species'] = 'N/A'
                    
                    if net.order() <= 2:
                        continue
                    
                    if CALCULATE_ROBUSTNESS:
                        rob_conn_a = net.robustness(criterion='conn', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        rob_conn_d = net.robustness(criterion='conn', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        rob_mass_a = net.robustness(criterion='mass', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        rob_mass_d = net.robustness(criterion='mass', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        rob_tp_a = net.robustness(criterion='trophic_position', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        rob_tp_d = net.robustness(criterion='trophic_position', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        rob_random = net.robustness(criterion='random', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                        
                        rob_conn_a_dyn = net.robustness(criterion='conn', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        rob_conn_d_dyn = net.robustness(criterion='conn', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        rob_mass_a_dyn = net.robustness(criterion='mass', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        rob_mass_d_dyn = net.robustness(criterion='mass', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        rob_tp_a_dyn = net.robustness(criterion='trophic_position', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        rob_tp_d_dyn = net.robustness(criterion='trophic_position', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        rob_random_dyn = net.robustness(criterion='random', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                        
                        fractions = [10, 20, 30, 40, 50, 60]
                        for fr in fractions:
                            value = 'N/A'
                            key = 'R_CA_'+str(fr)
                            if rob_conn_a.has_key(fr):
                                value = rob_conn_a[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'R_CD_'+str(fr)
                            if rob_conn_d.has_key(fr):
                                value = rob_conn_d[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'R_MA_'+str(fr)
                            if rob_mass_a.has_key(fr):
                                value = rob_mass_a[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'R_MD_'+str(fr)
                            if rob_mass_d.has_key(fr):
                                value = rob_mass_d[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'R_TPA_'+str(fr)
                            if rob_tp_a.has_key(fr):
                                value = rob_tp_a[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'R_TPD_'+str(fr)
                            if rob_tp_d.has_key(fr):
                                value = rob_tp_d[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'R_R_'+str(fr)
                            if rob_random.has_key(fr):
                                value = rob_random[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_CA_'+str(fr)
                            if rob_conn_a_dyn.has_key(fr):
                                value = rob_conn_a_dyn[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_CD_'+str(fr)
                            if rob_conn_d_dyn.has_key(fr):
                                value = rob_conn_d_dyn[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_MA_'+str(fr)
                            if rob_mass_a_dyn.has_key(fr):
                                value = rob_mass_a_dyn[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_MD_'+str(fr)
                            if rob_mass_d_dyn.has_key(fr):
                                value = rob_mass_d_dyn[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_TPA_'+str(fr)
                            if rob_tp_a_dyn.has_key(fr):
                                value = rob_tp_a_dyn[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_TPD_'+str(fr)
                            if rob_tp_d_dyn.has_key(fr):
                                value = rob_tp_d_dyn[fr]
                            
                            out_row[key] = value
                            
                            value = 'N/A'
                            key = 'Rdyn_R_'+str(fr)
                            if rob_random_dyn.has_key(fr):
                                value = rob_random_dyn[fr]
                            
                            out_row[key] = value
                    
                    if invasion:
                        if CALCULATE_ROBUSTNESS:
                            rob_conn_a = net_inv.robustness(criterion='conn', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            rob_conn_d = net_inv.robustness(criterion='conn', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            rob_mass_a = net_inv.robustness(criterion='mass', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            rob_mass_d = net_inv.robustness(criterion='mass', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            rob_tp_a = net_inv.robustness(criterion='trophic_position', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            rob_tp_d = net_inv.robustness(criterion='trophic_position', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            rob_random = net_inv.robustness(criterion='random', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                            
                            rob_conn_a_dyn = net_inv.robustness(criterion='conn', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            rob_conn_d_dyn = net_inv.robustness(criterion='conn', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            rob_mass_a_dyn = net_inv.robustness(criterion='mass', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            rob_mass_d_dyn = net_inv.robustness(criterion='mass', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            rob_tp_a_dyn = net_inv.robustness(criterion='trophic_position', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            rob_tp_d_dyn = net_inv.robustness(criterion='trophic_position', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            rob_random_dyn = net_inv.robustness(criterion='random', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                            
                            for fr in fractions:
                                #robustness by connectance
                                #ascending
                                value = 'N/A'
                                key = 'R_CA_'+str(fr)
                                if rob_conn_a.has_key(fr):
                                    value = rob_conn_a[fr]
                                
                                out_row_invasive[key] = value
                                
                                #descending
                                value = 'N/A'
                                key = 'R_CD_'+str(fr)
                                if rob_conn_d.has_key(fr):
                                    value = rob_conn_d[fr]
                                
                                out_row_invasive[key] = value
                                
                                #robustness by biomass
                                #ascending
                                value = 'N/A'
                                key = 'R_MA_'+str(fr)
                                if rob_mass_a.has_key(fr):
                                    value = rob_mass_a[fr]
                                
                                out_row_invasive[key] = value
                                
                                #descending
                                value = 'N/A'
                                key = 'R_MD_'+str(fr)
                                if rob_mass_d.has_key(fr):
                                    value = rob_mass_d[fr]
                                
                                out_row_invasive[key] = value
                                
                                #robustness by trophic position
                                #ascending
                                value = 'N/A'
                                key = 'R_TPA_'+str(fr)
                                if rob_tp_a.has_key(fr):
                                    value = rob_tp_a[fr]
                                
                                out_row_invasive[key] = value
                                
                                #descending
                                value = 'N/A'
                                key = 'R_TPD_'+str(fr)
                                if rob_tp_d.has_key(fr):
                                    value = rob_tp_d[fr]
                                
                                out_row_invasive[key] = value
                                
                                #robustness random
                                value = 'N/A'
                                key = 'R_R_'+str(fr)
                                if rob_random.has_key(fr):
                                    value = rob_random[fr]
                                
                                out_row_invasive[key] = value
                                
                                #dynamic robustness by connectance
                                #ascending
                                value = 'N/A'
                                key = 'Rdyn_CA_'+str(fr)
                                if rob_conn_a_dyn.has_key(fr):
                                    value = rob_conn_a_dyn[fr]
                                
                                out_row_invasive[key] = value
                                
                                #descending
                                value = 'N/A'
                                key = 'Rdyn_CD_'+str(fr)
                                if rob_conn_d_dyn.has_key(fr):
                                    value = rob_conn_d_dyn[fr]
                                
                                out_row_invasive[key] = value
                                
                                #dynamic robustness by biomass
                                #ascending
                                value = 'N/A'
                                key = 'Rdyn_MA_'+str(fr)
                                if rob_mass_a_dyn.has_key(fr):
                                    value = rob_mass_a_dyn[fr]
                                
                                out_row_invasive[key] = value
                                
                                #descending
                                value = 'N/A'
                                key = 'Rdyn_MD_'+str(fr)
                                if rob_mass_d_dyn.has_key(fr):
                                    value = rob_mass_d_dyn[fr]
                                
                                out_row_invasive[key] = value
                                
                                #dynamic robustness by trophic position
                                #ascending
                                value = 'N/A'
                                key = 'Rdyn_TPA_'+str(fr)
                                if rob_tp_a_dyn.has_key(fr):
                                    value = rob_tp_a_dyn[fr]
                                
                                out_row_invasive[key] = value
                                
                                #descending
                                value = 'N/A'
                                key = 'Rdyn_TPD_'+str(fr)
                                if rob_tp_d_dyn.has_key(fr):
                                    value = rob_tp_d_dyn[fr]
                                
                                #dynamic robustness random
                                out_row_invasive[key] = value
                                value = 'N/A'
                                key = 'Rdyn_R_'+str(fr)
                                if rob_random_dyn.has_key(fr):
                                    value = rob_random_dyn[fr]
                                
                                out_row_invasive[key] = value
                            
                        new_row = out_row.copy()
                        new_row['invader'] = 0
                        
                        new_row['mean_ratio'] = ''
                        new_row['variance_ratio'] = ''
                        new_row['skew_ratio'] = ''
                        new_row['kurt_ratio'] = ''
                        
                        out_invasive.writerow(new_row)
                        out_invasive.writerow(out_row_invasive)
                    
                     
                    indegrees = net.in_degree()
                    outdegrees = net.out_degree()
                    total_in = sum(indegrees.values())
                    total_out = sum(outdegrees.values())
                    betweeness = nx.algorithms.centrality.betweenness_centrality(net)
                    deg_centrality = nx.algorithms.centrality.degree_centrality(net)
                    
                    #clustering coefficients
                    ccs = nx.clustering(net_und)
                    
                    out_row_node['net_name'] = str(year)+'-'+str(habitat)+'-'+str(elevation)+'-'+str(season)
                    
                    links = net.size()
                    for vertex in net.nodes():
                        out_row_node['node_name'] = vertex
                        
                        #if species_mass.has_key(vertex):
                        m = net.node[vertex]['biomass']
                        if m == ' ':
                            out_row_node['species_mass(g)'] = 'N/A'
                        else:
                            out_row_node['species_mass(g)'] = m
                        #else:
                        #    out_row_node['species_mass(g)'] = 'N/A'
                        
                        indeg = indegrees[vertex]
                        out_row_node['indegree'] = indeg
                        
                        outdeg = outdegrees[vertex]
                        out_row_node['outdegree'] = outdeg
                        
                        out_row_node['norm_indegree'] = float(indeg)/links
                        out_row_node['norm_outdegree'] = float(outdeg)/links
                        out_row_node['betweeness'] = betweeness[vertex]
                        out_row_node['deg_centrality'] = deg_centrality[vertex]
                        out_row_node['cc'] = ccs[vertex]
                        
                        if vertex in basal_sps or indeg == 0:
                            out_row_node['TL'] = 'B'
                        elif vertex in top_sps:
                            out_row_node['TL'] = 'T'
                        elif vertex in inter_sps:
                            out_row_node['TL'] = 'I'
                        else: 
                            out_row_node['TL'] = 'N/A'
                        
                        if vertex in omni_sps:
                            out_row_node['O'] = '1'
                        else:
                            out_row_node['O'] = '0'
                        
                        
                        out_row_node['TP'] = tps[vertex]
                        out_row_node['NumberOfPaths'] = nops[vertex]
                        
                        if RGRAPH_MOD and net.node[vertex].has_key('module'):
                            out_row_node['module_rgraph'] = net.node[vertex]['module']
                            out_row_node['role_rgraph'] = net.node[vertex]['role']
                        else:
                            out_row_node['module_rgraph'] = 'N/A'
                            out_row_node['role_rgraph'] = 'N/A'
                        
                        out_per_node.writerow(out_row_node)
                    
                    net.obtain_interactions_strengths()
                    out_row_edge['net_name'] = str(year)+'-'+str(habitat)+'-'+str(elevation)+'-'+str(season)
                    
                    ratio_dist = []
                    
                    for prey, predator, atts in net.edges(data=True):
                        #header_names = ['net_name', 'predator', 'predator_mass(g)', 'prey', 'prey_mass(g)', 'strength']
                        out_row_edge['predator'] = predator
                        out_row_edge['predator_mass(g)'] = net.node[predator]['biomass']
                        
                        out_row_edge['prey'] = prey
                        out_row_edge['prey_mass(g)'] = net.node[prey]['biomass']
                        
                        out_row_edge['strength'] = atts['weight']
                        
                        if atts.has_key('mass_ratio'):
                            out_row_edge['mass_ratio'] = atts['mass_ratio']
                            ratio_dist.append(math.log10(atts['mass_ratio']))
                        else:
                            out_row_edge['mass_ratio'] = 'N/A'
                        
                        out_per_edge.writerow(out_row_edge)
                    
                    
                    out_row['mean_ratio'] = scipy.mean(ratio_dist)
                    out_row['variance_ratio'] = scipy.var(ratio_dist)
                    out_row['skew_ratio'] = scipy.stats.stats.skew(ratio_dist)
                    out_row['kurt_ratio'] = scipy.stats.stats.kurtosis(ratio_dist)
                    
                    out.writerow(out_row)
    
    if INVASION_SCENARIOS:
        for t in threads:
            t.join()
    
    return data_for_plots

def create_compare_output(name_a, network_a, name_b, network_b):
    #for the 'per node' values
    header_names = ['net_name', 'node_name', 'species_mass(g)', 'indegree', 'outdegree', 'norm_indegree', 'norm_outdegree', 'betweeness', 'deg_centrality', 'cc', 'TL', 'TP', 'NumberOfPaths', 'O']
    
    filename = name_a+'_vs_'+name_b+'.csv'
    out_vs = csv.DictWriter(open('../output/'+filename, 'w'), header_names, delimiter=',')
    
    headers_dict = dict()
    for n in header_names:
        headers_dict[n] = n
        
    out_vs.writerow(headers_dict)
    
    out_row_node_a = dict()
    out_row_node_b = dict()
    
    out_row_node_a['net_name'] = 'The species that are in both networks are: '
    out_vs.writerow(out_row_node_a)
    out_row_node_a.clear()
    
    sps_a = set(network_a.nodes())
    sps_b = set(network_b.nodes())
    both_species = sorted(sps_a & sps_b)
    
    net_a_und = network_a.to_undirected()
    indegrees_a = network_a.in_degree()
    outdegrees_a = network_a.out_degree()
    total_in_a = sum(indegrees_a.values())
    total_out_a = sum(outdegrees_a.values())
    betweeness_a = nx.algorithms.centrality.betweenness_centrality(network_a)
    deg_centrality_a = nx.algorithms.centrality.degree_centrality(network_a)
    
    #clustering coefficients
    ccs_a = nx.clustering(net_a_und)
    
    network_a.longest_path_length()
    x, top_sps_a = network_a.top_predators()
    x, basal_sps_a = network_a.basal(heterotrophs=True)
    x, inter_sps_a = network_a.intermediate(heterotrophs=True)
    x, omni_sps_a = network_a.omnivory()
    
    tps_a, nops_a, mean_length_a = network_a.find_trophic_positions()
    
    net_b_und = network_b.to_undirected()
    indegrees_b = network_b.in_degree()
    outdegrees_b = network_b.out_degree()
    total_in_b = sum(indegrees_b.values())
    total_out_b = sum(outdegrees_b.values())
    betweeness_b = nx.algorithms.centrality.betweenness_centrality(network_b)
    deg_centrality_b = nx.algorithms.centrality.degree_centrality(network_b)
    
    #clustering coefficients
    ccs_b = nx.clustering(net_b_und)
    
    network_b.longest_path_length()
    x, top_sps_b = network_b.top_predators()
    x, basal_sps_b = network_b.basal(heterotrophs=True)
    x, inter_sps_b = network_b.intermediate(heterotrophs=True)
    x, omni_sps_b = network_b.omnivory()
    
    tps_b, nops_b, mean_length_b = network_b.find_trophic_positions()
    
    #species that are in both networks
    for vertex in both_species:
        out_row_node_a['net_name'] = name_a
        out_row_node_b['net_name'] = name_b
    
        out_row_node_a['node_name'] = vertex
        out_row_node_b['node_name'] = vertex
        
        if not network_a.node[vertex].has_key('biomass'):
            out_row_node_a['species_mass(g)'] = 'N/A'
        else:
            out_row_node_a['species_mass(g)'] = network_a.node[vertex]['biomass']
        
        if not network_b.node[vertex].has_key('biomass'):
            out_row_node_b['species_mass(g)'] = 'N/A'
        else:
            out_row_node_b['species_mass(g)'] = network_b.node[vertex]['biomass']
        
        indeg_a = indegrees_a[vertex]
        out_row_node_a['indegree'] = indeg_a
        
        indeg_b = indegrees_b[vertex]
        out_row_node_b['indegree'] = indeg_b
        
        outdeg_a = outdegrees_a[vertex]
        out_row_node_a['outdegree'] = outdeg_a
        
        outdeg_b = outdegrees_b[vertex]
        out_row_node_b['outdegree'] = outdeg_b
        
        out_row_node_a['norm_indegree'] = float(indeg_a)/float(total_in_a)
        out_row_node_b['norm_indegree'] = float(indeg_b)/float(total_in_b)
        
        out_row_node_a['norm_outdegree'] = float(outdeg_a)/float(total_out_a)
        out_row_node_b['norm_outdegree'] = float(outdeg_b)/float(total_out_b)
        
        out_row_node_a['betweeness'] = betweeness_a[vertex]
        out_row_node_b['betweeness'] = betweeness_b[vertex]
        
        out_row_node_a['deg_centrality'] = deg_centrality_a[vertex]
        out_row_node_b['deg_centrality'] = deg_centrality_b[vertex]
        
        out_row_node_a['cc'] = ccs_a[vertex]
        out_row_node_b['cc'] = ccs_b[vertex]
        
        if vertex in basal_sps_a or indeg_a == 0:
            out_row_node_a['TL'] = 'B'
        elif vertex in top_sps_a:
            out_row_node_a['TL'] = 'T'
        elif vertex in inter_sps_a:
            out_row_node_a['TL'] = 'I'
        else: 
            out_row_node_a['TL'] = 'N/A'
        
        if vertex in omni_sps_a:
            out_row_node_a['O'] = '1'
        else:
            out_row_node_a['O'] = '0'
            
        
        if vertex in basal_sps_b or indeg_b == 0:
            out_row_node_b['TL'] = 'B'
        elif vertex in top_sps_b:
            out_row_node_b['TL'] = 'T'
        elif vertex in inter_sps_b:
            out_row_node_b['TL'] = 'I'
        else: 
            out_row_node_b['TL'] = 'N/A'
        
        if vertex in omni_sps_b:
            out_row_node_b['O'] = '1'
        else:
            out_row_node_b['O'] = '0'
        
        
        out_row_node_a['TP'] = tps_a[vertex]
        out_row_node_b['TP'] = tps_b[vertex]
        
        out_row_node_a['NumberOfPaths'] = nops_a[vertex]
        out_row_node_b['NumberOfPaths'] = nops_b[vertex]
        
        out_vs.writerow(out_row_node_a)
        out_vs.writerow(out_row_node_b)
        
        out_row_node_a.clear()
        out_row_node_b.clear()
    
    #species that are in network a and not in b
    exclusive_a = sorted(sps_a - sps_b)
    out_row_node_a['net_name'] = 'The species in network '+name_a+' that are not in network '+name_b+' are: '
    out_vs.writerow(out_row_node_a)
    out_row_node_a.clear()
    
    for vertex in exclusive_a:
        out_row_node_a['net_name'] = name_a
        out_row_node_a['node_name'] = vertex
            
        if not network_a.node[vertex].has_key('biomass'):
            out_row_node_a['species_mass(g)'] = 'N/A'
        else:
            out_row_node_a['species_mass(g)'] = network_a.node[vertex]['biomass']
        
        indeg_a = indegrees_a[vertex]
        out_row_node_a['indegree'] = indeg_a
        
        outdeg_a = outdegrees_a[vertex]
        out_row_node_a['outdegree'] = outdeg_a
        
        out_row_node_a['norm_indegree'] = float(indeg_a)/float(total_in_a)
        out_row_node_a['norm_outdegree'] = float(outdeg_a)/float(total_out_a)
        
        out_row_node_a['betweeness'] = betweeness_a[vertex]
        out_row_node_a['deg_centrality'] = deg_centrality_a[vertex]
        
        out_row_node_a['cc'] = ccs_a[vertex]
        
        if vertex in basal_sps_a or indeg_a == 0:
            out_row_node_a['TL'] = 'B'
        elif vertex in top_sps_a:
            out_row_node_a['TL'] = 'T'
        elif vertex in inter_sps_a:
            out_row_node_a['TL'] = 'I'
        else: 
            out_row_node_a['TL'] = 'N/A'
        
        if vertex in omni_sps_a:
            out_row_node_a['O'] = '1'
        else:
            out_row_node_a['O'] = '0'
        
        out_row_node_a['TP'] = tps_a[vertex]
        
        out_row_node_a['NumberOfPaths'] = nops_a[vertex]
        
        out_vs.writerow(out_row_node_a)
        
        out_row_node_a.clear()
        
    
    #species that are in network b and not in a
    exclusive_b = sorted(sps_b - sps_a)
    out_row_node_b['net_name'] = 'The species in network '+name_a+' that are not in network '+name_b+' are: '
    out_vs.writerow(out_row_node_b)
    out_row_node_b.clear()
    
    
    for vertex in exclusive_b:
        out_row_node_b['net_name'] = name_b
        out_row_node_b['node_name'] = vertex
            
        if not network_b.node[vertex].has_key('biomass'):
            out_row_node_b['species_mass(g)'] = 'N/A'
        else:
            out_row_node_b['species_mass(g)'] = network_b.node[vertex]['biomass']
        
        indeg_b = indegrees_b[vertex]
        out_row_node_b['indegree'] = indeg_b
        
        outdeg_b = outdegrees_b[vertex]
        out_row_node_b['outdegree'] = outdeg_b
        
        out_row_node_b['norm_indegree'] = float(indeg_b)/float(total_in_b)
        out_row_node_b['norm_outdegree'] = float(outdeg_b)/float(total_out_b)
        
        out_row_node_b['betweeness'] = betweeness_b[vertex]
        out_row_node_b['deg_centrality'] = deg_centrality_b[vertex]
        
        out_row_node_b['cc'] = ccs_b[vertex]
        
        if vertex in basal_sps_b or indeg_b == 0:
            out_row_node_b['TL'] = 'B'
        elif vertex in top_sps_b:
            out_row_node_b['TL'] = 'T'
        elif vertex in inter_sps_b:
            out_row_node_b['TL'] = 'I'
        else: 
            out_row_node_b['TL'] = 'N/A'
        
        if vertex in omni_sps_b:
            out_row_node_b['O'] = '1'
        else:
            out_row_node_b['O'] = '0'
        
        out_row_node_b['TP'] = tps_b[vertex]
        out_row_node_b['NumberOfPaths'] = nops_b[vertex]
        
        out_vs.writerow(out_row_node_b)
        
        out_row_node_b.clear()
    
    return
    

def degree_per_periods(nets_dict, masses, sp_habitats):
    
    pre_years = [1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992]
    post_years = [1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001]
    
    elevations = [0,1,2]
    in_degs = dict.fromkeys(elevations, None)
    out_degs = dict.fromkeys(elevations, None)
    trophic_levels = dict.fromkeys(elevations, None)
    trophic_positions = dict.fromkeys(elevations, None)
    ratios = dict.fromkeys(elevations, None)
    
    
    sps_in_pre = dict.fromkeys(elevations, None)
    sps_in_post = dict.fromkeys(elevations, None)
    
    ratios_in_pre = dict.fromkeys(elevations, None)
    ratios_in_post = dict.fromkeys(elevations, None)
    
    for e in elevations:
        in_degs[e] = dict()
        out_degs[e] = dict()
        trophic_levels[e] = dict()
        trophic_positions[e] = dict()
        ratios[e] = dict()
        
        sps_in_pre[e] = set()
        sps_in_post[e] = set()
        
        ratios_in_pre[e] = set()
        ratios_in_post[e] = set()
        
    mean_species = dict.fromkeys(elevations, 0)
    
    count_species = dict.fromkeys(['Bird', 'BirdPrey', 'Mammal', 'MammalCarn', 'Amphibian', 'Reptile', 'Fish', 'N/A'], None)
    for k in count_species.keys():
        count_species[k] = set()
    
    
    for y in pre_years:
        for e in elevations:
            net = nets_dict[y][1][e]['all'].get_copy_removed_nodes(REMOVED_NODES)
            
            mean_species[e] += net.order()
            
            tls = net.get_trophic_levels()
            tps, np, ml = net.find_trophic_positions()
            
            if not in_degs[e].has_key('pre'):
                in_degs[e]['pre'] = dict()

            if not out_degs[e].has_key('pre'):
                out_degs[e]['pre'] = dict()
                
            if not trophic_levels[e].has_key('pre'):
                trophic_levels[e]['pre'] = dict()
                
            if not trophic_positions[e].has_key('pre'):
                trophic_positions[e]['pre'] = dict()

            if not ratios[e].has_key('pre'):
                ratios[e]['pre'] = []

            for n in net.nodes():
                
                ###this is for counting the species of each taxonomic group
                count_species[net.node[n]['group']].add(n)
                
                if not in_degs[e]['pre'].has_key(n):
                    in_degs[e]['pre'][n] = dict.fromkeys(['original', 'standardised'], None)
                    in_degs[e]['pre'][n]['original'] = []
                    in_degs[e]['pre'][n]['standardised'] = []
                    
                    
                in_degs[e]['pre'][n]['original'].append(net.in_degree(n))
                in_degs[e]['pre'][n]['standardised'].append(float(net.in_degree(n))/net.size())
                
                if not out_degs[e]['pre'].has_key(n):
                    out_degs[e]['pre'][n] = dict.fromkeys(['original', 'standardised'], None)
                    out_degs[e]['pre'][n]['original'] = []
                    out_degs[e]['pre'][n]['standardised'] = []
                
                out_degs[e]['pre'][n]['original'].append(net.out_degree(n))
                out_degs[e]['pre'][n]['standardised'].append(float(net.out_degree(n))/net.size())
                
                if not trophic_levels[e]['pre'].has_key(n):
                    trophic_levels[e]['pre'][n] = []
                
                trophic_levels[e]['pre'][n].append(tls[n])
                
                if not trophic_positions[e]['pre'].has_key(n):
                    trophic_positions[e]['pre'][n] = []
                
                trophic_positions[e]['pre'][n].append(tps[n])
                
                sps_in_pre[e].add(n)
            
            net.mean_body_mass_ratio()
            for u,v in net.edges():
                if net[u][v].has_key('mass_ratio') and not str(u+v) in ratios_in_pre[e]:
                    ratios[e]['pre'].append([u,v,math.log10(net[u][v]['mass_ratio'])])
                    ratios_in_pre[e].add(str(u+v))
    
    for e in elevations:
        print e, mean_species[e]
             
    for y in post_years:
        for e in elevations:
            net = nets_dict[y][1][e]['all'].get_copy_removed_nodes(REMOVED_NODES)
            tls = net.get_trophic_levels()
            tps, np, ml = net.find_trophic_positions()
        
            if not in_degs[e].has_key('post'):
                in_degs[e]['post'] = dict()
            
            if not out_degs[e].has_key('post'):
                out_degs[e]['post'] = dict()
            
            if not trophic_levels[e].has_key('post'):
                trophic_levels[e]['post'] = dict()
                
            if not trophic_positions[e].has_key('post'):
                trophic_positions[e]['post'] = dict()
            
            if not ratios[e].has_key('post'):
                ratios[e]['post'] = []
            
            for n in net.nodes():
                if n in sps_in_pre[e]:
                    continue
                
                count_species[net.node[n]['group']].add(n)
                
                if not in_degs[e]['post'].has_key(n):
                    in_degs[e]['post'][n] = dict.fromkeys(['original', 'standardised'], None)
                    in_degs[e]['post'][n]['original'] = []
                    in_degs[e]['post'][n]['standardised'] = []
                
                in_degs[e]['post'][n]['original'].append(net.in_degree(n))
                in_degs[e]['post'][n]['standardised'].append(float(net.in_degree(n))/net.size())
                
                if not out_degs[e]['post'].has_key(n):
                    out_degs[e]['post'][n] = dict.fromkeys(['original', 'standardised'], None)
                    out_degs[e]['post'][n]['original'] = []
                    out_degs[e]['post'][n]['standardised'] = []
                
                out_degs[e]['post'][n]['original'].append(net.out_degree(n))
                out_degs[e]['post'][n]['standardised'].append(float(net.out_degree(n))/net.size())
                
                if not trophic_levels[e]['post'].has_key(n):
                    trophic_levels[e]['post'][n] = []
                
                trophic_levels[e]['post'][n].append(tls[n])
                
                if not trophic_positions[e]['post'].has_key(n):
                    trophic_positions[e]['post'][n] = []
                
                trophic_positions[e]['post'][n].append(tps[n])
                
                sps_in_post[e].add(n)
            
            net.mean_body_mass_ratio()
            for u,v in net.edges():
                if net[u][v].has_key('mass_ratio') and not str(u+v) in ratios_in_post[e] and not str(u+v) in ratios_in_pre[e]:
                    ratios[e]['post'].append([u,v,math.log10(net[u][v]['mass_ratio'])])
                    ratios_in_post[e].add(str(u+v))
    
    for k in count_species.keys():
        print k, len(count_species[k]), count_species[k]
    
    return sps_in_post
            
#    for e in elevations:
#        for sp in in_degs[e]['pre'].keys():
#            degs = in_degs[e]['pre'][sp]['original']
#            avg = float(sum(degs))/len(degs)
#            in_degs[e]['pre'][sp]['original'] = avg
#            
#            degs = in_degs[e]['pre'][sp]['standardised']
#            avg = float(sum(degs))/len(degs)
#            in_degs[e]['pre'][sp]['standardised'] = avg
#            
#            degs = out_degs[e]['pre'][sp]['original']
#            avg = float(sum(degs))/len(degs)
#            out_degs[e]['pre'][sp]['original'] = avg
#            
#            degs = out_degs[e]['pre'][sp]['standardised']
#            avg = float(sum(degs))/len(degs)
#            out_degs[e]['pre'][sp]['standardised'] = avg
#            
#            trophic_levels[e]['pre'][sp] = scipy.mean(trophic_levels[e]['pre'][sp])
#            trophic_positions[e]['pre'][sp] = scipy.mean(trophic_positions[e]['pre'][sp])
#            
#        for sp in in_degs[e]['post'].keys():
#            degs = in_degs[e]['post'][sp]['original']
#            avg = float(sum(degs))/len(degs)
#            in_degs[e]['post'][sp]['original'] = avg
#            
#            degs = in_degs[e]['post'][sp]['standardised']
#            avg = float(sum(degs))/len(degs)
#            in_degs[e]['post'][sp]['standardised'] = avg
#            
#            degs = out_degs[e]['post'][sp]['original']
#            avg = float(sum(degs))/len(degs)
#            out_degs[e]['post'][sp]['original'] = avg
#            
#            degs = out_degs[e]['post'][sp]['standardised']
#            avg = float(sum(degs))/len(degs)
#            out_degs[e]['post'][sp]['standardised'] = avg
#            
#            trophic_levels[e]['post'][sp] = scipy.mean(trophic_levels[e]['post'][sp])
#            trophic_positions[e]['post'][sp] = scipy.mean(trophic_positions[e]['post'][sp])
#    
#    headers = ['period', 'elevation', 'species', 'group', 'in_out', 'L', 'L_standardised', 'TL', 'TP', 'mass', 'prey_mass', 'coming_from_0', 'coming_from_1', 'coming_from_2', 'coming_from_3', 'n_habitats']
#    
#    out_file = open('../output/degree_distributions.csv', 'w')
#    out = csv.DictWriter(out_file, headers, delimiter=',')
#
#    out.writeheader()
#    
#    out_row_in = dict()
#    out_row_out = dict()
#    
#    
#    headers = ['elevation', 'species', 'group', 'degree', 'TL', 'TP', 'mass', 'n_habitats']
#    
#    out_file_newcomers = open('../output/newcomers.csv', 'w')
#    out_newcomers = csv.DictWriter(out_file_newcomers, headers, delimiter=',')
#
#    out_newcomers.writeheader()
#    
#    out_row_newcomers = dict()
#    
#    
#    headers = ['period', 'elevation', 'mean', 'variance', 'skewness', 'kurtosis']
#    
#    out_file_stats = open('../output/distributions_stats.csv', 'w')
#    out_stats = csv.DictWriter(out_file_stats, headers, delimiter=',')
#
#    out_stats.writeheader()
#    
#    out_row_st = dict()
#    
#    
#    headers = ['period', 'elevation', 'value', 'prey', 'prey_group', 'prey_tl', 'predator', 'predator_group', 'predator_tl']
#    
#    out_file_mass = open('../output/mass_ratio_distributions.csv', 'w')
#    out_mass = csv.DictWriter(out_file_mass, headers, delimiter=',')
#
#    out_mass.writeheader()
#    
#    out_row_mass = dict()
#    
#    out_row_in['in_out'] = 'in'
#    out_row_out['in_out'] = 'out'
#    
#    out_row_out['coming_from_0'] = ''
#    out_row_out['coming_from_1'] = ''
#    out_row_out['coming_from_2'] = ''
#    out_row_out['coming_from_3'] = ''
#    
#    out_row_out['n_habitats'] = ''
#    
#    for e in elevations:
#        out_row_in['elevation'] = e
#        out_row_out['elevation'] = e
#        out_row_st['elevation'] = e
#        out_row_newcomers['elevation'] = e
#        
#        out_row_in['period'] = 'pre'
#        out_row_out['period'] = 'pre'
#        out_row_st['period'] = 'pre'
#        
#        mass_dist = []
#        in_dist_plot = []
#        y = []
#        
#        for sp in in_degs[e]['pre'].keys():
#            out_row_in['species'] = sp
#            out_row_out['species'] = sp
#            
#            ind = in_degs[e]['pre'][sp]['original']
#            out_row_in['L'] = ind
#            out_row_in['L_standardised'] = in_degs[e]['pre'][sp]['standardised']
#            
#            out_row_in['TP'] = trophic_positions[e]['pre'][sp]
#            out_row_in['TL'] = trophic_levels[e]['pre'][sp]
#            
#            out_row_in['group'] = masses[sp]['group']
#            
#            if math.ceil(ind) > 0:
#                out_row_in['mass'] = masses[sp]['biomass']
#                biomass = float(masses[sp]['biomass'])
#                if biomass > 0.0:
#                    mass_dist.append(math.log10(biomass))
#                    in_dist_plot.append(ind)
#            else:
#                out_row_in['mass'] = ''
#            
#            
#            outd = out_degs[e]['pre'][sp]['original']
#            out_row_out['L'] = outd
#            out_row_out['L_standardised'] = out_degs[e]['pre'][sp]['standardised']
#            
#            if outd > 0.0:
#                out_row_in['prey_mass'] = masses[sp]['biomass']
#                mass = float(masses[sp]['biomass']) 
#                if  mass > 0.0:
#                    y.append(math.log10(mass))
#            else:
#                out_row_in['prey_mass'] = ''
#            
#            out_row_in['n_habitats'] = len(sp_habitats[sp])
#           
#            out.writerow(out_row_in)
#            out.writerow(out_row_out)
#        
#        out_row_st['mean'] = scipy.mean(mass_dist)
#        out_row_st['variance'] = scipy.var(mass_dist)
#        out_row_st['skewness'] = scipy.stats.stats.skew(mass_dist)
#        out_row_st['kurtosis'] = scipy.stats.stats.kurtosis(mass_dist)
#        
#        out_stats.writerow(out_row_st)
#        
#        out_row_mass['elevation'] = e
#        out_row_mass['period'] = 'pre'
#        
#        #x = []
#        for v in ratios[e]['pre']:
#            out_row_mass['prey'] = v[0]
#            out_row_mass['predator'] = v[1]
#            out_row_mass['value'] = v[2]
#            
#            out_row_mass['prey_group'] = masses[v[0]]['group']
#            out_row_mass['prey_tl'] = trophic_levels[e]['pre'][v[0]]
#            
#            out_row_mass['predator_group'] = masses[v[1]]['group']
#            out_row_mass['predator_tl'] = trophic_levels[e]['pre'][v[1]]
#            
#            
#            #x.append(v)
#            out_mass.writerow(out_row_mass)
#        
#        fig = plt.figure()
#        ax = fig.add_subplot(111)
#        fig.suptitle('Prey mass - elevation '+str(e))
#        
#        # the histogram of the data
#        a, bins, patches = ax.hist(y, 10, normed=1, facecolor='green', alpha=0.5, label=['pre'])
#
#        
#        #plt.hist(mass_dist)
##        plt.plot(mass_dist, in_dist_plot, 'o', label='pre'+str(e))
##        plt.show()
#        
#        out_row_in['period'] = 'post'
#        out_row_out['period'] = 'post'
#        out_row_st['period'] = 'post'
#        mass_dist = []
#        in_dist_plot = []
#        
#        y=[]
#        
#        for sp in in_degs[e]['post'].keys():
#            out_row_in['species'] = sp
#            out_row_out['species'] = sp
#            
#            ind = in_degs[e]['post'][sp]['original']
#            out_row_in['L'] = ind
#            out_row_in['L_standardised'] = in_degs[e]['post'][sp]['standardised']
#            
#            out_row_in['TP'] = trophic_positions[e]['post'][sp]
#            out_row_in['TL'] = trophic_levels[e]['post'][sp]
#            
#            out_row_in['group'] = masses[sp]['group']
#            
#            if math.ceil(ind) > 0:
#                out_row_in['mass'] = masses[sp]['biomass']
#                biomass = float(masses[sp]['biomass'])
#                if biomass > 0.0:
#                    mass_dist.append(math.log10(biomass))
#                    in_dist_plot.append(ind)
#            else:
#                out_row_in['mass'] = ''
#                
#            outd = out_degs[e]['post'][sp]['original']
#            out_row_out['L'] = outd
#            out_row_out['L_standardised'] = out_degs[e]['post'][sp]['standardised']
#            
#            if outd > 0.0:
#                out_row_in['prey_mass'] = masses[sp]['biomass']
#                mass = float(masses[sp]['biomass']) 
#                if  mass > 0.0:
#                    y.append(math.log10(mass))
#            else:
#                out_row_in['prey_mass'] = ''
#            
#            
#            out_row_in['n_habitats'] = len(sp_habitats[sp])
#            
#            out_row_in['coming_from_0'] = '0'
#            out_row_in['coming_from_1'] = '0'
#            out_row_in['coming_from_2'] = '0'
#            out_row_in['coming_from_3'] = '0'
#            new = True
#            for e_prev in elevations:
#                if sp in sps_in_pre[e_prev]:
#                    out_row_in['coming_from_'+str(e_prev)] = '1'
#                    new = False
#            
#            if new:
#                out_row_in['coming_from_3'] = '1'
#            
#            out.writerow(out_row_in)
#            out.writerow(out_row_out)
#            
#            out_row_newcomers['species'] = sp
#            out_row_newcomers['group'] = masses[sp]['group'] 
#            out_row_newcomers['degree'] = out_degs[e]['post'][sp]['original'] + in_degs[e]['post'][sp]['original']
#            out_row_newcomers['TP'] = trophic_positions[e]['post'][sp]
#            out_row_newcomers['TL'] = trophic_levels[e]['post'][sp]
#            out_row_newcomers['mass'] = masses[sp]['biomass']
#            out_row_newcomers['n_habitats'] = len(sp_habitats[sp])
#            
#            out_newcomers.writerow(out_row_newcomers)
#            
#        out_row_st['mean'] = scipy.mean(mass_dist)
#        out_row_st['variance'] = scipy.var(mass_dist)
#        out_row_st['skewness'] = scipy.stats.stats.skew(mass_dist)
#        out_row_st['kurtosis'] = scipy.stats.stats.kurtosis(mass_dist)
#        
#        out_stats.writerow(out_row_st)
#        
#        out_row_mass['period'] = 'post'
#        #x = []
#        for v in ratios[e]['post']:
#            out_row_mass['prey'] = v[0]
#            out_row_mass['predator'] = v[1]
#            out_row_mass['value'] = v[2]
#            
#            out_row_mass['prey_group'] = masses[v[0]]['group']
#            
#            if trophic_levels[e]['post'].has_key(v[0]):
#                out_row_mass['prey_tl'] = trophic_levels[e]['post'][v[0]]
#            else:
#                out_row_mass['prey_tl'] = trophic_levels[e]['pre'][v[0]]
#            
#            out_row_mass['predator_group'] = masses[v[1]]['group']
#            
#            if trophic_levels[e]['post'].has_key(v[1]):
#                out_row_mass['predator_tl'] = trophic_levels[e]['post'][v[1]]
#            else:
#                out_row_mass['predator_tl'] = trophic_levels[e]['pre'][v[1]]    
#            
#            #x.append(v)
#            out_mass.writerow(out_row_mass)
#        
#        # the histogram of the data
#        b, bins, patches = ax.hist(y, 10, normed=1, facecolor='red', alpha=0.5, label=['post'])
#        
#        #ax.legend((a,b), ('pre','post'))
#        plt.legend()
#         
##        plt.plot(mass_dist, in_dist_plot, 'ro', label='post'+str(e))
##        plt.show()
#    
#    out_file.close()
#    out_file_stats.close()
#    out_file_mass.close()
#    out_file_newcomers.close()
#    
#    plt.show()
#    
#    return sps_in_post
        
def calculate_migrant_species(nets_dict):
    headers = ['elev', 'type']
    
    min_year = 1980
    years = sorted(nets_dict.keys())
    current_years = []
    for y in years:
        if y >= min_year:
            headers.append(str(y))
            current_years.append(y)
    
    max_year = max(current_years)
    elevations = nets_dict[min_year][1].keys()
    
    species_presence = dict()
    years_presence = dict()
    
    for y in current_years:
        if not years_presence.has_key(y):
            years_presence[y] = dict()
            
            
        for e in elevations:
            if not years_presence[y].has_key(e):
                years_presence[y][e] = dict()
                years_presence[y][e]['enter'] = []
                years_presence[y][e]['leave'] = []
            
            net = nets_dict[y][1][e]['all']
            
            for n in net.nodes():
                if not species_presence.has_key(n):
                    species_presence[n] = dict()
                
                if not species_presence[n].has_key(e):
                    species_presence[n][e] = []
                    
                species_presence[n][e].append(y)
                
    
    n_sp_year = dict()
    
    for y in current_years:
        
        unique_species = set()
        
        for e in elevations:
            species = nets_dict[y][1][e]['all'].nodes()
            
            
            for sp in species:
                if min(species_presence[sp][e]) == y:
                    years_presence[y][e]['enter'].append(sp)
                
                if max(species_presence[sp][e]) == y:
                    years_presence[y][e]['leave'].append(sp)
                    
                    
                unique_species.add(sp)
                
        n_sp_year[y] = len(unique_species)
    
    
    
    out_file = open('../output/migration_matrix.csv', 'w')
    out = csv.DictWriter(out_file, headers, delimiter=',')

    #out.writeheader()
    
    headers_dict = dict()
    for n in headers:
        headers_dict[n] = n
        
    out.writerow(headers_dict)
    
    out_row_enter = dict()
    out_row_leave = dict()
    
    out_row_enter['type'] = 'enter'
    out_row_leave['type'] = 'leave'
        
    for e in elevations:
        out_row_enter['elev'] = e
        out_row_leave['elev'] = e
        
        for y in current_years:
            string_enter = ''
            string_leave = ''
            
            for s in sorted(years_presence[y][e]['enter']):
                string_enter += s+'\n'
            
            for s in sorted(years_presence[y][e]['leave']):
                string_leave += s+'\n'
                
            out_row_enter[str(y)] = string_enter
            out_row_leave[str(y)] = string_leave 
            
        out.writerow(out_row_enter)
        out.writerow(out_row_leave)
            
    out_file.close()
    
    headers = []
    for y in n_sp_year.keys():
        headers.append(str(y))
    
    out_file = open('../output/species_year.csv', 'w')
    out = csv.DictWriter(out_file, headers, delimiter=',')

    headers_dict = dict()
    for n in headers:
        headers_dict[n] = n
        
    out.writerow(headers_dict)
    
    #out.writeheader()    
    out_row = dict()
    
    for y in n_sp_year.keys():
        out_row[str(y)] = n_sp_year[y]
    
    out.writerow(out_row)
    
    out_file.close()
    
    
def calculate_jaccard_matrix(nets_dict):
    """
    Calculates the Jaccard similarity index for each pair of networks in nets_dict and creates an output file:
    'jaccard_similarity.csv' with a matrix reporting the indexes calculated.
    """
    #we write the headers
    headers_dict = dict()
    headers_dict['name'] = 'name'
    
    networks = []
    
    years = sorted(nets_dict.keys())
    for y in years:
        habitats = sorted(nets_dict[y].keys())
        for h in habitats:
            altitudes = sorted(nets_dict[y][h].keys())
            for a in altitudes:
                seasons = sorted(nets_dict[y][h][a].keys())
                for s in seasons:
                    net_name = str(y)+'-'+str(h)+'-'+str(a)+'-'+str(s)
                    headers_dict[net_name] = net_name
                    
                    networks.append(net_name)
    
    out = csv.DictWriter(open('../output/jaccard_similarity.csv', 'w'), ['name']+networks, delimiter=',')

    out.writerow(headers_dict)    
    out_row = dict()
    
    for n in networks:
        y, h, a, s = n.split('-')
        net1 = nets_dict[int(y)][int(h)][int(a)][s]
        nodes1 = set(net1.nodes())
        
        out_row['name'] = n
        for n2 in networks:
            y2, h2, a2, s2 = n2.split('-')
            net2 = nets_dict[int(y2)][int(h2)][int(a2)][s2]
            nodes2 = set(net2.nodes())
            
            union = nodes1 | nodes2
            intersection = nodes1 & nodes2
            
            if len(union) == 0:
                jaccard_index = 'N/A'
            else:
                jaccard_index = 1 - (float(len(intersection))/float(len(union)))
            
            out_row[n2] = jaccard_index
        
        out.writerow(out_row)
            


class ThreadInvadedNets(threading.Thread):
    def __init__(self, net, year, habitat, elevation, season):
        threading.Thread.__init__(self)
        self.net = net.copy()
        self.year = year
        self.habitat = habitat
        self.elevation = elevation
        self.season = season
    
    def run(self):
        self.generate_invaded_networks(self.net,self.year,self.habitat,self.elevation,self.season)

    
    def generate_invaded_networks(self, net, year, habitat, elevation, season):
        name = str(year)+'-'+str(habitat)+'-'+str(elevation)+'-'+str(season)
        
        header_names = ['net_name', 'invasion', 'failure_reason']
        general_file = open('../output_inv/'+name+'_general.csv', 'w')
        general_out = csv.DictWriter(general_file, header_names, delimiter=',')
        
        #we write the headers
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
        
        general_out.writerow(headers_dict)
        
        #output file for the invaded networks
        header_names = ['type', 'group', 'gen', 'vul', 'ext', 'invader', 'lost', 'S', 'L', 'L/S','C', 'T', 'B', 'I', 'Ca', 'Loop', 'NCycles', 'O', 'T-B', 'T-I', 'I-I', 'I-B', 'GenSD', 'VulSD', 'MxSim', 'MaxChainLength', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'MeanShortestPathLength', 'complexity', 'components', 'cc', 'cc_rnd', 'diameter', 'diameter_rnd', 'radius', 'radius_rnd', 'compartmentalisation', 'modularity_rgraph', 'modules_rgraph', 'modularity_rgraph_randoms', 'sd_modularity_rgraph_randoms', 'chi_square', 'p_value', 'mean_mass_ratio', 'min_biomass', 'max_biomass', 'mean_biomass', 'R_removed', 'removed_species', 'R_CA_10', 'R_CA_20', 'R_CA_30', 'R_CA_40', 'R_CA_50', 'R_CA_60', 'R_CD_10', 'R_CD_20', 'R_CD_30', 'R_CD_40', 'R_CD_50', 'R_CD_60', 'R_MA_10', 'R_MA_20', 'R_MA_30', 'R_MA_40', 'R_MA_50', 'R_MA_60', 'R_MD_10', 'R_MD_20', 'R_MD_30', 'R_MD_40', 'R_MD_50', 'R_MD_60', 'R_TPA_10', 'R_TPA_20', 'R_TPA_30', 'R_TPA_40', 'R_TPA_50', 'R_TPA_60', 'R_TPD_10', 'R_TPD_20', 'R_TPD_30', 'R_TPD_40', 'R_TPD_50', 'R_TPD_60', 'R_R_10', 'R_R_20', 'R_R_30', 'R_R_40', 'R_R_50', 'R_R_60', 'Rdyn_CA_10', 'Rdyn_CA_20', 'Rdyn_CA_30', 'Rdyn_CA_40', 'Rdyn_CA_50', 'Rdyn_CA_60', 'Rdyn_CD_10', 'Rdyn_CD_20', 'Rdyn_CD_30', 'Rdyn_CD_40', 'Rdyn_CD_50', 'Rdyn_CD_60', 'Rdyn_MA_10', 'Rdyn_MA_20', 'Rdyn_MA_30', 'Rdyn_MA_40', 'Rdyn_MA_50', 'Rdyn_MA_60', 'Rdyn_MD_10', 'Rdyn_MD_20', 'Rdyn_MD_30', 'Rdyn_MD_40', 'Rdyn_MD_50', 'Rdyn_MD_60', 'Rdyn_TPA_10', 'Rdyn_TPA_20', 'Rdyn_TPA_30', 'Rdyn_TPA_40', 'Rdyn_TPA_50', 'Rdyn_TPA_60', 'Rdyn_TPD_10', 'Rdyn_TPD_20', 'Rdyn_TPD_30', 'Rdyn_TPD_40', 'Rdyn_TPD_50', 'Rdyn_TPD_60', 'Rdyn_R_10', 'Rdyn_R_20', 'Rdyn_R_30', 'Rdyn_R_40', 'Rdyn_R_50', 'Rdyn_R_60']
        invasive_file = open('../output_inv/'+name+'.csv', 'w')
        out_invasive = csv.DictWriter(invasive_file, header_names, delimiter=',')
        
        #we write the headers
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
        
        out_invasive.writerow(headers_dict)
    
        #for the 'per node' values
        header_names = ['net_name', 'node_name', 'species_mass(g)', 'indegree', 'outdegree', 'norm_indegree', 'norm_outdegree', 'betweeness', 'deg_centrality', 'cc', 'TL', 'TP', 'NumberOfPaths', 'O', 'module_rgraph', 'role_rgraph']
        per_node_file = open('../output_inv/'+name+'_per_node.csv', 'w')
        out_per_node = csv.DictWriter(per_node_file, header_names, delimiter=',')
        
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
            
        out_per_node.writerow(headers_dict)
        
        #for the 'per edge' values
        header_names = ['net_name', 'predator', 'predator_mass(g)', 'prey', 'prey_mass(g)', 'strength']
        per_edge_file = open('../output_inv/'+name+'_per_edge.csv', 'w')
        out_per_edge = csv.DictWriter(per_edge_file, header_names, delimiter=',')
        
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
            
        out_per_edge.writerow(headers_dict)
        
        #output file for the invaded networks
        header_names = ['net_name', 'group', 'gen', 'vul', 'ext', 'species_lost']
        lost_file = open('../output_inv/'+name+'_lost.csv', 'w')
        out_invasive_lost = csv.DictWriter(lost_file, header_names, delimiter=',')
        
        #we write the headers
        headers_dict = dict()
        for n in header_names:
            headers_dict[n] = n
        
        out_invasive_lost.writerow(headers_dict)
        
        out_row_general = dict()
        
        out_row_node = dict()
        out_row_edge = dict()
        out_row_invasive = dict()
        out_row_lost = dict()
        
        original_nodes = set(net.nodes())
        
        for gr in INTRODUCED_GROUPS:
            for gen in arange(0.1, 1.0, PARAMETERS_OFFSET):
                for vul in arange(0.1, 1.0, PARAMETERS_OFFSET):
                    for ext in arange(0.1, 1.0, PARAMETERS_OFFSET):
                        name_for_invasion = name+'-'+str(gr)+'-'+str(gen)+'-'+str(vul)+'-'+str(ext)
                        out_row_general['net_name'] = name_for_invasion
                        try:
                            net_pre, net_inv = net.get_invaded_network(introduced_group=gr, invasive_name=INTRODUCED_SPECIES, biomass_introduced=BIOMASS_INTRODUCED, generalism=gen, predators=vul, ext_threshold=ext)
                            out_row_general['invasion'] = True
                            out_row_general['failure_reason'] = ''
                            general_out.writerow(out_row_general)
                        
                        except NotInvadableNetwork as e:
                            out_row_general['invasion'] = False
                            out_row_general['failure_reason'] = e.value
                            general_out.writerow(out_row_general)
                            
                            net_pre = e.net
                            net_inv = None
                            
                        networks = dict.fromkeys(['pre', 'inv'], None)
                        networks['pre'] = net_pre
                        networks['inv'] = net_inv
                        
                        for k in networks.keys():
                            net_inv = networks[k]
                            
                            if net_inv == None:
                                continue
                        
                            net_inv.longest_path_length()
                            out_row_invasive['type'] = k
                            out_row_invasive['group'] = gr
                            out_row_invasive['gen'] = gen
                            out_row_invasive['vul'] = vul
                            out_row_invasive['ext'] = ext
                            out_row_invasive['invader'] = INTRODUCED_SPECIES
                            out_row_invasive['lost'] = ''
                            
                            if k == 'inv':
                                lost_species = original_nodes - set(net_inv.nodes())
                                out_row_invasive['lost'] = lost_species
                                
                                if len(lost_species) > 0:
                                    out_row_lost['net_name'] = name
                                    out_row_lost['group'] = gr
                                    out_row_lost['gen'] = gen
                                    out_row_lost['vul'] = vul
                                    out_row_lost['ext'] = ext
                                    for sp in lost_species:
                                        out_row_lost['species_lost'] = sp
                                        out_invasive_lost.writerow(out_row_lost)
                                
                            out_row_invasive['S'] = net_inv.number_of_nodes()
                            out_row_invasive['L'] = net_inv.size()
                            out_row_invasive['L/S'] = net_inv.linkage_density()
                            out_row_invasive['C'] = net_inv.connectance()
                            out_row_invasive['T'], top_sps = net_inv.top_predators()
                            out_row_invasive['B'], basal_sps = net_inv.basal(heterotrophs=False)
                            out_row_invasive['I'], inter_sps = net_inv.intermediate(heterotrophs=False)
                            
                            fractions = net_inv.get_links_fractions_between_levels(basal_sps, inter_sps, top_sps)
                            out_row_invasive['T-B'] = fractions['tb']
                            out_row_invasive['T-I'] = fractions['ti']
                            out_row_invasive['I-I'] = fractions['ii']
                            out_row_invasive['I-B'] = fractions['ib']
                            
                            gen_sd, vul_sd = net_inv.generality_vulnerability_sd()
                            out_row_invasive['GenSD'] = gen_sd
                            out_row_invasive['VulSD'] = vul_sd
                            
                            out_row_invasive['MxSim'] = net_inv.maximum_similarity()
                            
                            out_row_invasive['Ca'] = net_inv.cannibalism()
                            out_row_invasive['Loop'], out_row_invasive['NCycles'] = net_inv.fraction_in_loops()
                            out_row_invasive['MaxChainLength'] = net_inv.longest_path_length()
                            out_row_invasive['MeanShortestPathLength'] = net_inv.mean_path_length()
                            out_row_invasive['O'], omni_sps = net_inv.omnivory()
                            out_row_invasive['complexity'] = net_inv.complexity()
                            out_row_invasive['components'] = net_inv.components()
            
                            net_und_inv = net_inv.to_undirected()
                            cc = 0.0
                            cc_rnd = 0.0
                            diameter = 0
                            diameter_rnd = 0
                            radius = 0
                            radius_rnd = 0
                            if net_und_inv.number_of_nodes() > 0:
                                net_rnd = net_inv.generate_random_graph()
                                net_rnd_und = net_rnd.to_undirected()
                                
                                cc = nx.average_clustering(net_und_inv)
                                cc_rnd = nx.average_clustering(net_rnd_und)
                                
                                #we calculate structural measures of diameter and radius for the original network
                                try:
                                    eccs_net = nx.eccentricity(net_und_inv)
                                    diameter = max(eccs_net.values())
                                    radius = min(eccs_net.values())
                                except:
                                    pass
                                #we calculate structural measures of diameter and radius for the random network
                                try:
                                    eccs_net_rnd = nx.eccentricity(net_rnd_und)
                                    diameter_rnd = max(eccs_net_rnd.values())
                                    radius_rnd = min(eccs_net_rnd.values())
                                except:
                                    pass
                                    
                            out_row_invasive['cc'] = cc
                            out_row_invasive['cc_rnd'] = cc_rnd
                            out_row_invasive['diameter'] = diameter
                            out_row_invasive['diameter_rnd'] = diameter_rnd
                            out_row_invasive['radius'] = radius
                            out_row_invasive['radius_rnd'] = radius_rnd
                            
                            tps_inv, nops_inv, mean_length_inv = net_inv.find_trophic_positions()
                            out_row_invasive['MeanFoodChainLength'] = mean_length_inv
                            
                            length_var_inv, length_sd_inv, paths_no_inv = net_inv.get_path_length_feats()
                            out_row_invasive['ChnSD'] = length_sd_inv
                            out_row_invasive['ChnNo'] = paths_no_inv
                            
                            out_row_invasive['compartmentalisation'] = net_inv.degree_of_compartmentalization()
                            
                            if RGRAPH_MOD:
                                if net_inv.order() < 7:
                                    randoms = 0
                                else:
                                    randoms = RGRAPH_RANDOMIZATIONS 
            
                                modularity, no_modules = net_inv.modularity_rgraph(seed=RGRAPH_SEED, iter_factor=RGRAPH_ITERATION_F, cooling_factor=RGRAPH_COOLING_F, randoms=randoms)
                                out_row_invasive['modularity_rgraph'] = modularity
                                out_row_invasive['modules_rgraph'] = no_modules
                                out_row_invasive['modularity_rgraph_randoms'] = net_inv.modularity_of_randomizations
                                out_row_invasive['sd_modularity_rgraph_randoms'] = net_inv.sd_mod_of_randomizations
                            else:
                                out_row_invasive['modularity_rgraph'] = 'N/A'
                                out_row_invasive['modules_rgraph'] = 'N/A'
                                out_row_invasive['modularity_rgraph_randoms'] = 'N/A'
                                out_row_invasive['sd_modularity_rgraph_randoms'] = 'N/A'
                                
                            chi_sq_inv, p_val_inv = net_inv.chi_square_deg_freq()
                            out_row_invasive['chi_square'] = chi_sq_inv
                            out_row_invasive['p_value'] = p_val_inv
                            
                            out_row_invasive['mean_mass_ratio'] = net_inv.mean_body_mass_ratio()
                            
                            min_bm_inv, max_bm_inv, mean_bm_inv = net_inv.get_min_max_biomasses()
                            out_row_invasive['min_biomass'] = min_bm_inv
                            out_row_invasive['max_biomass'] = max_bm_inv
                            out_row_invasive['mean_biomass'] = mean_bm_inv
                            
                            if REMOVAL:
                                sps_to_remove = dr.get_species_to_remove()
                                if len(sps_to_remove) == 0:
                                    sps_to_remove = set(REMOVED_SPECIES)
                                if len(sps_to_remove) == 0:
                                    out_row_invasive['R_removed'] = 'N/A'
                                    out_row_invasive['removed_species'] = 'N/A'
                                else:
                                    n_removed, r_removed = net_inv.robustness_to_removal(to_remove=sps_to_remove, beta=REMOVAL_BETA, ext_threshold=REMOVAL_EXT_THRESHOLD)
                                    out_row_invasive['R_removed'] = r_removed
                                    out_row_invasive['removed_species'] = sorted(set(net_inv.nodes()) - set(n_removed.nodes()))
                            else:
                                out_row_invasive['R_removed'] = 'N/A'
                                out_row_invasive['removed_species'] = 'N/A'
        
        
                            #robustness calculations
                            if CALCULATE_ROBUSTNESS:
                                rob_conn_a = net_inv.robustness(criterion='conn', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                rob_conn_d = net_inv.robustness(criterion='conn', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                rob_mass_a = net_inv.robustness(criterion='mass', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                rob_mass_d = net_inv.robustness(criterion='mass', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                rob_tp_a = net_inv.robustness(criterion='trophic_position', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                rob_tp_d = net_inv.robustness(criterion='trophic_position', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                rob_random = net_inv.robustness(criterion='random', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, cumulative=ROBUSTNESS_CUMULATIVE)
                                
                                rob_conn_a_dyn = net_inv.robustness(criterion='conn', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                rob_conn_d_dyn = net_inv.robustness(criterion='conn', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                rob_mass_a_dyn = net_inv.robustness(criterion='mass', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                rob_mass_d_dyn = net_inv.robustness(criterion='mass', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                rob_tp_a_dyn = net_inv.robustness(criterion='trophic_position', ordering='asc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                rob_tp_d_dyn = net_inv.robustness(criterion='trophic_position', ordering='desc', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                rob_random_dyn = net_inv.robustness(criterion='random', interval=ROBUSTNESS_INTERVAL, max_removed=ROBUSTNESS_MAX_REMOVED, weighted=True, beta=ROBUSTNESS_BETA, ext_threshold=ROBUST_EXT_THRESHOLD)
                                
                                fractions = [10, 20, 30, 40, 50, 60]
                                for fr in fractions:
                                    #robustness by connectance
                                    #ascending
                                    value = 'N/A'
                                    key = 'R_CA_'+str(fr)
                                    if rob_conn_a.has_key(fr):
                                        value = rob_conn_a[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #descending
                                    value = 'N/A'
                                    key = 'R_CD_'+str(fr)
                                    if rob_conn_d.has_key(fr):
                                        value = rob_conn_d[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #robustness by biomass
                                    #ascending
                                    value = 'N/A'
                                    key = 'R_MA_'+str(fr)
                                    if rob_mass_a.has_key(fr):
                                        value = rob_mass_a[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #descending
                                    value = 'N/A'
                                    key = 'R_MD_'+str(fr)
                                    if rob_mass_d.has_key(fr):
                                        value = rob_mass_d[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #robustness by trophic position
                                    #ascending
                                    value = 'N/A'
                                    key = 'R_TPA_'+str(fr)
                                    if rob_tp_a.has_key(fr):
                                        value = rob_tp_a[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #descending
                                    value = 'N/A'
                                    key = 'R_TPD_'+str(fr)
                                    if rob_tp_d.has_key(fr):
                                        value = rob_tp_d[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #robustness random
                                    value = 'N/A'
                                    key = 'R_R_'+str(fr)
                                    if rob_random.has_key(fr):
                                        value = rob_random[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #dynamic robustness by connectance
                                    #ascending
                                    value = 'N/A'
                                    key = 'Rdyn_CA_'+str(fr)
                                    if rob_conn_a_dyn.has_key(fr):
                                        value = rob_conn_a_dyn[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #descending
                                    value = 'N/A'
                                    key = 'Rdyn_CD_'+str(fr)
                                    if rob_conn_d_dyn.has_key(fr):
                                        value = rob_conn_d_dyn[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #dynamic robustness by biomass
                                    #ascending
                                    value = 'N/A'
                                    key = 'Rdyn_MA_'+str(fr)
                                    if rob_mass_a_dyn.has_key(fr):
                                        value = rob_mass_a_dyn[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #descending
                                    value = 'N/A'
                                    key = 'Rdyn_MD_'+str(fr)
                                    if rob_mass_d_dyn.has_key(fr):
                                        value = rob_mass_d_dyn[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #dynamic robustness by trophic position
                                    #ascending
                                    value = 'N/A'
                                    key = 'Rdyn_TPA_'+str(fr)
                                    if rob_tp_a_dyn.has_key(fr):
                                        value = rob_tp_a_dyn[fr]
                                    
                                    out_row_invasive[key] = value
                                    
                                    #descending
                                    value = 'N/A'
                                    key = 'Rdyn_TPD_'+str(fr)
                                    if rob_tp_d_dyn.has_key(fr):
                                        value = rob_tp_d_dyn[fr]
                                    
                                    #dynamic robustness random
                                    out_row_invasive[key] = value
                                    value = 'N/A'
                                    key = 'Rdyn_R_'+str(fr)
                                    if rob_random_dyn.has_key(fr):
                                        value = rob_random_dyn[fr]
                                    
                                    out_row_invasive[key] = value
                            
                            out_invasive.writerow(out_row_invasive)
        
                            
                            #output per node
                            indegrees = net_inv.in_degree()
                            outdegrees = net_inv.out_degree()
                            total_in = sum(indegrees.values())
                            total_out = sum(outdegrees.values())
                            betweeness = nx.algorithms.centrality.betweenness_centrality(net_inv)
                            deg_centrality = nx.algorithms.centrality.degree_centrality(net_inv)
                            
                            #clustering coefficients
                            ccs = nx.clustering(net_und_inv)
                            
                            name_net = str(gr)+'-'+str(gen)+'-'+str(vul)+'-'+str(ext)
                            if k == 'pre':
                                name_net += '-pre'
                            
                            out_row_node['net_name'] = name_net
                            for vertex in net_inv.nodes():
                                out_row_node['node_name'] = vertex
                                
                                #if species_mass.has_key(vertex):
                                m = net_inv.node[vertex]['biomass']
                                if m == ' ':
                                    out_row_node['species_mass(g)'] = 'N/A'
                                else:
                                    out_row_node['species_mass(g)'] = m
                                #else:
                                #    out_row_node['species_mass(g)'] = 'N/A'
                                
                                indeg = indegrees[vertex]
                                out_row_node['indegree'] = indeg
                                
                                outdeg = outdegrees[vertex]
                                out_row_node['outdegree'] = outdeg
                                
                                out_row_node['norm_indegree'] = float(indeg)/float(total_in)
                                out_row_node['norm_outdegree'] = float(outdeg)/float(total_out)
                                out_row_node['betweeness'] = betweeness[vertex]
                                out_row_node['deg_centrality'] = deg_centrality[vertex]
                                out_row_node['cc'] = ccs[vertex]
                                
                                if vertex in basal_sps or indeg == 0:
                                    out_row_node['TL'] = 'B'
                                elif vertex in top_sps:
                                    out_row_node['TL'] = 'T'
                                elif vertex in inter_sps:
                                    out_row_node['TL'] = 'I'
                                else: 
                                    out_row_node['TL'] = 'N/A'
                                
                                if vertex in omni_sps:
                                    out_row_node['O'] = '1'
                                else:
                                    out_row_node['O'] = '0'
                                
                                
                                out_row_node['TP'] = tps_inv[vertex]
                                out_row_node['NumberOfPaths'] = nops_inv[vertex]
                                
                                if RGRAPH_MOD and net_inv.node[vertex].has_key('module'):
                                    out_row_node['module_rgraph'] = net_inv.node[vertex]['module']
                                    out_row_node['role_rgraph'] = net_inv.node[vertex]['role']
                                else:
                                    out_row_node['module_rgraph'] = 'N/A'
                                    out_row_node['role_rgraph'] = 'N/A'
                                
                                out_per_node.writerow(out_row_node)
                            
                            net_inv.obtain_interactions_strengths()
                            
                            out_row_edge['net_name'] = name_net 
                            for prey, predator, atts in net_inv.edges(data=True):
                                out_row_edge['predator'] = predator
                                out_row_edge['predator_mass(g)'] = net_inv.node[predator]['biomass']
                                
                                out_row_edge['prey'] = prey
                                out_row_edge['prey_mass(g)'] = net_inv.node[prey]['biomass']
                                
                                out_row_edge['strength'] = atts['weight']
                                
                                out_per_edge.writerow(out_row_edge)
                        
        general_file.close()
        invasive_file.close()
        per_edge_file.close()
        per_node_file.close()
        lost_file.close()