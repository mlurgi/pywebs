

from random import Random
import csv
import numpy
from scipy.stats import beta, uniform
import networkx as nx

from meliae import scanner

from guppy import hpy

from web import Network

class NetworkCreator():
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self.rnd = Random()
        self.rnd_uniform = uniform()
        self.net = Network()
    
#        self.invaders = []
    
    def reset_state(self):
        self.rnd = Random()
        self.rnd_uniform = uniform()
        self.net.clear()
        
    def create_niche_model_network(self, S, C, threshold):
        '''
        This is an implementation of the niche model
        '''
        ns = self.rnd_uniform.rvs(size=S)
        mean = (1/(2*C)) - 1
        self.beta_dist = beta(1, mean)
        
        lower_c = C - (C*threshold)
        upper_c = C + (C*threshold)
        
        while self.net.size() == 0 or not nx.is_connected(self.net.to_undirected()) or self.net.connectance() > upper_c or self.net.connectance() <lower_c:
            self.net.clear()

            basal = None
            smallest_n = None
            
            #here we obtain the fundamental niche values for each one of the species in the network
            # n, c, r
            for i in xrange(1,S+1):
                self.net.add_node(i)
                self.net.node[i]['n'] = float(ns[i-1]) #self.rnd_uniform.rvs()
                self.net.node[i]['r'] = self.beta_dist.rvs() * self.net.node[i]['n']
                self.net.node[i]['c'] = self.rnd.uniform((self.net.node[i]['r']/2), min(self.net.node[i]['n'], (1-(self.net.node[i]['r']/2))))
                
                #self.net.node[i]['c'] = self.rnd.uniform(self.net.node[i]['r']/2, self.net.node[i]['n'])
                
                #the original value of c (commented bit) as originally presented in the Nature 2000 paper
                #was changed after reading the niche model specification presented in the JAE 2008 paper
                #the new specification ensures that the r of all the species always lies within the niche interval [0,1]
                
                if smallest_n == None or self.net.node[i]['n'] < smallest_n:
                    smallest_n = self.net.node[i]['n']
                    basal = i
                          
            self.net.node[basal]['r'] = 0.0
            self._create_links_based_on_fundamental_niche()
            
        return self.net
            
    def _create_links_based_on_fundamental_niche(self, nodes_to_link=None):
        #based on the fundamental niche values obtained above we construct the network by adding
        #the corresponding links according to the species' niche and feeding range
        
        if nodes_to_link == None:
            nodes_to_link = set(self.net.nodes())
        
        #we iterate overall the nodes and assign links when necessary
        for i in self.net.nodes():
            if self.net.node[i]['r'] == 0.0:
                continue
            
            r_lower_bound = self.net.node[i]['c'] - (self.net.node[i]['r']/2)
            r_upper_bound = self.net.node[i]['c'] + (self.net.node[i]['r']/2) 
            for j in self.net.nodes():
                if i not in nodes_to_link and j not in nodes_to_link:
                    continue
                
                if self.net.node[j]['n'] >= r_lower_bound and self.net.node[j]['n'] <= r_upper_bound:  
                    self.net.add_edge(j,i)
        
        #for disconnected or duplicated nodes
        disc_nodes = set()        
        self_loops = set(self.net.selfloop_edges())
        
        for i in nodes_to_link:
            disconnected = False
            
            if self.net.degree(i) == 0:
                disconnected = True
            elif self.net.in_degree(i) == 1 and (i,i) in self_loops:
#                print 'producer with selfloop'
                if self.net.out_degree(i) > 1:
                    self.net.remove_edge(i,i)
                else:
                    disconnected = True
                    self.net.remove_node(i)
                    self.net.add_node(i)
            else:
                i_succs = set(self.net.successors(i))
                i_predecs = set(self.net.predecessors(i))
                for j in self.net.nodes():
                    j_succs = self.net.successors(j)
                    j_predecs = self.net.predecessors(j)
                    
                    if i_succs == j_succs and i_predecs == j_predecs:
                        disconnected = True
                        self.net.remove_node(i)
                        self.net.add_node(i)
                        print 'duplicated node'
                        break
                    
            #we reassign the fundamental niche values to nodes that are disconnected or duplicated        
            if disconnected:
                self.net.node[i]['n'] = self.rnd_uniform.rvs()
                self.net.node[i]['r'] = self.beta_dist.rvs() * self.net.node[i]['n']
                self.net.node[i]['c'] = self.rnd.uniform(self.net.node[i]['r']/2, self.net.node[i]['n'])
                disc_nodes.add(i)
        
        
        if len(disc_nodes) > 0:
            self._create_links_based_on_fundamental_niche(disc_nodes)
        
        return

if __name__ == '__main__':
    
    h = hpy()
    
    replicates = 1000
    c_thresh = 0.03
    nc = NetworkCreator()
    
    in_file = open('../input/study_networks_agg_0.75.csv', 'rb')
    nets_atts = csv.DictReader(in_file, delimiter=',')
    
    properties = ['S', 'L' , 'L/S', 'C', 'T', 'B', 'I', 'GenSD', 'VulSD', 'MxSim', 'MeanFoodChainLength', 'ChnSD', 'ChnNo', 'Ca', 'Loop', 'O']
    
    header_names = ['year', 'habitat', 'elevation', 'season']
    suffixes = ['(mean)', '(median)', '(error)']
    for p in properties:
        for s in suffixes:
            h_string = p+s
            header_names.append(h_string)
    
    out_file = open('../output/output_niche.csv', 'w')
    out = csv.DictWriter(out_file, header_names)
    
    out.writeheader()
    out_row = dict()
    
    for row in nets_atts:
        empirical_values = dict()
        model_values = dict()
        for p in properties:
            empirical_values[p] = float(row[p])
            model_values[p] = []
        
        for i in range(replicates):
            print row['year'], row['elevation'], i
            nc.reset_state()
            S = int(row['S'])
            C = float(row['C'])
            
            net = nc.create_niche_model_network(S, C, c_thresh)
            
            model_values['S'].append(net.number_of_nodes())
            model_values['L'].append(net.size())
            model_values['L/S'].append(net.linkage_density())
            
            conn = net.connectance()
            #print conn
            model_values['C'].append(conn)
            
            a,b = net.top_predators()
            model_values['T'].append(a)
            
            a,b = net.basal()
            model_values['B'].append(a)
            
            a,b = net.intermediate()
            model_values['I'].append(a) 
            
            gen_sd, vul_sd = net.generality_vulnerability_sd()
            model_values['GenSD'].append(gen_sd)
            model_values['VulSD'].append(vul_sd)
            model_values['MxSim'].append(net.maximum_similarity())
            model_values['Ca'].append(net.cannibalism())
            
            a,b = net.fraction_in_loops()
            model_values['Loop'].append(a)
            
            net.longest_path_length()
            a,b =  net.omnivory()
            model_values['O'].append(a) 
            
            tps, nops, mean_length = net.find_trophic_positions()
            model_values['MeanFoodChainLength'].append(mean_length)
            
            length_var, length_sd, paths_no = net.get_path_length_feats()
            model_values['ChnSD'].append(length_sd)
            model_values['ChnNo'].append(paths_no)
            
        out_row['year'] = row['year']
        out_row['habitat'] = row['habitat']
        out_row['elevation'] = row['elevation']
        out_row['season'] = row['season']
        
        print h.heap().more
        
        scanner.dump_all_objects( 'out_meliae' )
        
        for p in properties:
            array_prop = model_values[p]
            mean = numpy.mean(array_prop)
            median = numpy.median(array_prop)
            std_dev = numpy.std(array_prop)
            
            if std_dev == 0.0:
                calc_error = 0.0
            else:
                calc_error = (mean - empirical_values[p])/std_dev
            
            out_row[p+'(mean)'] = mean
            out_row[p+'(median)'] = median
            out_row[p+'(error)'] = calc_error
        
        out.writerow(out_row)
    
    in_file.close()
    out_file.close()
    
    
    
    
    