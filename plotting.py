
"""
:mod:`~plotting` is a module that provides plotting functionalities for displaying networks
created using the :class:`~web.Network` class for ecological networks representations; and
for displaying 2D data plots of 1, 2 and 3 series.
"""

import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import axes3d

from scipy.interpolate import interp1d
from scipy.optimize import leastsq

import networkx as nx

from utils import plfit_lsq

from config import RGRAPH_SEED, RGRAPH_ITERATION_F, RGRAPH_COOLING_F

class NetworkPlotter():
    """
    The NetworkPlotter class implements a method for plotting a network of the type :class:`~web.Network` using the plotting
    functionalities provided by the NetworkX library (http://networkx.lanl.gov), which in turns utilizes the 
    matplotlib API for displaying graphs (http://matplotlib.sourceforge.net/). 
    
    It also implements the logic for a basic set of user interactions with the network, via the mouse and the keyboard,
    for edge highlighting and node picking for visualising trophic relationships and second level effects of these
    interactions.
    """
    def plot_network(self, n, highlight=None, node_diameter='degree'):
        """
        This is the main method of this class, where the algorithms for drawing the network are implemented, including
        the assignment of different colours to different kinds of nodes and also to different kinds of edges. The nodes
        are classified according to the species groups defined for the species divisions and which each node possesses
        as an attribute in the ecological network; nodes are thus drawn in different colours according to the group
        they belong to. The group-colour dictionary is as follows::
        
            'Bird':'#87CEFA', 'BirdPrey':'#0000CD', 'Reptile':'y', 'Mammal':'#DEB887', 'MammalCarn':'#A52A2A', 'Fish':'#A9A9A9', 'Amphibian':'#ADFF2F', 'Invertebrate':'w', 'Other':'k'
        
        Also the diameter of the nodes are dependent on one of three possible features: the degree of the node within the
        network (mainly a topological feature), the biomass of the species represented by the node (obtained from the
        attribute *biomass* of the network's vertices), or the clustering coefficient of each node (a network metric
        that measures the extent to which a given node belongs to a cluster within the network). 
        
        Before the drawing process finishes, two kind of events over the produced figure are registered, one for the
        mouse and another for the keyboard actions, and the corresponding listeners for each are also assigned.
        
        call signature::
    
            plot_network(n, highlight=None)
        
        *n*
            an instance of the type :class:`Network` which is the network to be drawn 
        
        Optional keyword arguments:
        
        *highlight*
            a string specifying the id of a node in the network that will be highlighted (drawn in red) when
            the network is displayed. The default value for this argument is None, in which case no nodes are
            highlighted
            
        *node_diameter*
            this argument is a string specifying the criterion on which the diameter of the nodes in the network
            is going to be calculated. As mentioned above, for the current implementation there are three possible
            values for this argument: *'degree'* for the nodes' degree (this is the default value for this argument),
            *'biomass'* for the biomass of the species the node represents, or *'clustering'* a measure of the 
            nodes' clustering in the network. If the value given for this argument is different from all of these
            all the nodes in the network will be of the same size.
        """
        self.network = n.copy()
        self.fig = plt.figure()
        self.network_plot = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=-0.04, bottom=-0.01, top=1, right=1.03)
        self.network_plot.xaxis.set_ticks_position('none')
        self.network_plot.yaxis.set_ticks_position('none')
        self.network_plot.set_xticklabels([])
        self.network_plot.set_yticklabels([])
        
        
        tps = self.network.get_trophic_levels()
        
        print 'tps ', tps
        
        ranks = dict()
        for n in self.network.nodes():
            if not ranks.has_key(tps[n]):
                ranks[tps[n]] = []
                
            ranks[tps[n]].append(n)
        
        print 'ranks', ranks.values()
        
        print self.network.nodes(data=True)
        
#        self.network.remove_edges_from(self.network.selfloop_edges())
#        
#        removed_edges = []
#        for u,v in self.network.edges():
#            try:
#                self.network.remove_edge(u,v)
#                self.layout = nx.graphviz_layout(self.network, prog="dot", args='-Grankdir=BT -Granksep=.5 -Gnodesep=.01 -Nlabel=' ' -Nwidth=.05 -Nheight=.05', ranks=ranks.values())
#                removed_edges.append((u,v))
#                break
#        #self.layout = nx.circular_layout(self.network.to_undirected())
#            except:
#                removed_edges.append((u,v))
#                pass
#        
#        print 'removed_edges', removed_edges
#        self.network.add_edges_from(removed_edges)
        
        
        ##this bit implements a hand-made layout based on the ranking of species
        
        minx = 1.0
        maxx = 105.0
        miny = 1.0
        maxy = 85.0
        
        middlex = (maxx - minx)/2
        middley = (maxy - miny)/2
        
        levels = len(ranks.keys())
        level_size = (maxy - miny) / levels
        current_y = miny
        ys = []
        for i in range(levels):
            ys.append(current_y)
            current_y += level_size
            
        self.layout = dict.fromkeys(self.network.nodes())
        
        for r in ranks.keys():
            n_nodes = len(ranks[r])
            offset = (maxx - minx) / n_nodes
            #current_x = minx
            current_x = middlex - ((float(n_nodes-1)/2) *offset)
            for n in ranks[r]:
                self.layout[n] = (current_x, ys[r])
                current_x += offset
        
#        xs = dict()
#        minx = None
#        maxx = None
#        
#        for n in self.layout.keys():
#            if not xs.has_key(tps[n]):
#                xs[tps[n]] = []
#            
#            xs[tps[n]].append(self.layout[n][0])
#            
#            if minx == None or self.layout[n][0] < minx:
#                minx = self.layout[n][0]
#            
#            if maxx == None or self.layout[n][0] > maxx:
#                maxx = self.layout[n][0]
#        
#        middle = (maxx - minx)/2
#        
#        print 'minx, maxx ', minx, maxx
        
#        for n in xs.keys():
#            offset = (max(xs[n]) - min(xs[n])) / len(xs[n])
#            
#            #current_x = min(xs[n])
#            current_x = middle-((float(len(xs[n]))/2) *offset)
#            
#            if n == 0:
#                count = 0
#            for v in ranks[n]:
#                if n == 0 and count % 2 == 0:
#                    self.layout[v] = (current_x, self.layout[v][1])
#                    
#                else:
#                    self.layout[v] = (current_x, self.layout[v][1]+10.0)
#                
#                #self.layout[v] = (current_x, self.layout[v][1])
#                
#                count += 1
#                current_x += offset
                
                 
        
        self.color_reference = dict({'Bird':'#87CEFA', 'BirdPrey':'#0000CD', 'Reptile':'y', 'Mammal':'#DEB887', 'MammalCarn':'#A52A2A', 'Fish':'#A9A9A9', ' ':'k', 'N/A':'k', '':'k', 'Amphibian':'#ADFF2F', 'Invertebrate':'w'})
        #birds: blue family mammals:brown family, fish:grey, reptile:yellow, amphibian:green, and others ok
        
        #greyscale
        #self.color_reference = dict({'Bird':'#FFFFFF', 'BirdPrey':'#000000', 'MammalCarn':'#BABABA', 'Amphibian':'#878787', 'Reptile':'#545454', 'Fish':'#A9A9A9', ' ':'k', 'N/A':'k', '':'k', 'Mammal':'#DEDEDE', 'Invertebrate':'w'})

###### This piece of code is for painting nodes using grayscales...

        #grayscale = ['#FFFFFF', '#DEDEDE', '#BABABA', '#A9A9A9', '#878787', '#545454', '#000000', '#000000']
#        group_values = set()
#        for n,atts in self.network.nodes(data=True):
#            group_values.add(atts['group'])
#        
#        self.color_reference = dict.fromkeys(group_values)
#        current_colour = 0
#        n_groups = len(group_values)
#        
#        if n_groups == 1:
#            self.color_reference[0] = grayscale[current_colour]
#        else:
#            for g in group_values:
#                self.color_reference[g] = grayscale[current_colour]
#                current_colour += 1
#                if current_colour == n_groups:
#                    self.color_reference[g] = grayscale[len(grayscale)-1]
#            
        #self.color_reference = dict({0:'w', 1:'k'})
        
        post_nodes = []
        other_nodes = []
        for n, atts in self.network.nodes(data=True):
            try:
                if atts['in_post']:
                    post_nodes.append(n)
                else:
                    other_nodes.append(n)
            except:
                other_nodes = self.network.nodes()
                break
        
        if node_diameter == 'degree':
            sizs = nx.algorithms.centrality.degree_centrality(self.network)
            self.sizes = [1000*sizs[v] for v in self.network]
            
            sizes_others = [1000*sizs[v] for v in other_nodes]
            sizes_posts = [1000*sizs[v] for v in post_nodes]
            
        elif node_diameter == 'biomass':
            self.sizes = []
            
            sizes_others = []
            sizes_posts = []
            
            print 'nodes with no biomass: '
            
            for v in self.network:
                str_biomass = self.network.node[v]['biomass']
                biomass = 0.0
                if str_biomass == ' ':
                    biomass = 0.0
                else:
                    biomass = float(str_biomass)
                
                #biomass += 2.0
                
#                if biomass != 0.0:
#                    biomass = math.log10(biomass)
#                else:
#                    biomass = 2.0
                
                self.sizes.append(biomass*0.1)
                
                if self.network.node[v]['in_post']:
                    sizes_posts.append(biomass*60.0)
                else:
                    sizes_others.append(biomass*60.0)
                
        elif node_diameter == 'clustering':
            net_und = self.network.to_undirected()
            ccs = nx.clustering(net_und)
            self.sizes = [(300*ccs[v])+20 for v in self.network]
        else:
            print 'None of the possible values for nodes diameters specification was selected'
            self.sizes = 500
        
        self.colors = []
#        colors_posts = []
#        colors_others = []
        for v in self.network:
            #if highlighting by colour ('red') uncomment this
            if v == highlight or ('in_post' in self.network.node[v].keys() and self.network.node[v]['in_post']):
                self.colors.append('r')
                continue
           # else:
            #    self.colors.append(self.color_reference[self.network.node[v]['group']])
        
#            if self.network.node[v]['in_post']:
#                colors_posts.append(self.color_reference[self.network.node[v]['group']])
#            else:
#                colors_others.append(self.color_reference[self.network.node[v]['group']])
            
            if(tps[v] == 0):
                self.colors.append('#aac46b')
#                 if(self.network.node[v]['mut_prod'] == True):
#                     self.colors.append('#edea00')
#                 else:
#                     self.colors.append('#aac46b')
            
                
            if(tps[v] == 1):
                self.colors.append('#5f95c9')                    
#                 if(self.network.node[v]['mut'] == True):
#                     self.colors.append('#937ab1')
#                 else:
#                     self.colors.append('#5f95c9')

            
            if(tps[v] == 2):
                self.colors.append('#faa756')
            
            if(tps[v] == 3):
                self.colors.append('#cd665f')
            
                                                    
        
        
        self.nodes = nx.draw_networkx_nodes(self.network, self.layout, 
                                            ax=self.network_plot,
                                            nodelist=self.network.nodes(), 
                                            node_size=self.sizes,
                                            node_color=self.colors,
                                            node_shape='o', 
                                            with_labels=True)
        
#        other_nodes = nx.draw_networkx_nodes(self.network, self.layout, 
#                                            ax=self.network_plot,
#                                            nodelist=post_nodes, 
#                                            node_size=sizes_posts,
#                                            node_color=colors_posts,
#                                            node_shape='^', 
#                                            with_labels=False)
        
        self.nodes.set_picker(True)
        
        
#        self.network.obtain_interactions_strengths(normalise=True)
#        edge_color_set = ['c','b','m']
#        self.edges_colors = []
#        self.edge_styles = []
#        number_of_intervals = 3
#        
#        for n in self.network.nodes():
#            preys = self.network.predecessors(n)
#            max_weight = None
#            min_weight = None
#            if len(preys) == 0:
#                continue
#            for prey in preys:
#                w = self.network[prey][n]['weight']
#                if w == 0.0:
#                    continue
#                if min_weight == None or w < min_weight:
#                    min_weight = w
#                if max_weight == None or w > max_weight:
#                    max_weight = w
#        
##        for u,v,atts in self.network.edges(data=True):
##            if min_weight == None or atts['weight'] < min_weight:
##                min_weight = atts['weight']
##            if max_weight == None or atts['weight'] > max_weight:
##                max_weight = atts['weight']
#        
#        
#            offset = ((max_weight - min_weight)/number_of_intervals)
#        
#            weights_range = [min_weight]
#            next_value = min_weight
#            for i in xrange(number_of_intervals-1):
#                next_value += offset
#                weights_range.append(next_value)
#            
#            weights_range.append(max_weight)
#                
#            for prey in preys:
#                w = self.network[prey][n]['weight']
#                if w == 0.0:
#                    self.network[prey][n]['color'] = 'k'
#                else:
#                    for i in xrange(number_of_intervals):
#                        if w >= weights_range[i] and w <= weights_range[i+1]:
#                            self.network[prey][n]['color'] = edge_color_set[i]
                
        #this part of the algorithm finds the colours to be assigned to each edge based on the 
        #differences in biomass between the species that interact via that edge
        
#        for u,v,atts in self.network.edges(data=True):
#            if atts['weight'] == 0.0:
#                self.edge_styles.append('dotted')
#                self.edges_colors.append('k')
#                continue
#            else:
#                self.edge_styles.append('solid')
#            
#            prey_mass = self.network.node[u]['biomass']
#            predator_mass = self.network.node[v]['biomass']
#            
#            if prey_mass == None or predator_mass == None:
#                self.edges_colors.append('k')
#                continue
#            
#            prey_mass = float(prey_mass)
#            predator_mass = float(predator_mass)
#            if prey_mass > predator_mass*0.05:
#                self.edges_colors.append('c')
#            elif prey_mass > predator_mass*0.01:
#                self.edges_colors.append('b')
#            else:
#                self.edges_colors.append('m')
#            
        
        #self.edges = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, style=self.edge_styles, arrows=False, edge_color=self.edges_colors)
        
        self.edges = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, arrows=True, width=0.3, edge_color='grey')
        
        self.current_node = None
        self.old_node = None
        self.current_preds = []
        self.current_sucs = []
        self.edges_selected = None
        self.edges_second_level = None
        self.txt = self.network_plot.text(0, 0, " ")
        self.last_key = 0
        self.text_x = -5
        self.text_y = -5
        
        #here we declare the events that we want to capture on the canvas
        #and the corresponding listeners for each of the events
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        #plt.draw()
        #plt.show()
        #nx.draw_networkx(self.network, self.layout, ax=self.network_plot, node_size=[1000*sizes[v] for v in self.network], style='solid', with_labels=True, arrows=False)    

    
    def plot_network_no_classes(self, n):
        self.network = n.copy()
        self.fig = plt.figure()
        self.network_plot = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=-0.04, bottom=-0.01, top=1, right=1.03)
        self.network_plot.xaxis.set_ticks_position('none')
        self.network_plot.yaxis.set_ticks_position('none')
        self.network_plot.set_xticklabels([])
        self.network_plot.set_yticklabels([])
        
        self.layout = nx.graphviz_layout(self.network, prog="dot", args='-Gnodesep=1, -Granksep=.1, -Grankdir=BT')
        #layout = nx.spring_layout(self.network)
      
        sizs = nx.algorithms.centrality.degree_centrality(self.network)
        self.sizes = [1000*sizs[v] for v in self.network]
        
        
        self.nodes = nx.draw_networkx_nodes(self.network, self.layout, 
                                            ax=self.network_plot,
                                            with_labels=False)
        self.nodes.set_picker(True)
        
        self.edges = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, arrows=False)
        
        self.current_node = None
        self.old_node = None
        self.current_preds = []
        self.current_sucs = []
        self.edges_selected = None
        self.edges_second_level = None
        self.txt = self.network_plot.text(0, 0, " ")
        self.last_key = 0
        self.text_x = 30
        self.text_y = 10
        
        #here we declare the events that we want to capture on the canvas
        #and the corresponding listeners for each of the events
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.show()
        #nx.draw_networkx(self.network, self.layout, ax=self.network_plot, node_size=[1000*sizes[v] for v in self.network], style='solid', with_labels=True, arrows=False)    


    def plot_network_modules(self, n):        
        self.network = n.copy()
        mod, n_mods = self.network.modularity_rgraph(RGRAPH_SEED, RGRAPH_ITERATION_F, RGRAPH_COOLING_F, 0)
        self.fig = plt.figure()
        self.network_plot = self.fig.add_subplot(111)
        self.fig.subplots_adjust(left=-0.05, bottom=-0.01, top=1, right=1.01)
        self.network_plot.xaxis.set_ticks_position('none')
        self.network_plot.yaxis.set_ticks_position('none')
        self.network_plot.set_xticklabels([])
        self.network_plot.set_yticklabels([])
        self.fig.suptitle('Network modules')
#        
#        #we try to position important nodes of each module in a particular area of the plot
#        #to see how the spring layout works on that basis
#        divisions = math.ceil(n_mods/2)
#        offset = 1.0/(divisions+1)
#        modules_dict = dict.fromkeys(range(1,n_mods+1), None)
#        for n, atts in self.network.nodes(data=True):
#            current_module = atts['module'] 
#            if modules_dict[current_module] != None: 
#                if self.network.degree(n) > self.network.degree(modules_dict[current_module]):
#                    modules_dict[current_module] = n
#            else:
#                modules_dict[current_module] = n
#        
#        init_positions = dict()
#        y_pos = 0.25
#        second_offset = offset/2
#        current_div = 0
#        for k in sorted(modules_dict.keys()):
#            x_pos = offset*(current_div) + second_offset 
#            init_positions[modules_dict[k]] = [x_pos, y_pos]
#            if int(k % (divisions+1)) == 0:
#                y_pos = 0.75
#                current_div = 0
#            current_div += 1
#        
#        for n in self.network.nodes():
#            if not init_positions.has_key(n):
#                init_positions[n] = [np.random.random_sample(), np.random.random_sample()]
#        
#        fixed_nodes = init_positions.keys()
#        #self.layout = nx.graphviz_layout(self.network, prog='neato', args='-Gnodesep=5, -Granksep=.1, -Grankdir=BT')
#        self.layout = nx.spring_layout(self.network, pos=init_positions, fixed=fixed_nodes, iterations=1, scale=0.1)
        
        
        divisions = math.ceil(n_mods/2)
        offset = 1.0/(divisions+1)
        
        self.layout = dict()
        modules_dict = dict.fromkeys(range(1,n_mods+1))
        for n, atts in self.network.nodes(data=True):
            current_module = atts['module']
            if modules_dict[current_module] == None:
                modules_dict[current_module] = []
             
            modules_dict[current_module].append(n)
        
        min_y_pos = 0.0
        max_y_pos = 0.5 
        second_offset = offset/2
        current_div = 0
        for m in sorted(modules_dict.keys()):
            temp_net = self.network.subgraph(modules_dict[m])
            temp_layout = nx.spring_layout(temp_net)
            
            min_x_pos = offset*(current_div)
            max_x_pos = (offset*(current_div+1))-0.02
            
            min_x = max_x = min_y = max_y = None
            for n in temp_layout.keys():
                [x,y] = temp_layout[n]
                if min_x == None or x < min_x:
                    min_x = x
                if max_x == None or x > max_x:
                    max_x = x
                if min_y == None or y < min_y:
                    min_y = y
                if max_y == None or y > max_y:
                    max_y = y
            
            for n in temp_layout.keys():                
                [x,y] = temp_layout[n]
                x = float(x-min_x)/float(max_x-min_x)
                y = float(y-min_y)/float(max_y-min_y)
                
                x_pos = min_x_pos+(x*(max_x_pos-min_x_pos))
                y_pos = min_y_pos+(y*(max_y_pos-min_y_pos))
                
                temp_layout[n] = [x_pos,y_pos]
            
            current_div += 1
            if int(m % (divisions+1)) == 0:
                min_y_pos = 0.55
                max_y_pos = 1.0
                current_div = 0
            
            self.layout.update(temp_layout)
        
        self.color_reference = dict({1:'r', 2:'m', 3:'g', 4:'y', 5:'c', 6:'b', 7:'k', 8:'#ADFF2F', 9:'w'})

        self.colors = []
        self.node_labels = dict()
        for v in self.network:
            self.colors.append(self.color_reference[self.network.node[v]['module']])
            self.node_labels[v] = str(self.network.node[v]['role'])
        
#        sizs = nx.algorithms.centrality.degree_centrality(self.network)
#        self.sizes = [1000*sizs[v] for v in self.network]
        
        self.nodes = nx.draw_networkx_nodes(self.network, self.layout, 
                                            ax=self.network_plot,
                                            node_color=self.colors,
                                            node_size=150,
                                            with_labels=False)
        self.nodes.set_picker(True)
        
        self.edges = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, arrows=False)
        
        nx.draw_networkx_labels(self.network, self.layout, 
                                            ax=self.network_plot,
                                            labels=self.node_labels,
                                            font_size=10)
        
        self.current_node = None
        self.old_node = None
        self.current_preds = []
        self.current_sucs = []
        self.edges_selected = None
        self.edges_second_level = None
        self.txt = self.network_plot.text(0, 0, " ")
        self.last_key = 0
        self.text_x = 0.0
        self.text_y = 1.1
        
        #here we declare the events that we want to capture on the canvas
        #and the corresponding listeners for each of the events
        self.fig.canvas.mpl_connect('key_press_event', self.onkey)
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.show()

    def onpick(self, event):
        """
        The onpick method is called when the *pick_event* event occurs on the figure displaying the network. It is
        a private method that is used only inside the :meth:`plot_network` method when registering this event
        with the action listener. Once the pick event is captured on the figure and this method is called, some
        modifications are performed on the network graphic.
        
        In the current implementation only nodes' picking is supported. When a node is picked in the first instance
        all its relationships/links (with predators and prey alike) are highlighted in green, also the name
        of the picked species is shown in the lower left corner of the figure. If a node is picked in this way
        and then another node connected to it is also picked, the selection of edges does not change but the
        name of the second species is shown next to the first species name and the direction of the relationship
        is reported, in order to illustrate the interaction and its participants. 
        
        If a species that is currently connected to the selected one wants to be in turn picked in order to show
        its relationships, it has to be picked twice; so the links of the first species will be cleared and the 
        second species and its links will be shown.
        """
        ind=event.ind[0]
        redraw = self.current_node == self.network.nodes()[ind]
        self.current_node = self.network.nodes()[ind]
        
        self.last_key = 0
        
        if not (self.current_node in (self.current_preds + self.current_sucs)) or redraw:
            if self.edges_selected != None:
                try:
                    self.edges_selected.remove()
                except:
                    pass
            
            if self.edges_second_level != None:
                try:
                    self.edges_second_level.remove()
                except:
                    pass
            
            self.old_node = self.current_node
            edges_new = self.network.out_edges(self.current_node)+self.network.in_edges(self.current_node)
    
            self.current_sucs = self.network.successors(self.current_node)
            self.current_preds = self.network.predecessors(self.current_node)
    
            print edges_new
    
            self.edges_selected = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, edgelist=edges_new, width=4, style='solid', arrows=False, edge_color='g')


        self.txt.remove()
        if self.current_node in self.current_preds:
            self.txt = self.network_plot.text(self.text_x, self.text_y, self.old_node+' <- '+self.current_node)
        elif self.current_node in self.current_sucs:
            self.txt = self.network_plot.text(self.text_x, self.text_y, self.old_node+' -> '+self.current_node)
        else:
            self.txt = self.network_plot.text(self.text_x, self.text_y, self.old_node)

        self.fig.canvas.draw()

        print self.network.nodes()[ind]
        return True

    def onkey(self, event):
        """
        This method implements the action listener for the keyboard *key_press_event* event which currently captures the action performed
        over three keys on the keyboard:
        
            ===    ===============================================================================================================
            Key    Action
            ===    ===============================================================================================================
            2      shows the second level interactions of the selected node in red. When a node is picked on the network, all its
                   direct links are shown in green (see :meth:`onpick`). After the network has been modified in this way it is
                   possible to show the second level interactions (the interactions of the nodes that interact with the originally
                   picked one) using the key *2* on the keyboard
            3      operates as the key above (*2*) but does not take into account the secondary interactions that pass through the 
                   nodes 'ProdPri' and 'Invertebrates' which in our ecological networks are used to group the primary producers
                   and the invertebrates into two single nodes
            0      once any of the keys above has been pressed or any of the nodes in the network picked (showing in this way its
                   interactions) the network can be returned to its original representation (no highlighted links) by pressing the
                   key *0* on the keyboard  
            ===    ===============================================================================================================
        
        
        """
        if event.key == '2':
            if self.last_key != 2:
                try:
                    self.edges_second_level.remove()
                except:
                    pass
                edges_new = []
                if len(self.current_sucs) > 0:
                    edges_new = edges_new + self.network.out_edges(self.current_sucs)
                    for x in (set(self.network.in_edges(self.current_sucs)) - set(self.network.out_edges(self.current_node))):
                        edges_new.append(x) 
                
                if len(self.current_preds) > 0: 
                    edges_new = edges_new + self.network.in_edges(self.current_preds)
                    for x in (set(self.network.out_edges(self.current_preds)) - set(self.network.in_edges(self.current_node))):
                        edges_new.append(x) 
                
                if len(edges_new) > 0:
                    print edges_new
                    self.edges_second_level = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, edgelist=edges_new, width=2, style='solid', arrows=False, edge_color='r')
                    self.fig.canvas.draw()
                self.last_key = 2
            else:
                if self.edges_second_level != None and len(self.edges_second_level._paths) > 0:
                    self.edges_second_level.remove()
                    self.fig.canvas.draw()
                    self.last_key = 0
        
        if event.key == '3':
            if self.last_key != 3:
                try:
                    self.edges_second_level.remove()
                except:
                    pass
                edges_new = []
                exclude = set(['ProdPri','Invertebrates'])
                temp_sucs = set(self.current_sucs) - exclude
                if len(temp_sucs) > 0:
                    edges_new = edges_new + self.network.out_edges(temp_sucs)
                    for x in (set(self.network.in_edges(temp_sucs)) - set(self.network.out_edges(self.current_node))):
                        edges_new.append(x)
                
                temp_preds = set(self.current_preds) - exclude 
                if len(temp_preds) > 0:
                    edges_new = edges_new + self.network.in_edges(temp_preds)
                    for x in (set(self.network.out_edges(temp_preds)) - set(self.network.in_edges(self.current_node))):
                        edges_new.append(x) 
                
                if len(edges_new) > 0:
                    print edges_new
                    self.edges_second_level = nx.draw_networkx_edges(self.network, self.layout, ax=self.network_plot, edgelist=edges_new, width=2, style='solid', arrows=False, edge_color='r')
                    self.fig.canvas.draw()
                self.last_key = 3
            else:
                if self.edges_second_level != None and len(self.edges_second_level._paths) > 0:
                    self.edges_second_level.remove()
                    self.fig.canvas.draw()
                    self.last_key = 0
        
        if event.key == '0':
            try:
                self.txt.remove()
                self.edges_selected.remove()
                self.edges_second_level.remove()
            except:
                pass
            self.current_node = None
            self.old_node = None
            self.current_preds = []
            self.current_sucs = []
            self.edges_selected = None
            self.edges_second_level = None
            self.txt = plt.text(0, 0, " ")
            self.fig.canvas.draw()
            
        return
        
class DataStatsPlotter():
    """
    This class was implemented with the intention to provide an object that could encapsulate the functionalities of plotting
    several data series on 2-dimensional graphs in a generic way. It provides methods for plotting 3, 2 and 1 temporal data
    series on the same graph, in latter case also offering the option of fitting the data points with a power law curve
    calculated based on the least squares method.
    
    There is also a function that allows plotting the number of observations coming from the database per year, in which case
    the data received by the method needs to be provided specifically for years and observations reported on those years
    for species and number of individuals of each species.
    
    All the plots and graphics drawn through this class are instances of plotting objects provided by the matplotlib 
    library (http://matplotlib.sourceforge.net/)
    """
    #3d plot
    #    fig_conn = plt.figure()    
    #    conn = fig_conn.add_subplot(111, projection='3d')
    #    #conn.plot_wireframe(yrs_plot, nsp_plot, con_plot, rstride=10, cstride=10)
    #    conn.scatter(yrs_plot, nsp_plot, con_plot)
    #    conn.set_zlim3d(0,0.04, 'b-')

    def plot_series_points(self, years, x1, x2, y1, y2, x_label, y_label, title, style_lines=False, point_labels=False, legend_labels=['ss','aw']):
        """
        This method provides a generic way of plotting two 3D data time series where the third dimension is represented
        as labels on the plots and in our specific case corresponds to the years in which the other two parameters
        have been extracted. 
        
        call signature::
    
            plot_series_points(years, x1, x2, y1, y2, x_label, y_label, title, style_lines=False, point_labels=False, legend_labels=['ss','aw'])
        
        *years*
            an array with the years in which the data to be plotted where registered, this array can actually contain
            any series but it was originally thought for the years in our time series. The values in this array
            are represented as labels on the 2d plot.
        *x1*
            the x's values for the first data series (an array of numbers)
        *x2*
            the x's values for the second data series (an array of numbers)
        *y1*
            the y's values for the first data series (an array of numbers)
        *y2*
            the y's values for the second data series (an array of numbers)
        *x_label*
            a string: the label for the x axis of the graph
        *y_label*
            a string: the label for the y axis of the graph
        *title*
            a string: the title to be shown on the top of the plot
            
        Optional keyword arguments:
        
        *style_lines*
            a boolean that specifies whether to join the data points using a line, i.e. whether to represent the data 
            only as points or as points and lines. Its default value is False (i.e. only points are displayed on the graph)
            
        *point_labels*
            a boolean stating whether to display labels on the plots, if True, the labels are taken from the array of values
            given as the first argument of the function (*years*). It is set to False by default.
            
        *legend_labels*
            an array of two strings containing the labels of each of the time series displayed. The default value is:
            ['ss', 'aw'], as this method was originally conceived for plotting the data from the spring-summer (ss) and 
            autumn-winter (aw) networks.
        """
        fig = plt.figure()    
        plot = fig.add_subplot(111)
        plot.hold(True)
        
        if style_lines:
            style1 = 'ro-'
            style2 = 'bo-'
        else:
            style1 = 'ro'
            style2 = 'bo'
        
        plot.plot(x1, y1, style1)
        plot.plot(x2, y2, style2)
        plot.legend(legend_labels)
        
        if point_labels:
            for i, label in enumerate(years):
                plt.text(x1[i], y1[i], label)
                plt.text(x2[i], y2[i], label)
        
        plot.set_xlabel(x_label)
        plot.set_ylabel(y_label)
        plot.set_title(title)
        plot.hold(False)
    
    
    def plot_three_series_points(self, years, x1, x2, x3, y1, y2, y3, x_label, y_label, title, style_lines=False, point_labels=False, legend_labels=['elevation 0','elevation 1', 'elevation 2']):
        """
        This method provides a generic way of plotting three 3D data time series where the third dimension is represented
        as labels on the plots and in our specific case corresponds to the years in which the other two parameters
        have been extracted
        
        call signature::
    
            plot_three_series_points(self, years, x1, x2, x3, y1, y2, y3, x_label, y_label, title, style_lines=False, point_labels=False, legend_labels=['elevation 0','elevation 1', 'elevation 2'])
        
        *years*
            an array with the years in which the data to be plotted where registered, this array can actually contain
            any series but it was originally thought for the years in our time series. The values in this array
            are represented as labels on the 2d plot.
        *x1*
            the x's values for the first data series (an array of numbers)
        *x2*
            the x's values for the second data series (an array of numbers)
        *x3*
            the x's values for the third data series (an array of numbers)
        *y1*
            the y's values for the first data series (an array of numbers)
        *y2*
            the y's values for the second data series (an array of numbers)
        *y3*
            the y's values for the third data series (an array of numbers)
        *x_label*
            a string: the label for the x axis of the graph
        *y_label*
            a string: the label for the y axis of the graph
        *title*
            a string: the title to be shown on the top of the plot
            
        Optional keyword arguments:
        
        *style_lines*
            a boolean that specifies whether to join the data points using a line, i.e. whether to represent the data 
            only as points or as points and lines. Its default value is False (i.e. only points are displayed on the graph)
            
        *point_labels*
            a boolean stating whether to display labels on the plots, if True, the labels are taken from the array of values
            given as the first argument of the function (*years*). It is set to False by default.
            
        *legend_labels*
            an array of three strings containing the labels of each of the time series displayed. The default value is:
            ['elevation 0','elevation 1', 'elevation 2'], as this method was originally conceived for plotting the data from
            the three elevations normally considered in our data analyses.
        """
        fig = plt.figure()    
        plot = fig.add_subplot(111)
        plot.hold(True)
        
        if style_lines:
            style1 = 'ro-'
            style2 = 'bo-'
            style3 = 'go-'
        else:
            style1 = 'ro'
            style2 = 'bo'
            style3 = 'go'
        
        plot.plot(x1, y1, style1, alpha=0.5)
        plot.plot(x2, y2, style2, alpha=0.5)
        plot.plot(x3, y3, style3, alpha=0.5)
        plot.legend(legend_labels)
        
        if point_labels:
            for i, label in enumerate(years):
                plt.text(x1[i], y1[i], label)
                plt.text(x2[i], y2[i], label)
                plt.text(x3[i], y3[i], label)
        
        plot.set_xlabel(x_label)
        plot.set_ylabel(y_label)
        plot.set_title(title)
        plot.hold(False)
    
    def plot_observations_per_year(self, years, species_nos, obs_nos, obs_nos_ps):
        """
        This method provides a way of automating the plotting of three key values of species observations in our data: 
        species numbers, number of observations (individuals), and number of observed individuals per species. It then
        calls the method :meth:`plot_single_series` of this same class for plotting the first to data series and
        for the third one it algorithmically obtains different colours for each and displays the result in another
        plot. This last plot actually displays a time series for each species so it will show as many data series
        as species are reported in the argument obs_nos_ps.
        
        call signature::
    
            plot_observations_per_year(years, species_nos, obs_nos, obs_nos_ps)
        
        *years*
            an array with the years in which the data to be plotted where registered, this array can actually contain
            any series but it was originally thought for the years in our time series.
        *species_nos*
            an array containing the number of species seen for each year in *years*
        *obs_nos*
            an array of integers containing the number of observations registered for each year in *years*
        *obs_nos_ps*
            a dictionary containing N arrays of integers reporting the number of observations of each species
            for each year in *years*, where N is the number of species and the keys in the dictionary are the
            name of the reported species 
        """
        
        self.plot_single_series(years, species_nos, 'Years', 'S', 'Number of Species observed per year', 'b-')
        
        self.plot_single_series(years, obs_nos, 'Years', 'Number of observations', 'Number of observations per year', 'r-')
        
        colors = ['b','g','r','c','m','y','k'] #,'w'
        styles = ['.',',','o','v','^','<','>','1','2','3','4','s','p','*','h','H','+','x','D','d','|','_']
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.set_xlabel('Years')
        ax.set_ylabel('Number of observations')
        ax.set_title('Number of observations per species per year')
        
        qty_colors = len(colors)
        qty_styles = len(styles)
        names = []
        cant = 0
        for k2 in obs_nos_ps.keys():
            style = colors[cant%qty_colors] + styles[cant%qty_styles] + '-'
            ax.plot(years,obs_nos_ps[k2], style)
            names.append(k2)
            cant += 1
            if cant > 20:
                break

        leg = ax.legend(names)
        for t in leg.get_texts():
            t.set_fontsize('small')
        
    
    def plot_single_series(self, xdata, ydata, x_label, y_label, title, style, log_scale=False, fitting=None):
        """
        This method implements the plotting of 2D data series including some additional features for the
        presentation of the data and the configuration of the plot, such as log scaling of the axes and
        the points fitting using power-law curves. It displays the plots of the given series using the 
        matplotlib libraries.
        
        call signature::
    
            plot_single_series(xdata, ydata, x_label, y_label, title, style, log_scale=False, fitting=None)
        
        *xdata*
            the x's values of the data series
        *ydata*
            the y's values of the data series
        *x_label*
            a string: the label for the x axis of the graph
        *y_label*
            a string: the label for the y axis of the graph
        *title*
            a string: the title to be shown on the top of the plot
        *style*
            any of the possible styles provided by matplotlib for displaying the plotted data
            
        Optional keyword arguments:
        
        *log_scale*
            a boolean specifying whether to log scale the axis of the plot
        *fitting*
            a string that determines the kind of curve fitting that will be used to fit the data in the series
            as a consequence of the fitting a line will be drawn on the plot representing the best curve
            fitting for the type of fitting specified. Currently there are two implemented curve fittings:
            *'powerlaw'* (a power-law fitting) and *'truncated'* (a truncated power-law); both of these
            are implemented using a least squares method for curve approximation. Another possible value
            for this argument is 'all', in which case both fittings will be calculated and shown on the plot.
            The default value for this argument is None, in which case no fitting will be calculated and
            consequently no curves are going to appear on the plot.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        
        if fitting == 'powerlaw' or fitting == 'all':
            xdata = array(xdata)
            ydata = array(ydata)
            
            #yerr = 0.2 * ydata
            
            # Power-law fitting is best done by first converting
            # to a linear equation and then fitting to a straight line.
            #
            #  y = a * x^b
            #  log(y) = log(a) + b*log(x)
            
            # truncated power law: y = x^b * e(x/c)
            
            logx = log(xdata)
            logy = log(ydata)
            
            #logyerr = yerr / ydata
            
            # define our (line) fitting function
            fitfunc = lambda p, x: p[0] + p[1] * x
            
            # line fitting function for the truncated power law
            #fitfunc_trunc = lambda p, x, x2: p[0] + p[1] * x + ((x/p[2])/log(10))
            #fitf_trunc = lambda p, x, x2, y: (y - fitfunc_trunc(p, x, x2))
            
            fitf = lambda p, x, y: (y - fitfunc(p, x))
            
            #errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
            
            pinit = [1.0, -1.0]
            #out = leastsq(errfunc, pinit, args=(logx, logy, logyerr), full_output=1)
            
            out = leastsq(fitf, pinit, args=(logx, logy), full_output=1)
            
            pfinal = out[0]
            covar = out[1]
            #print pfinal
            #print covar
            
            #powerlaw = amplitude * xdata**exponent
            powerlaw = lambda x, amp, index: amp * (x**index)
            
            #truncated_powerlaw = lambda x, amp, index, index2: amp * (x**index) * exp((x/index2))
            
            #index2 = pfinal[2]
            index = pfinal[1]       #exponent
            #amp = 10.0**pfinal[0]   #amplitude
            amp = exp(pfinal[0])   #amplitude
            
            indexErr = sqrt( covar[0][0] )          #error in the exponent index
            ampErr = sqrt( covar[1][1] ) * amp      #error in the amplitude
            
            print index, indexErr
            print amp, ampErr
            
            print plfit_lsq(xdata,ydata)
            
            ax.plot(xdata, powerlaw(xdata, amp, index))
            
        if fitting == 'truncated' or fitting == 'all':
            xdata = array(xdata)
            ydata = array(ydata)
            # Power-law fitting is best done by first converting
            # to a linear equation and then fitting to a straight line.
            #
            
            # truncated power law: y = a * x^b * e(x/c)
            # log(y) = log(a) + b*log(x) + x/c
            
            logx = log(xdata)
            logy = log(ydata)
            
            # define our (line) fitting function
            fitfunc = lambda p, x, xraw: p[0] + p[1] * x + (xraw/p[2])
            fitf = lambda p, x, xraw, y: (y - fitfunc(p, x, xraw))
            
            pinit = [1.0, -1.0, 1.0]
            out = leastsq(fitf, pinit, args=(logx, xdata, logy), full_output=1)
            
            pfinal = out[0]
            covar = out[1]
            
            truncpl = lambda amp, x, index, subindex: amp * (x**index) * exp(x/subindex)
            
            amp = exp(pfinal[0])
            index = pfinal[1]       #exponent
            subindex = pfinal[2]
            
            indexErr = sqrt( covar[1][1] )          #error in the exponent index
            subErr = sqrt( covar[2][2] )     #error in the subindex
            
            print index, indexErr
            print subindex, subErr
            
            ax.plot(xdata, truncpl(amp, xdata, index, subindex))
            
        if log_scale:
            ax.set_xscale('log')
            ax.set_yscale('log')
        
        ax.set_title(title)
        ax.plot(xdata,ydata,style)