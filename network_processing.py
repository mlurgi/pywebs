
"""
The module :mod:`~network_processing` implements a series of functions for the processing
of input data containing information about species, their habitats, spatial location,
and trophic relationships with other species in a given ecosystem and obtains the food
web representation of those interconnected species, represented as graphs and grouped
according to different criteria such as: seasonality, elevation, and habitat clusters
"""

from config import START_YEAR, END_YEAR, ELEVATION_SLOTS, INTERVAL
from web import Network
from utils import DatabaseReader

class NetworkCreator():
    """
    The NetworkCreator class provides an object through which the user can process data given as arguments to its
    methods and that contain information about species and their relationships in an ecosystem. Through
    the processing of these data the methods of this class are able to obtain the network of interactions
    of the species that comply with certain climatic, spatial, and temporal criteria.
    
    The class also possesses several attributes, data structures in particular, that are used to store data
    referent to the obtained networks and that can be used to plot the changes of several attibutes of the data
    along different dimensions. These structures are dynamically filled with the information obtained
    from the networks during the processing of the data.
    
    *dr*
        a data reader object which must possess certain methods that are used throughout the class for obtaining
        data from a data source. This was implemented in this way in order to allow for the development of data
        plug-ins for different data sources; making in this way transparent for this class any changes made on the
        way the data for the network creation is stored. We have developed data sources plug-ins for file-based
        databases and for relational databases. The default value for this argument is None, in which case
        an instance of the class :class:`~database_reader.DatabaseReader` is created and used for the NetworkCreator operations.
        This is the current status because our latest version of the data lives now in a MySQL relational database,
        even though it was previously file-based.
                
    .. seealso::

        :class:`~database_reader.DatabaseReader`
            for a relational database data source plug-in
        
        :class:`~files_data_reader.FilesDataReader` 
            for a files-based data source plug-in
    """
    def __init__(self, dr=None):
        #in this dictionary we store the obtained networks
        self.nets = dict()
        
        #the following arrays are useful for plotting statistics about the data
        #array keeping the considered years
        self.yrs = []
        
        #number of species in ss (spring-summer) networks for the years in yrs
        self.nsp_ss = []
        #number of species in aw (autumn-winter) networks for the years in yrs
        self.nsp_aw = []
        
        #connectances in ss networks for the years in yrs
        self.cons_ss = []
        #connectances in aw networks for the years in yrs
        self.cons_aw = []
        
        #number of links in ss networks for the years in yrs
        self.nlinks_ss = []
        #number of links in aw networks for the years in yrs
        self.nlinks_aw = []
        
        #linkage density in ss networks for the years in yrs
        self.ld_ss = []
        #linkage density in aw networks for the years in yrs
        self.ld_aw = []
        
        self.obs_ss = []
        self.obs_aw = []
        
        if dr == None:
            self.dr = DatabaseReader()
        else:
            self.dr = dr
    
    #the following are the getter methods for the data structures declared above
    #this data is generated through the create_networks_from_data method, and therefore
    #in order to obtain the data this method must be called before the getters can yield useful info
    def get_years_data(self):
        """
        Returns the array of the years from which the data on the networks has been considered 
        """
        return self.yrs
    
    def get_nsp_data(self):
        """
        Returns two arrays containing the number of species seen each year in self.yrs for the
        spring-summer and autumn-winter networks respectively 
        """
        return self.nsp_ss, self.nsp_aw

    def get_cons_data(self):
        """
        Returns two arrays containing the connectance for the spring-summer and autumn-winter 
        networks respectively (for each year in self.yrs) 
        """
        return self.cons_ss, self.cons_aw
    
    def get_nlinks_data(self):
        """
        Returns two arrays containing the number of links in the spring-summer and autumn-winter 
        networks respectively (for each year in self.yrs) 
        """
        return self.nlinks_ss, self.nlinks_aw
    
    def get_ld_data(self):
        """
        Returns two arrays containing the mean number of links per species in the spring-summer and autumn-winter 
        networks respectively (for each year in self.yrs) 
        """
        return self.ld_ss, self.ld_aw
    
    def get_obs_data(self):
        """
        Returns two arrays containing the number of individual observations for each year in self.yrs for the
        spring-summer and autumn-winter networks respectively 
        """
        return self.obs_ss, self.obs_aw
    
    def get_networks(self):
        """
        Returns a dictionary containing all the obtained networks from the given species and their links according
        to four different parameters: year, season, habitat, and elevation  
        """
        return self.nets
    
    def get_species_period(self, observations, current_yr):
        """
        This method takes as arguments the dictionary *observations*, which has a series of years as its keys and
        keeps record of the species observed each of those years; and the *current_yr* which is an integer representing the
        year for which we want to extract the species observations data. It then aggregates the information of the species
        observed during the interval of years from *current_yr* to *current_yr* + INTERVAL (which is a constant defined
        in the module :mod:`~config`) in order to return another dictionary with all the species reported as observed during that
        period of time and the number of times they were seen; and an array with the years that are represented in that interval.
        
        call signature::
    
            get_species_period(observations, current_yr)
        
        *observations*
            a dictionary containing the number of observations for each species classified by year. The years are thus
            the keys for that dictionary 
        *current_yr*
            the year starting from which the observed species are going to be aggregated in order to obtain the set of
            them that were reported for that time interval 
            
        This method returns a dictionary containing the name of species as keys and the number of observations as values
        for the period of time specified by *current_yr* + INTERVAL. An set of all the years considered for the interval
        is also returned in a second return variable.
        """
        years = set([current_yr])
        
        print '\n', current_yr, ' - ', current_yr+INTERVAL, '\n'
        
        species_period = observations[current_yr]
        if INTERVAL > 1:
            yr_tmp = current_yr + 1
            finish = current_yr + INTERVAL 
            if finish > END_YEAR:
                finish = END_YEAR+1
            
            while yr_tmp < finish:
                while not observations.has_key(yr_tmp):
                    yr_tmp += 1
                    
                sps_temp = observations[yr_tmp]
                for k in sps_temp.keys():
                    if species_period.has_key(k):
                        species_period[k] += sps_temp[k]
                    else:
                        species_period[k] = sps_temp[k]
                
                #print 'current year = ', current_yr, 'now visiting: ', yr_tmp
                years.add(yr_tmp)
                yr_tmp += 1 
            
        return species_period, years
    
    def create_nets_for_elevations(self, species_temp, nets_dict):
        """
        Given a set of elevation intervals, which are obtained in advance and kept as an attribute of the class, this method
        obtains the set of species that are available in those elevations according to the habitats in which they live
        and the elevations reported for those habitats in the spatial area considered. After obtaining the species that
        are present in a certain altitude range it calls the method :meth:`create_nets_from_seasonality` to obtain the 
        three networks that can be calculated for the different seasonality values and which contain the species obtained
        beforehand. At the end these obtained networks are stored in a dictionary (received as an argument) classified
        according to different seasonalities. This is done for each elevation, creating in this way a 2 level dictionary.
        
        call signature::
    
            create_nets_for_elevations(species_temp, nets_dict)
        
        *species_temp*
            a set containing the species that are going to be considered and classified according to altitudes for the 
            creation of the corresponding networks
        *nets_dict*
            a reference to the dictionary where the networks obtained for each seasonality and for each elevation are going
            to be stored after the calculations
        
        .. seealso::
    
            :meth:`create_nets_from_seasonality`
                for the way in which networks are obtained based on the seasonalities of the composing species
        """
        for k in self.elevation_intervals.keys():
            nets_current_alt = dict()
            current_min = self.elevation_intervals[k][0]
            current_max = self.elevation_intervals[k][1]
            
            print current_min, current_max
            species = set()
            
            for spec in species_temp:
                if not self.habitats.has_key(spec):
                    print 'Species '+spec+' does not have habitat'
                    continue
                sp_habs = self.habitats[spec]
                for hb in sp_habs:
                    if self.habs_elevs.has_key(hb):
                        current_half = current_min + ((current_max-current_min)/2)
                        #here we determine whether the habitat of the
                        hb_sp_min = self.habs_elevs[hb][0]
                        hb_sp_max = self.habs_elevs[hb][1]  
                        if ( hb_sp_min > current_min and hb_sp_max < current_max) or (hb_sp_min < current_half and hb_sp_min > current_min) or (hb_sp_max < current_max and hb_sp_max > current_half):
                            species.add(spec)
                            break 
                    
            print species 
        
            net_ss, net_aw, net_ay, net_all = self.create_nets_from_seasonality(species)
        
            nets_current_alt['ss'] = net_ss
            nets_current_alt['aw'] = net_aw
            nets_current_alt['ay'] = net_ay
            nets_current_alt['all'] = net_all
        
            nets_dict[k] = nets_current_alt
        
    def create_nets_from_seasonality(self, species):
        """
        By looking at the interactions and seasonalities between the set of species given as argument in *'species'*
        (the interactions and seasonalities are retrieved in advance and are stored as attributes of this class),
        this method creates instances of the object :class:`~web.Network` in order to build network representations
        of the trophic interactions between those species also related to the seasonalities in which the species are
        normally seen. Three networks are thus created: one for the spring-summer ('ss'), another for the autumn-winter
        ('aw') and the last one for the all-year ('ay') periods. Additionally, this method also assigns two attributes
        to each node/species in the networks: 'biomass' and 'group', which are useful species attributes that are
        used for different network features. 
        
        call signature::
    
            create_nets_from_seasonality(species)
        
        *species*
            the set of species that are going to be considered (through their trophic interactions and seasonalities)
            for the creation of the ecological networks produced by this method
        
        It returns three instances of the object :class:`~web.Network`: *net_ss*, *net_aw*, and *net_ay*; corresponding
        to the spring-summer, autumn-winter, and all-year networks respectively
        """
        net_ss = Network()
        net_aw = Network()
        net_ay = Network()
        net_all = Network()
        
       # leave_out = ['young', 'nests', 'eggs']
        
        for sp in species:
#             if any(suff in sp for suff in leave_out):
#                 print sp
#                 continue
#             ##for each species we look at its interactions
            if self.interactions.has_key(sp):
                preys_inter = self.interactions[sp]
                ##we iterate over the interactions to consider each separately
                for kp in preys_inter.keys():
                    food = preys_inter[kp]
                    for f in food:
                        
#                         if any(suff in f for suff in leave_out):
#                             continue
                        ##depending on what kind of interaction they are involved the appropriate node is added to the network
                        
                        ##we also confirm the species and its food co-occurr in the same season
                        sea_sp = self.seasons[sp]
                        
                        if self.seasons.has_key(f):
                            sea_f = self.seasons[f]
                        else:
                            if f in species or kp == 3:
                                print 'Node '+ f + ' does not have seasonality' 
                            sea_f = 'yr'
                        
                        if sea_sp == 'ay':
                            sea_sp = 'yr'
                        if sea_f == 'ay':
                            sea_f = 'yr'
                        
                        if (kp == 1 and (f in species)):
                            if sea_sp == 'yr':
                                if sea_f == 'ss':
                                    net_ss.add_edge(f,sp)
                                elif sea_f == 'aw':
                                    net_aw.add_edge(f,sp)
                                else:
                                    net_ss.add_edge(f,sp)
                                    net_aw.add_edge(f,sp)
                            elif sea_sp == 'ss':
                                if sea_f == 'ss' or sea_f == 'yr':
                                    net_ss.add_edge(f,sp)
                            elif sea_sp == 'aw':
                                if sea_f == 'aw' or sea_f == 'yr':
                                    net_aw.add_edge(f,sp)
                            #always add the edge to the control network
                            net_ay.add_edge(f,sp)
                            
                        elif kp == 3:
                            if sea_sp == 'aw' or sea_sp == 'yr':
                                net_aw.add_edge(f,sp)
                                net_ay.add_edge(f,sp)
                                if self.interactions.has_key(f):
                                    for ff in self.interactions[f][3]:
                                        net_aw.add_edge(ff, f)
                                        net_ay.add_edge(ff, f)
                            
                            if sea_sp == 'ss' or sea_sp == 'yr':
                                net_ss.add_edge(f,sp)
                                net_ay.add_edge(f,sp)
                                if self.interactions.has_key(f):
                                    for ff in self.interactions[f][3]:
                                        net_ss.add_edge(ff, f)
                                        net_ay.add_edge(ff, f)
        
        for v in net_ss:
            try:
                net_ss.node[v]['biomass'] = self.biomasses[v]['biomass']
            except:
                print 'Node '+ v + ' does not have biomass'
                net_ss.node[v]['biomass'] = ' '
            try:
                net_ss.node[v]['group'] = self.biomasses[v]['group']
            except:
                print 'Node '+ v + ' does not have group'
                net_ss.node[v]['group'] = 'N/A'
        for v in net_aw:
            try:
                net_aw.node[v]['biomass'] = self.biomasses[v]['biomass']
            except:
                print 'Node '+ v + ' does not have biomass'
                net_aw.node[v]['biomass'] = ' '
            try:
                net_aw.node[v]['group'] = self.biomasses[v]['group']
            except:
                print 'Node '+ v + ' does not have group'
                net_aw.node[v]['group'] = 'N/A'
        for v in net_ay:
            try:
                net_ay.node[v]['biomass'] = self.biomasses[v]['biomass']
            except:
                print 'Node '+ v + ' does not have biomass'
                net_ay.node[v]['biomass'] = ' '
            try:
                net_ay.node[v]['group'] = self.biomasses[v]['group']
            except:
                print 'Node '+ v + ' does not have group'
                net_ay.node[v]['group'] = 'N/A'
        
        net_all.add_nodes_from(net_ss.nodes(data=True))
        net_all.add_edges_from(net_ss.edges(data=True))
        
        net_all.add_nodes_from(net_aw.nodes(data=True))
        net_all.add_edges_from(net_aw.edges(data=True))
        
        net_all.add_nodes_from(net_ay.nodes(data=True))
        net_all.add_edges_from(net_ay.edges(data=True)) 
        
        return net_ss, net_aw, net_ay, net_all
    
    def create_networks_from_data(self, observations):
        """
        This method creates the networks representations of the interactions between species given as input through the argument
        *observations*. For each year between START_YEAR and END_YEAR (constants defined in the module :mod:`config`), it calls 
        the method :meth:`get_species_period` for obtaining the species that were observed in the current period as defined
        by the INTERVAL constant; after obtaining the set of species in this way, it calls on the method :meth:`create_nets_from_seasonality`
        in order to obtain the networks for each season (spring-summer = 'ss', autumn-winter = 'aw', and all-year = 'ay').
        It stores the obtained networks in the networks dictionary *self.nets* according to year and seasonality.
        Additionally, this method fills the arrays of number of species, observations and the network attributes such as 
        linkage density and connectance, that can be used for plotting the relationships of these parameters against
        time, and which can be obtained through the getters of this class.
        At the beginning of this method the information about species seasonalities, biomass and interactions are
        retrieved from the data reader and kept in data structures owned by this class.
        
        call signature::
    
            create_networks_from_data(observations)
        
        *observations*
            a dictionary with the observations of each species classified by years (the keys of the dictionary)
        
         .. seealso::
    
            :meth:`get_species_period`
                for the way in which the dictionary *observations* is processed
        """
        #for every year we obtain all the species seen on that year and the number of times they were seen
#        observations = get_species_observations()
        self.interactions, predators, preys = self.dr.get_interactions()
        self.seasons = self.dr.read_seasonality()
        self.biomasses = self.dr.read_species_mass()
       
        self.nets.clear()
        
        ##here we construct the interaction networks from the species/year records and the possible interactions
        for current_yr in xrange(START_YEAR, END_YEAR+1):
            if(INTERVAL < 1):
                break
            
            while not observations.has_key(current_yr):
                current_yr += 1
            
            species_period, years = self.get_species_period(observations, current_yr)
            
            species = species_period.keys()
            
            net_ss, net_aw, net_ay, net_all = self.create_nets_from_seasonality(species)
            
            self.nets[current_yr] = dict()
            self.nets[current_yr]['ss'] = net_ss
            self.nets[current_yr]['aw'] = net_aw
            self.nets[current_yr]['ay'] = net_ay
            self.nets[current_yr]['all'] = net_all
            
            ob_ss = 0
            for v in net_ss.nodes():
                if species_period.has_key(v):
                    ob_ss += species_period[v]
            
            ob_aw = 0
            for v in net_aw.nodes():
                if species_period.has_key(v):
                    ob_aw += species_period[v]
            
            #data for the plots
            self.yrs.append(current_yr)
            self.nsp_ss.append(net_ss.number_of_nodes())
            self.cons_ss.append(net_ss.connectance())
            self.nlinks_ss.append(net_ss.size())
            self.ld_ss.append(net_ss.linkage_density())
            self.obs_ss.append(ob_ss)
            
            self.nsp_aw.append(net_aw.number_of_nodes())
            self.cons_aw.append(net_aw.connectance())
            self.nlinks_aw.append(net_aw.size())
            self.ld_aw.append(net_aw.linkage_density())
            self.obs_aw.append(ob_aw)
    
    
    
    def create_networks_from_data_and_habitats(self, observations, clustering=False, high_res=True):
        """
        This method creates networks representations of the species in the *observations* data set in more or less the same way
        as detailed in :meth:`create_networks_from_data` but it also considers clusters of regions within the territory
        of study to create networks for different areas within it (clusters), based on the elevations, habitats for classifying 
        the species in order to obtain different networks for different elevations within each cluster of the territory.
        These clusters are defined in advance and are retrieved by the DataReader of this class through the method
        :meth:`~database_reader.read_habitat_clusters` which divides the area of study in georeferenced sections that
        possess many habitats. This results in a combinatorial number of networks in each cluster, according to the 
        habitats present in each cluster; and by elevation, according to the habitats encountered at each elevation
        within each cluster of land.
        Once the composing species based on clusters (or not) are obtained it calls the method :meth:`create_nets_for_elevations`
        in order to obtain the networks for each of the considered elevations.
        At the beginning of this method the information about species seasonalities, biomass, interactions and
        the habitats and their corresponding elevations within each grid are retrieved from the data reader and kept in 
        data structures owned either by this method or the class.
        
        call signature::
    
            create_networks_from_data_and_habitats(observations, clustering=False, high_res=True)
        
        *observations*
            a dictionary with the observations of each species classified by years (the keys of the dictionary)
        
        Optional keyword arguments:
            *clustering*
                a boolean stating whether the cluster classification of the territory is going to be considered for
                the production of the networks, its default value is False, in which case the networks are constructed
                only based on elevations and seasonalities per year.
            *high_res*
                a boolean that determines whether to consider each cluster homogeneously (low resolution) or heterogeneously
                (high_resolution). In the second case, which is the default value for this parameter, in order to determine 
                whether a species appear on a certain network, the presence of any of its habitat in the range of elevations
                considered within each cell of the territory grid is considered. Whereas the more coarse-grained resolution
                only considers the presence of the habitat in the whole cluster, not reaching the cell-in-the-grid detail.
                Since this is only possible when considering the territory in clusters, this parameter has effects only
                when the argument *clustering* is set to True. 
        """
        #for every year we obtain all the species seen on that year and the number of times they were seen
        self.interactions, predators, preys = self.dr.get_interactions()
        self.seasons = self.dr.read_seasonality()
        self.biomasses = self.dr.read_species_mass()
        
        #self.obs_mismatchs = dict()
        
        self.nets.clear()
        #here we set up the data structures that will help us obtain different nets
        #for different elevations based on the habitats present on them
        self.habitats = self.dr.read_habitats()
        self.habs_elevs, min_alt, max_alt = self.dr.read_habitats_min_max_elevations()
        
        offset = (max_alt - min_alt)/ELEVATION_SLOTS
        
        min_interval = min_alt 
        max_interval = min_alt + offset
        self.elevation_intervals = dict()
        for i in range(0,ELEVATION_SLOTS):
            self.elevation_intervals[i] = [min_interval, max_interval]
            min_interval = max_interval 
            max_interval = min_interval + offset
        
        #if we are doing clustering we need to read additional information from the database
        #to determine what are the squares in each habitat cluster and in which species
        #are present in each of the available quads
        #this feature makes the study georeferenced
        if clustering:
            clusters = self.dr.read_habitat_clusters()
            species_sqrs = self.dr.read_species_per_grid()
            
        ##here we construct the interaction networks from the species/year records and the possible interactions
        for current_yr in xrange(START_YEAR, END_YEAR+1, INTERVAL):
            nets_current_year = dict()
            
            if(INTERVAL < 1):
                break
            
            
#            while not observations.has_key(current_yr):
#                current_yr += 1
            
            #with this method we obtain the species that are present in the period being considered (current_yr + INTERVAL)
            #and also the set of years that compose the period
            species_period, years = self.get_species_period(observations, current_yr)
            
            species_temp = set(species_period.keys())
            #this is done for assuring that we are considering only species that are going to be in the final network
            #(because there is actually a trophic relationship reported between them)
            #this step is supposed to save us some time later on
            species_temp = species_temp & (set(predators) | set(preys))
            
            if current_yr == 2000:
                print 'species =', species_temp
            
            
            
            if clustering: 
                for cls in clusters.keys():
                    nets_current_cluster = dict()
                    
                    print '\n current cluster = ', cls, '\n'
                    
                    sqrs_in_cluster = clusters[cls]
                    species_cluster = set()
                    
                    grids_per_species = dict()
                    
                    for y in years:
                        sps_yr = species_sqrs[y]
                        for sp in species_temp:
                            if sp in sps_yr.keys():
                                common_grids = sqrs_in_cluster & sps_yr[sp] 
                                if len(common_grids) > 0:
                                    species_cluster.add(sp)
                                    
                                    grids_per_species[sp] = common_grids  
                    
                    if high_res:
                        self.create_nets_for_elevations_and_grids(species_cluster, nets_current_cluster, grids_per_species)
                    else:
                        self.create_nets_for_elevations(species_cluster, nets_current_cluster)
                    
                    nets_current_year[cls] = nets_current_cluster
                        
            else:  
                self.create_nets_for_elevations(species_temp, nets_current_year)
                
            self.nets[current_yr] = nets_current_year
            
            
        #return self.obs_mismatchs
        return   
            
    def create_nets_for_elevations_and_grids(self, species_temp, nets_dict, grids_per_species):
        """
        Given a set of elevation intervals, which are obtained in advance and kept as an attribute of the class, and the
        distribution of species by grid in the georeferenced territory being consirdered, (given as an argument
        in *grids_per_species*) this method obtains the set of species that are available in those elevations according 
        to the habitats in which they live and the elevations reported for those habitats in each one of the grids
        (considered one by one) that compose the spatial area being considered, i.e. a species will only be added
        to the set of species if any of the grids where it has been seen presents any of the habitats in which it lives
        and the elevation in which that habitat is found within that grid is between the limits of the elevation
        interval being considered.
        
        After obtaining the species that are present in a certain altitude range it calls the method 
        :meth:`create_nets_from_seasonality` to obtain the three networks that can be calculated for the different 
        seasonality values and which contain the species obtained beforehand. 
        
        At the end these obtained networks are stored in a dictionary (received as an argument) classified
        according to different seasonalities. This is done for each elevation, creating in this way a 2 level dictionary.
        
        This is a finer grained version of the method :meth:`create_nets_for_elevations`
        
        call signature::
    
            create_nets_for_elevations_and_grids(species_temp, nets_dict, grids_per_species)
        
        *species_temp*
            a set containing the species that are going to be considered and classified according to altitudes for the 
            creation of the corresponding networks
        *nets_dict*
            a reference to the dictionary where the networks obtained for each seasonality and for each elevation are going
            to be stored after the calculations
        *grids_per_species*
            a dictionary containing the set of grids in which each species has been reported to be seen, where the species
            names are the dictionary keys
        
        .. seealso::
    
            :meth:`create_nets_from_seasonality`
                for the way in which networks are obtained based on the seasonalities of the composing species
            :meth:`create_nets_for_elevations`
                for a simpler version of this method (not considering grids distributions of the species)
        """
        grids_features = self.dr.read_grids_features()
        
        for k in self.elevation_intervals.keys():
            nets_current_alt = dict()
            current_min = self.elevation_intervals[k][0]
            current_max = self.elevation_intervals[k][1]
            
            print current_min, current_max
            species = set()
            common_habitats = set()
            for spec in species_temp:
                if not self.habitats.has_key(spec):
                    print 'Species '+spec+' does not have habitat'
                    continue
                
                sp_habs = self.habitats[spec]
                
                sp_grids = grids_per_species[spec]
                
                for grid in sp_grids:
                    
                    habitats_grid = set(grids_features[grid].keys())
                
                    common_habitats.clear()
                    common_habitats = set(sp_habs) & set(habitats_grid)
                    
                    if len(common_habitats) > 0:
                        for habitat in common_habitats:        
                            current_half = current_min + ((current_max-current_min)/2)
                            
                            #here we determine whether the current habitat of the current species lies
                            hb_sp_min = grids_features[grid][habitat]['min_alt']
                            hb_sp_max = grids_features[grid][habitat]['max_alt']  
                            if ( hb_sp_min > current_min and hb_sp_max < current_max) or (hb_sp_min < current_half and hb_sp_min > current_min) or (hb_sp_max < current_max and hb_sp_max > current_half):
                                species.add(spec)
                                break
                            
                    else:
                        print 'Species ', spec,' was seen on grid ', grid, ' which does not contains any of its habitats'
                        
                        #uncomment this piece of code for recording the species that were seen in grids where none of their habitats are present
                        
#                        if not self.obs_mismatchs.has_key(spec):
#                            self.obs_mismatchs[spec] = dict()
#                        
#                        if not self.obs_mismatchs[spec].has_key('grids'):
#                            self.obs_mismatchs[spec]['grids'] = set()
#                        
#                        self.obs_mismatchs[spec]['grids'].add(grid)
#                        
#                        if not self.obs_mismatchs[spec].has_key('habitats'):
#                            self.obs_mismatchs[spec]['habitats'] = sp_habs
#                        
            net_ss, net_aw, net_ay, net_all = self.create_nets_from_seasonality(species)
        
            nets_current_alt['ss'] = net_ss
            nets_current_alt['aw'] = net_aw
            nets_current_alt['ay'] = net_ay
            nets_current_alt['all'] = net_all
        
            nets_dict[k] = nets_current_alt
            

    def create_whole_network(self, species):
        """
        By looking at the interactions and seasonalities between the set of species given as argument in *'species'*
        (the interactions and seasonalities are retrieved in advance and are stored as attributes of this class),
        this method creates an instance of the object :class:`~web.Network` and builds a network representation
        of the trophic interactions between those species. Interactions are only created when the seasonalities of the
        species coincide. Additionally, this method also assigns two attributes
        to each node/species in the network: 'biomass' and 'group', which are useful species attributes that are
        used for different network features. 
        
        call signature::
    
            create_whole_network(species)
        
        *species*
            the set of species that are going to be considered (through their trophic interactions and seasonalities)
            for the creation of the ecological networks produced by this method
        
        It returns an instance of the object :class:`~web.Network`: corresponding to the overall interaction network
        in the database
        """
        net = Network()
        self.interactions, predators, preys = self.dr.get_interactions()
        self.seasons = self.dr.read_seasonality()
        self.biomasses = self.dr.read_species_mass()
        
        for sp in species:
            ##for each species we look at its interactions
            if self.interactions.has_key(sp):
                preys_inter = self.interactions[sp]
                ##we iterate over the interactions to consider each separately
                for kp in preys_inter.keys():
                    food = preys_inter[kp]
                    for f in food:
                        ##depending on what kind of interaction they are involved the appropriate node is added to the network
                        
                        ##we also confirm the species and its food co-occurr in the same season
                        sea_sp = self.seasons[sp]
                        
                        if self.seasons.has_key(f):
                            sea_f = self.seasons[f]
                        else:
                            if f in species or kp == 3:
                                print 'Node '+ f + ' does not have seasonality' 
                            sea_f = 'yr'
                        
                        if sea_sp == 'ay':
                            sea_sp = 'yr'
                        if sea_f == 'ay':
                            sea_f = 'yr'
                        
                        if (kp == 1 and (f in species)):
                            if sea_sp == 'yr':
                                net.add_edge(f,sp)                                
                            elif sea_sp == 'ss':
                                if sea_f == 'ss' or sea_f == 'yr':
                                    net.add_edge(f,sp)
                            elif sea_sp == 'aw':
                                if sea_f == 'aw' or sea_f == 'yr':
                                    net.add_edge(f,sp)
                        elif kp == 3:
                            if sea_sp == 'aw' or sea_sp == 'yr':
                                net.add_edge(f,sp)
                                if self.interactions.has_key(f):
                                    for ff in self.interactions[f][3]:
                                        net.add_edge(ff, f)
                            if sea_sp == 'ss' or sea_sp == 'yr':
                                net.add_edge(f,sp)
                                if self.interactions.has_key(f):
                                    for ff in self.interactions[f][3]:
                                        net.add_edge(ff, f)
                                        
        for v in net:
            try:
                net.node[v]['biomass'] = self.biomasses[v]['biomass']
            except:
                print 'Node '+ v + ' does not have biomass'
                net.node[v]['biomass'] = ' '
            try:
                net.node[v]['group'] = self.biomasses[v]['group']
            except:
                print 'Node '+ v + ' does not have group'
                net.node[v]['group'] = 'N/A'
        
        return net
     
            
    def create_networks_from_cell_and_species(self, grid, species_temp):
        """
        This method creates networks representations of the species in the *observations* data set in more or less the same way
        as detailed in :meth:`create_networks_from_data` but it also considers clusters of regions within the territory
        of study to create networks for different areas within it (clusters), based on the elevations, habitats for classifying 
        the species in order to obtain different networks for different elevations within each cluster of the territory.
        These clusters are defined in advance and are retrieved by the DataReader of this class through the method
        :meth:`~database_reader.read_habitat_clusters` which divides the area of study in georeferenced sections that
        possess many habitats. This results in a combinatorial number of networks in each cluster, according to the 
        habitats present in each cluster; and by elevation, according to the habitats encountered at each elevation
        within each cluster of land.
        Once the composing species based on clusters (or not) are obtained it calls the method :meth:`create_nets_for_elevations`
        in order to obtain the networks for each of the considered elevations.
        At the beginning of this method the information about species seasonalities, biomass, interactions and
        the habitats and their corresponding elevations within each grid are retrieved from the data reader and kept in 
        data structures owned either by this method or the class.
        
        call signature::
    
            create_networks_from_data_and_habitats(observations, clustering=False, high_res=True)
        
        *observations*
            a dictionary with the observations of each species classified by years (the keys of the dictionary)
        
        Optional keyword arguments:
            *clustering*
                a boolean stating whether the cluster classification of the territory is going to be considered for
                the production of the networks, its default value is False, in which case the networks are constructed
                only based on elevations and seasonalities per year.
            *high_res*
                a boolean that determines whether to consider each cluster homogeneously (low resolution) or heterogeneously
                (high_resolution). In the second case, which is the default value for this parameter, in order to determine 
                whether a species appear on a certain network, the presence of any of its habitat in the range of elevations
                considered within each cell of the territory grid is considered. Whereas the more coarse-grained resolution
                only considers the presence of the habitat in the whole cluster, not reaching the cell-in-the-grid detail.
                Since this is only possible when considering the territory in clusters, this parameter has effects only
                when the argument *clustering* is set to True. 
        """
        #for every year we obtain all the species seen on that year and the number of times they were seen
        self.interactions, predators, preys = self.dr.get_interactions()
        self.seasons = self.dr.read_seasonality()
        self.biomasses = self.dr.read_species_mass()
       
        #self.nets.clear()
        nets_current_grid = dict()
        #here we set up the data structures that will help us obtain different nets
        #for different elevations based on the habitats present on them
        self.habitats = self.dr.read_habitats()
        self.habs_elevs, min_alt, max_alt = self.dr.read_habitats_min_max_elevations()
       
        #this is done for assuring that we are considering only species that are going to be in the final network
        #(because there is actually a trophic relationship reported between them)
        #this step is supposed to save us some time later on
        species_temp = species_temp & (set(predators) | set(preys))
        
        
        grids_features = self.dr.read_grids_features()
        habitats_grid = set(grids_features[grid].keys())
        
#         print grids_features
        
        species = set()
        common_habitats = set()
        for spec in species_temp:
            if not self.habitats.has_key(spec):
                print 'Species '+spec+' does not have habitat'
                continue
            
            sp_habs = self.habitats[spec]
            
            common_habitats.clear()
            common_habitats = set(sp_habs) & set(habitats_grid)
            
#             print common_habitats
            
            if len(common_habitats) > 0:
#                     for habitat in common_habitats:        
#                         current_half = current_min + ((current_max-current_min)/2)
#                         
#                         #here we determine whether the current habitat of the current species lies
#                         hb_sp_min = grids_features[grid][habitat]['min_alt']
#                         hb_sp_max = grids_features[grid][habitat]['max_alt']  
#                         if ( hb_sp_min > current_min and hb_sp_max < current_max) or (hb_sp_min < current_half and hb_sp_min > current_min) or (hb_sp_max < current_max and hb_sp_max > current_half):
                species.add(spec)
#                             break
                    
            else:
                print 'Species ', spec,' was seen on grid ', grid, ' which does not contains any of its habitats'
       
#                        
        net_ss, net_aw, net_ay, net_all = self.create_nets_from_seasonality(species)
    
        nets_current_grid['ss'] = net_ss
        nets_current_grid['aw'] = net_aw
        nets_current_grid['ay'] = net_ay
        nets_current_grid['all'] = net_all
        
                
        self.nets[grid] = nets_current_grid
            
        
        return   
            