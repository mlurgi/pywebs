
import MySQLdb
import math

import csv

from config import DB_HOST, DB_USER, DB_PASSWD, DB_NAME
from config import START_YEAR, END_YEAR, SPS_SUFFIXES, CLUSTERING
from config import DATA_AGGREGATION, DATA_AGGREGATION_INTERVAL

class DatabaseReader():
    """
    The DatabaseReader implements an interface with a MySQL relational database for retrieving data about the studied
    species and other features related to the ecosystem and provides methods for providing access to particular
    objects of the database. 
    
    For the implementation of the connections and interface with the database it takes advantage of the functionalities
    provided by the MySQLdb library for python. 
    """
    def __init__(self):
        self.taxons = None
        self.species_to_remove = None
        self.seasons = None
        self.biomasses = None
        self.habitats = None
        
        self.children_species = None
        
        self.interactions = None
        self.predators = None
        self.preys = None
        
        self.grids = None
        self.study_grids = None
        self.grids_features = None

    def get_db_instance(self):
        """
        Gets an instance of the database
        """
        db = MySQLdb.connect(host=DB_HOST, # your host, usually localhost
                             user=DB_USER, # your username
                             passwd=DB_PASSWD, # your password
                             db=DB_NAME) # name of the data base
        return db

    def get_interactions(self):
        """
        Returns the trophic interactions between the species obtained from the database, the predators and the preys.
        This information (the interactions) is kept in data structures that belong to this class.
        """
        if self.interactions == None:
            self.read_interactions()
            
        return self.interactions, self.predators, self.preys
    
    def read_taxons(self, new=True):
        """
        Retrieves the species from the database, including their biomass, seasonality and group they belong to.
        """
        self.taxons = dict()
        self.seasons = dict()
        self.biomasses = dict()
        self.species_to_remove = set()
        
        db = self.get_db_instance()

        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM species")
        species = cur.fetchall()
        
        cur.execute("SELECT * FROM seasonality")
        db_seasonalities = cur.fetchall()
        
        cur.execute("SELECT * FROM species_group")
        db_groups = cur.fetchall()
        
        cur.close()
        db.close()
        
        dict_seasons = dict()
        for row in db_seasonalities:
            dict_seasons[row['id']] = row['seasonality_name']
        
        dict_groups = dict()
        for row in db_groups:
            dict_groups[row['id']] = row['group_name']
        
        for row in species:
            sp_name = row['species_name']
            self.taxons[row['id']] = sp_name
            
            season_id = row['seasonality_id']
            if season_id == None:
                self.seasons[sp_name] = None
            else:
                self.seasons[sp_name] = dict_seasons[season_id]
                
            self.biomasses[sp_name] = dict()
            mass = row['species_mass']
            if mass == None:
                self.biomasses[sp_name]['biomass'] = '0.0'
            else:
                self.biomasses[sp_name]['biomass'] = mass
             
            group_id = row['group_id']
            if group_id == None:
                self.biomasses[sp_name]['group'] = 'N/A'
            else:
                self.biomasses[sp_name]['group'] = dict_groups[group_id]
            
            if row['remove'] == 1:
                self.species_to_remove.add(sp_name)
            
        return
    
    def read_species_mass(self):
        """
        Returns a reference to the data structure where the biomasses of the species are stored, which are read in the :meth:`read_taxons` method.
        """
        if self.biomasses == None:
            self.read_taxons()
        
        return self.biomasses
    
    ##seasonality - we obtain the group of species that occurr together based on their seasonality
    def read_seasonality(self):
        """
        Returns a reference to the data structure where the seasonalities of the species are stored, which are read in the :meth:`read_taxons` method.
        """
        if self.seasons == None:
            self.read_taxons()
        
        return self.seasons
    
    def read_grids(self):
        """
        This method reads the information about the grids (locations) that are kept in the database
        """
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM locations")
    
        grids = cur.fetchall()
        cur.close()
        db.close()
        
        self.grids = dict()
        for row in grids:
            self.grids[row['id']] = row['location_name']
        
        return
    
    def read_pyrenees_quads(self):
        """
        This method reads exclusively the grids or locations that belong to the study to be performed.
        It should be a subset of the total locations and must be stored in the study_locations table.
        """
        if self.grids == None:
            self.read_grids()
        
        self.study_grids = set()
        
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM study_locations")
        
        db_grids = cur.fetchall()
        cur.close()
        db.close()
        
        for row in db_grids:
            self.study_grids.add(self.grids[row['location_id']])
            
        return
    
    
    def read_interactions(self):
        """
        This method reads the interactions between the species from the database.
        """
        self.interactions = dict()
        self.predators = set()
        self.preys = set()
        self.children_species = set()
        
        if self.taxons == None:
            self.read_taxons()
        
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM interactions")
        
        db_inters = cur.fetchall()
        cur.close()
        db.close()
        
        
        for row in db_inters:
            predator = self.taxons[row['predator_id']]
            prey = self.taxons[row['prey_id']]
            group = int(row['type_id'])
            if  predator in self.interactions.keys():
                inter_sp = self.interactions[predator]      
                if group in inter_sp.keys():
                    if not prey in inter_sp[group]:
                        inter_sp[group].append(prey)
                else:
                    inter_sp[group] = [prey]
            else:
                self.interactions[predator] = dict({group:[prey]})
            
            self.predators.add(predator)
            self.preys.add(prey)
            
            
            for sf in SPS_SUFFIXES:
                if sf in predator:
                    self.children_species.add(predator)
                if sf in prey:
                    self.children_species.add(prey)
        
        return
    
    
    #this method was re-implemented in its current form (information retrieved from database) for
    #historical reasons, but it is deprecated. In its place read_species_per_grid should be employed (see below)
    def read_species(self, geo_ref=True):
        """
        This method reads the observations made for each species by year and by location/grid
        
        Optional keyword argument:
            *geo_ref*
                a boolean stating whether to report all of the observations available or only
                those belonging to the grids of the current study (which are stored in
                `self.study_grids` and retrieved using the method :meth:`read_pyrenees_quads`.
                The default value for this argument is True.
        """
        #here we read the species data from database
        if self.taxons == None:
            self.read_taxons()
            
        if self.grids == None:
            self.read_grids()

        if self.study_grids == None:
            self.read_pyrenees_quads()
        
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM observations")
        
        db_observations = cur.fetchall()
        cur.close()
        db.close()
        
        years_sps = dict()
        #we read the data from the dataset and store it in memory
        for row in db_observations:
            current_grid = self.grids[row['location_id']]
            if (geo_ref and current_grid in self.study_grids) or not geo_ref:
                if not self.taxons.has_key(row['species_id']):
                    continue
                sp_name = self.taxons[row['species_id']]
                if sp_name in years_sps.keys():
                    years_sp = years_sps[sp_name]
                    try:
                        yr = int(row['year'])
                        if yr in years_sp.keys():
                            years_sp[yr] += 1
                        else: 
                            if yr >= START_YEAR and yr <= END_YEAR:
                                years_sp[yr] = 1
                    except:
                        continue
                else:
                    try:
                        year = int(row['year'])
                        if year >= START_YEAR and year <= END_YEAR:
                            years_sps[sp_name] = dict({year : 1})
                        else:
                            years_sps[sp_name] = dict()
                    except:
                        continue
        
        if self.interactions == None:
            self.read_interactions()
                
        for sp_child in self.children_species:
            parent_sp = sp_child.rpartition(' ')[0]
            if years_sps.has_key(parent_sp):
                years_sps[sp_child] = years_sps[parent_sp].copy()
        
        
        ### the folowing piece of code is used to add the observations that were detected missing from the BIOCAT's database
        ### according to which several species (those on the list below) got extinct during the 90's and 00's
        ### this was done in this way because it was found a better solution than adding this list to the database
        species_to_complete = ['Apodemus flavicollis', 'Apodemus sylvaticus', 'Arvicola sapidus', 'Arvicola terrestris', 
                               'Coronella girondica', 'Crocidura russula', 'Eliomys quercinus', 'Erinaceus europaeus', 
                               'Genetta genetta', 'Genetta genetta young', 'Glis glis', 'Hieraaetus fasciatus',
                               'Lepus europaeus', 'Lepus europaeus young', 'Malpolon monspessulanus', 'Meles meles', 
                               'Meles meles young', 'Microtus agrestis', 'Microtus arvalis', 'Microtus duodecimcostatus', 
                               'Microtus nivalis', 'Mus musculus', 'Mus spretus', 'Mustela nivalis', 'Mustela putorius', 
                               'Myodes glareolus', 'Neomys fodiens', 'Oryctolagus cuniculus', 'Oryctolagus cuniculus young', 
                               'Rallus aquaticus', 'Rattus norvegicus', 'Rattus rattus', 'Sorex araneus', 'Sorex coronatus',
                               'Sorex minutus', 'Suncus etruscus', 'Talpa europaea', 'Vulpes vulpes', 
                               'Neomys anomalus', 'Rhinechis scalaris', 'Tarentola mauritanica', 'Chalcides striatus', 
                               'Mustela vison', 'Barbus graellsii', 'Pelobates cultripes eggs', 'Triturus marmoratus eggs', 
                               'Bufo calamita', 'Triturus marmoratus', 'Pelobates cultripes', 'Bufo calamita eggs', 
                               'Psammodromus algirus', 'Falco naumanni', 'Timon lepidus', 'Aegypius monachus', 
                               'Acrocephalus schoenobaenus', 'Natrix natrix', 'Tringa nebularia',
                               'Ardea cinerea', 'Lacerta agilis', 'Sylvia borin', 'Numenius arquata', 'Lissotriton helveticus']         
        
        
        
        years_to_complete = [1995,1996,1997,1998,1999,2000,2001,2002]
        
        
        for sp in species_to_complete:
            years_sp = years_sps[sp]
            extinction_year = years_to_complete[0]
            for y in years_to_complete:
                if years_sp.has_key(y):
                    extinction_year = y+1
            
            print sp, 'extinction year = ', extinction_year
                        
            missing_years = range(extinction_year-3, 2003)
            for y in missing_years:
                 if not years_sp.has_key(y):
                     years_sp[y] = 1
        
        
        #this piece of code is used to duplicate the entries for species that were not
        #reported in some years but the gap between consecutive observations is less
        #than a certain number of years
        if DATA_AGGREGATION:
            added_years = dict()
            species = years_sps.keys()
            for sp in species:
                years = sorted(years_sps[sp].keys())
                gap = 0
                for i in xrange(len(years)-1):
                    gap = (years[i+1] - years[i])-1
                    if gap > 0 and gap <= DATA_AGGREGATION_INTERVAL:
                        min = years[i]
                        max = years[i+1]
                        half_gap = int(math.floor(gap/2))
                        max_half_gap = max+half_gap
#                        if max_half_gap > END_YEAR:
#                            max_half_gap = END_YEAR
                        for new_year in xrange(min-half_gap, max_half_gap+1):
                            if not years_sps[sp].has_key(new_year):
                                years_sps[sp][new_year] = 1
                                
                                if added_years.has_key(sp):
                                    added_years[sp].add(new_year)
                                else:
                                    added_years[sp] = set([new_year])
            
            #this is for saving the modifications done over the original
            #dataset on a file
            header_names = ['species','years']
            out = csv.DictWriter(open('../output/added_in_years.csv', 'w'), header_names, delimiter=',')
    
            #we write the headers
            headers_dict = dict()
            for n in header_names:
                headers_dict[n] = n
    
            out.writerow(headers_dict)
            
            out_row = dict()
            
            for sp in sorted(added_years.keys()):
                out_row['species'] = sp
                out_row['years'] = sorted(added_years[sp])
                
                out.writerow(out_row)
        
        print 'Number of documented species from', START_YEAR, 'to', END_YEAR,'=', len(years_sps)
        return years_sps
    
    def read_species_per_grid(self):
        """
        This methods reads the distribution of species per grid/location for each year.
        """
        if self.taxons == None:
            self.read_taxons(new=True)
        
        if self.grids == None:
            self.read_grids()
        
        if self.study_grids == None:
            self.read_pyrenees_quads()
        
        if self.interactions == None:
            self.read_interactions()
        
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM observations")
        
        db_observations = cur.fetchall()
        cur.close()
        db.close()
        
        years_sps = dict()
        for row in db_observations:
            sq = self.grids[row['location_id']]
            if sq in self.study_grids:
                if not self.taxons.has_key(row['species_id']):
                    continue
                sp_name = self.taxons[row['species_id']]
                
                current_children = set()
                for sp_child in self.children_species:
                    parent_sp = sp_child.rpartition(' ')[0]
                    
                    if parent_sp == sp_name:
                        current_children.add(sp_child)
                                        
                try:
                    year = int(row['year'])
                    if year in years_sps.keys():
                        species_year = years_sps[year]
                        if sp_name in species_year:
                            species_year[sp_name].add(sq)
                        else:
                            species_year[sp_name] = set([sq])
                        
                        for sp_child in current_children:
                            if sp_child in species_year:
                                species_year[sp_child].add(sq)
                            else:
                                species_year[sp_child] = set([sq])
                            
                    else:
                        if year >= START_YEAR and year <= END_YEAR:
                            years_sps[year] = dict({sp_name : set([sq])})
                            for sp_child in current_children:
                                years_sps[year] = dict({sp_child : set([sq])})
                                
                except:
                    continue
            
        ### the folowing piece of code is used to add the observations that were detected missing from the BIOCAT's database
        ### according to which several species (those on the list below) got extinct during the 90's and 00's
        ### this was done in this way because it was found a better solution than adding this list to the database
        species_to_complete = ['Apodemus flavicollis', 'Apodemus sylvaticus', 'Arvicola sapidus', 'Arvicola terrestris', 
                               'Coronella girondica', 'Crocidura russula', 'Eliomys quercinus', 'Erinaceus europaeus', 
                               'Genetta genetta', 'Genetta genetta young', 'Glis glis', 'Hieraaetus fasciatus',
                               'Lepus europaeus', 'Lepus europaeus young', 'Malpolon monspessulanus', 'Meles meles', 
                               'Meles meles young', 'Microtus agrestis', 'Microtus arvalis', 'Microtus duodecimcostatus', 
                               'Microtus nivalis', 'Mus musculus', 'Mus spretus', 'Mustela nivalis', 'Mustela putorius', 
                               'Myodes glareolus', 'Neomys fodiens', 'Oryctolagus cuniculus', 'Oryctolagus cuniculus young', 
                               'Rallus aquaticus', 'Rattus norvegicus', 'Rattus rattus', 'Sorex araneus', 'Sorex coronatus',
                               'Sorex minutus', 'Suncus etruscus', 'Talpa europaea', 'Vulpes vulpes', 
                               'Neomys anomalus', 'Rhinechis scalaris', 'Tarentola mauritanica', 'Chalcides striatus', 
                               'Mustela vison', 'Barbus graellsii', 'Pelobates cultripes eggs', 'Triturus marmoratus eggs', 
                               'Bufo calamita', 'Triturus marmoratus', 'Pelobates cultripes', 'Bufo calamita eggs', 
                               'Psammodromus algirus', 'Falco naumanni', 'Timon lepidus', 'Aegypius monachus', 
                               'Acrocephalus schoenobaenus', 'Natrix natrix', 'Tringa nebularia',
                               'Ardea cinerea', 'Lacerta agilis', 'Sylvia borin', 'Numenius arquata', 'Lissotriton helveticus']         
        
        
        
        years_to_complete = [1995,1996,1997,1998,1999,2000,2001,2002]


        added_years = dict()
        extinction_years = dict.fromkeys(species_to_complete)
        for sp in species_to_complete:
            extinction_year = years_to_complete[0]
            new_grids = set()
            for y in years_to_complete:
                if years_sps[y].has_key(sp):
                    extinction_year = y+1
            
            print sp, 'extinction year = ', extinction_year
            extinction_years[sp] = extinction_year
            
            for y in range(extinction_year-3, extinction_year):
                if years_sps[y].has_key(sp):
                    new_grids |= years_sps[y][sp]
            
            missing_years = range(extinction_year-3, 2003)
            for y in missing_years:
                years_sps[y][sp] = new_grids
            
                print 'species', sp, 'year', y, 'new grids', new_grids
            
                if not added_years.has_key(sp):
                    added_years[sp] = dict({y: new_grids})
                else:
                    if not added_years[sp].has_key(y):
                        added_years[sp][y] = new_grids
                    
        #this piece of code is used to duplicate the entries for species that were not
        #reported in some years but the gap between consecutive observations is less
        #than a certain number of years
        if DATA_AGGREGATION:

            new_years_sps = dict()
            years = years_sps.keys()
            for y in years:
                species = years_sps[y].keys()
                
                for sp in species:
                    grids = years_sps[y][sp]
                    
                    for g in grids:
                        if not new_years_sps.has_key(sp):
                            new_years_sps[sp] = dict({g: set([y])})
                        else:
                            if new_years_sps[sp].has_key(g):
                                new_years_sps[sp][g].add(y)
                            else:
                                new_years_sps[sp][g] = set([y])
            
            species = new_years_sps.keys()
            for sp in species:
                grids = new_years_sps[sp].keys()
                for g in grids:
                    years = sorted(new_years_sps[sp][g])
                    gap = 0
                    for i in xrange(len(years)-1):
                        gap = (years[i+1] - years[i])-1
                        if gap > 0 and gap <= DATA_AGGREGATION_INTERVAL:
                            min = years[i]
                            max = years[i+1]
                            half_gap = int(math.floor(gap/2))
                            max_half_gap = max+half_gap
                            if max_half_gap > END_YEAR:
                                max_half_gap = END_YEAR
                            for new_year in xrange(min-half_gap, max_half_gap+1):
                                added = False
                                if not years_sps.has_key(new_year):
                                    years_sps[new_year] = dict({sp: set([g])})
                                    added = True
                                    #print 'added :', g, 'in year: ', new_year, ' for species: ', sp
                                else:
                                    if not years_sps[new_year].has_key(sp): 
                                        years_sps[new_year][sp] = set([g])
                                        added = True
                                        #print 'added :', g, 'in year: ', new_year, ' for species: ', sp
                                    else:
                                        if g not in years_sps[new_year][sp]:
                                            years_sps[new_year][sp].add(g)
                                            added = True
                                            #print 'added :', g, 'in year: ', new_year, ' for species: ', sp
                                    
                                
                                if added:
                                    if not added_years.has_key(sp):
                                        added_years[sp] = dict({new_year: set([g])})
                                    else:
                                        if not added_years[sp].has_key(new_year):
                                            added_years[sp][new_year] = set([g])
                                        else:
                                            added_years[sp][new_year].add(g)
                                    
                                            
            
            #this is for saving the modifications done over the original
            #dataset on a file
            header_names = ['species','type','year','grids']
            out = csv.DictWriter(open('../output/added_in_grids_per_year.csv', 'w'), header_names, delimiter=',')
    
            #we write the headers
            headers_dict = dict()
            for n in header_names:
                headers_dict[n] = n
    
            out.writerow(headers_dict)
            
            out_row = dict()
            
            for sp in sorted(added_years.keys()):
                out_row['species'] = sp
                for y in sorted(added_years[sp].keys()):
                    
                    if sp in species_to_complete and (y >= extinction_years[sp] and y < 2003):
                        out_row['type'] = 0
                    else:
                        out_row['type'] = 1
                    
                    out_row['year'] = y
                    out_row['grids'] = sorted(added_years[sp][y])
                    
                    out.writerow(out_row)
                    
        return years_sps
    
    
    def read_habitats(self):
        """
        This method obtains the sets of habitats for each species from the database
        """
        if self.taxons == None:
            self.read_taxons()
        
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM species_habitats")
        
        db_sp_habitats = cur.fetchall()
        cur.close()
        db.close()
        
        self.habitats = dict()
        for row in db_sp_habitats:
            species = self.taxons[row['species_id']]
            habitat_code = row['habitat_id']
            if self.habitats.has_key(species):
                self.habitats[species].add(habitat_code)
            else:
                self.habitats[species] = set([habitat_code])
                
        return self.habitats
    
    def read_habitats_min_max_elevations(self):
        """
        This method obtains the minimum and maximum altitudes for each habitat in the database. It also obtains the overall min and max altitudes.
        """
        hab_elevations = dict()
        self.grids_features = dict()
        
        if self.grids == None:
            self.read_grids()
        
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM locations_features")
        
        db_grids_features = cur.fetchall()
        cur.close()
        db.close()
        
        min_alt = None
        max_alt = None
        
        for row in db_grids_features:
            #this part of the loop extracts the features of the grids
            grid_id = self.grids[row['location_id']]
            if not self.grids_features.has_key(grid_id):
                self.grids_features[grid_id] = dict()
            
            habitat_id = row['habitat_id']
            if not self.grids_features[grid_id].has_key(habitat_id):
                self.grids_features[grid_id][habitat_id] = dict()
                
            self.grids_features[grid_id][habitat_id]['area'] = row['area']
            self.grids_features[grid_id][habitat_id]['min_alt'] = row['min_altitude']
            self.grids_features[grid_id][habitat_id]['max_alt'] = row['max_altitude']
            
            
            #while this part of the loop determines the max and min altitudes for
            #each habitat and over all the grids
            tmp_min = float(row['min_altitude'])
            tmp_max = float(row['max_altitude'])
            
            if tmp_min != 0.0 and (min_alt == None or tmp_min < min_alt):
                min_alt = tmp_min
            if max_alt == None or tmp_max > max_alt:
                max_alt = tmp_max
            
            hab_code = row['habitat_id']
            if hab_elevations.has_key(hab_code):
                if hab_elevations[hab_code][0] == 0.0 or (tmp_min != 0.0 and tmp_min < hab_elevations[hab_code][0]):
                    hab_elevations[hab_code][0] = tmp_min
                if tmp_max > hab_elevations[hab_code][1]:
                    hab_elevations[hab_code][1] = tmp_max
            else:
                hab_elevations[hab_code] = [tmp_min,tmp_max]
            
        return hab_elevations, min_alt, max_alt
    
    
    def read_grids_features(self):
        """
        Obtains the features of the grids
        """
        if self.grids_features == None:
            self.read_habitats_min_max_elevations()
            
        return self.grids_features
    
    def read_habitat_clusters(self):
        """
        Retrieves the clusters of grids, specifying which grids belong to which cluster
        """
        if self.grids == None:
            self.read_grids()
        
        clusters = dict()
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM clusters")
        db_clusters = cur.fetchall()
        
        cur.close()
        db.close()
        
        for row in db_clusters:
            cluster = int(row['cluster_'+str(CLUSTERING)])
            grid = self.grids[row['location_id']].strip()
            if clusters.has_key(cluster):
                clusters[cluster].add(grid)
            else:
                clusters[cluster] = set([grid])
            
        return clusters
        
    
    def get_species_observations(self):
        """
        Obtains the number of observations for each species per year
        """
        years_sps = self.read_species()
        
        #observations_yr = dict{year: dict{sp_name: number_of_observations}}
        observations_yr = dict()
        for k in years_sps.keys():
            years = years_sps[k]
            for y in years.keys():
                if observations_yr.has_key(y):
                    temp = observations_yr[y]
                    temp[k] = years[y]
                else:
                    observations_yr[y] = dict({k : years[y]})
        
        return observations_yr, years_sps

    def read_habitats_names(self):
        """
        Reads the names of the habitats
        """
        habitats = dict()
        db = self.get_db_instance()
        cur = db.cursor(MySQLdb.cursors.DictCursor)
        cur.execute("SELECT * FROM habitats")
        db_habitats = cur.fetchall()
        
        cur.close()
        db.close()
        
        for row in db_habitats:
            habitats[row['id']] = row['habitat_name']
    
        return habitats
    
    def get_species_to_remove(self):
        """
        Returns a reference to the data structure where the species to be removed from the network are stored, 
        which are read in the :meth:`read_taxons` method.
        """
        if self.species_to_remove == None:
            self.read_taxons()
            
        return self.species_to_remove
    
    
#if __name__ == '__main__':
#    dr = DatabaseReader()
#    print dr.read_habitat_clusters()

