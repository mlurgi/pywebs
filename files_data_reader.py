

import csv
from config import START_YEAR, END_YEAR, SP_SUFFIX

class FilesDataReader():
    """
    The FilesDataReader implements an interface with a series of csv files for retrieving data about the studied
    species and other features related to the ecosystem and provides methods for providing access to particular
    objects of the dataset.
    
    For the implementation of the files parsing we use the csv module of the Python API. 
    """
    def __init__(self):
        self.interactions = None
        self.predators = None
        self.preys = None
        self.youngs = None
        self.taxons = None
    
    def get_interactions(self):
        """
        Returns the trophic interactions between the species obtained from the source files, the predators and the preys.
        This information (the interactions) is kept in data structures that belong to this class.
        """
        if self.interactions == None:
            self.read_interactions()
            
        return self.interactions, self.predators, self.preys
    
    def read_taxons(self, new=True):
        """
        Retrieves the species from the database and keeps this info in a dictionary where the codes given for each species
        in the dataset are the keys and the names are the values
        """
        self.taxons = dict()
        
        if new:
            taxs = csv.DictReader(open('../input/TAXONS_VERTE.DEL', 'rb'), delimiter=',', quotechar='"')
            for row in taxs:
                sp_name = row["Genre"].strip()+" "+row["Species"].strip()
                self.taxons[row["Code"]] = sp_name
            
            print 'Number of species in TAXONS_VERTE =', len(self.taxons)
            
        else:
            taxs = csv.DictReader(open('../input/TaxonsExportacio.csv', 'rb'), delimiter='\t', quotechar='"')
            for row in taxons:
                self.taxons[row["CodiSps"]] = row["Sps"]
        
            print 'Number of species in TaxonsExportacio =', len(self.taxons)
        
        return
    
    def read_species(self, new=True, geo_ref=True):
        """
        This method reads the observations made for each species by year and by location/grid.
        
        Optional keyword arguments:
            *new*
                a boolean value stating whether to use the newest data set available. Default to True
            *geo_ref*
                a boolean value stating whether to retrieve only the data belonging to the grids subject of study.
                In this case the ones for the Pyrenees region. Default value = True
        """
        #here we read the species data from biocat database
        if new:
            avail_species = csv.DictReader(open('../input/vertebrats_2.txt', 'rb'), delimiter=',', quotechar='"')
            self.read_taxons(new=True)
        else:
            avail_species = csv.DictReader(open('../input/DadesVertebratsExportacio.csv', 'rb'), delimiter=';', quotechar='"')
            self.read_taxons(new=False)
        
        #additionally, we also include the information from the ICO (Catalan Institute of Ornithology) dataset in order to obtain
        #a better resolver network 
        avail_birds = csv.DictReader(open('../input/ICO_data_in_plots.csv', 'rb'), delimiter=';', quotechar='"')
        
        squares = self.read_pyrenees_quads()
        
        years_sps = dict()
        #we read the data from the dataset and store it in memory
        for row in avail_species:
            if (geo_ref and row["Quad"] in squares) or not geo_ref:
                if not self.taxons.has_key(row["CodiSps"]):
                    continue
                sp_name = self.taxons[row["CodiSps"]]
                if sp_name in years_sps.keys():
                    years_sp = years_sps[sp_name]
                    try:
                        yr = int(row["Any"])
                        if yr in years_sp.keys():
                            years_sp[yr] += 1
                        else: 
                            if yr >= START_YEAR and yr <= END_YEAR:
                                years_sp[yr] = 1
                    except:
                        continue
                else:
                    try:
                        year = int(row["Any"])
                        if year >= START_YEAR and year <= END_YEAR:
                            years_sps[sp_name] = dict({year : 1})
                        else:
                            years_sps[sp_name] = dict()
                    except:
                        continue
        
        #to this dataset we add more information for the bird species (ICO's dataset)
        for row in avail_birds:
            if (geo_ref and row["plot"] in squares) or not geo_ref:
                sp_name = row["species"]
                if sp_name in years_sps.keys():
                    years_sp = years_sps[sp_name]
                    try:
                        yr = int(row["year"])
                        if yr in years_sp.keys():
                            years_sp[yr] += 1
                        else: 
                            if yr >= START_YEAR and yr <= END_YEAR:
                                years_sp[yr] = 1
                    except:
                        continue
                else:
                    try:
                        year = int(row["year"])
                        if year >= START_YEAR and year <= END_YEAR:
                            years_sps[sp_name] = dict({year : 1})
                        else:
                            years_sps[sp_name] = dict()
                    except:
                        continue
        
        
        if self.interactions == None:
            self.read_interactions()
        
        if len(self.youngs) > 0:
            for species in years_sps.keys():
                new_spec = species + ' ' + SP_SUFFIX 
                if new_spec in self.youngs:
                    years_sps[new_spec] = years_sps[species].copy() 
        
        print 'Number of documented species from', START_YEAR, 'to', END_YEAR,'=', len(years_sps)
        return years_sps
    
    def read_species_per_grid(self):
        """
        This methods reads the distribution of species per grid/location for each year.
        """
        #we read the information from the biocat database for the presence of species
        avail_species = csv.DictReader(open('../input/vertebrats_2.txt', 'rb'), delimiter=',', quotechar='"')
        
        #and we also incorporate the information from the ICO's database into our dataset for more detailed bird information
        avail_birds = csv.DictReader(open('../input/ICO_data_in_plots.csv', 'rb'), delimiter=';', quotechar='"')
        
        if self.taxons == None:
            self.read_taxons(new=True)
        squares = self.read_pyrenees_quads()
        
        if self.interactions == None:
            self.read_interactions()
        
        years_sps = dict()
        for row in avail_species:
            sq = str(row["Quad"]) 
            if sq in squares:
                if not self.taxons.has_key(row["CodiSps"]):
                    continue
                sp_name = self.taxons[row["CodiSps"]]
                sp_name_new = sp_name + ' ' + SP_SUFFIX
                ys = False
                if sp_name_new in self.youngs:
                    ys = True
                try:
                    year = int(row['Any'])
                    if year in years_sps.keys():
                        species_year = years_sps[year]
                        if sp_name in species_year:
                            species_year[sp_name].add(sq)
                        else:
                            species_year[sp_name] = set([sq])
                        
                        if ys:
                            if sp_name_new in species_year:
                                species_year[sp_name_new].add(sq)
                            else:
                                species_year[sp_name_new] = set([sq])
                            
                    else:
                        if year >= START_YEAR and year <= END_YEAR:
                            years_sps[year] = dict({sp_name : set([sq])})
                            if ys:
                                years_sps[year] = dict({sp_name_new : set([sq])})
                                
                except:
                    continue
                
        #we do the same thing for the birds (ICO)
        for row in avail_birds:
            sq = str(row["plot"]) 
            if sq in squares:
                sp_name = row["species"]
                sp_name_new = sp_name + ' ' + SP_SUFFIX
                ys = False
                if sp_name_new in self.youngs:
                    ys = True
                try:
                    year = int(row['year'])
                    if year in years_sps.keys():
                        species_year = years_sps[year]
                        if sp_name in species_year:
                            species_year[sp_name].add(sq)
                        else:
                            species_year[sp_name] = set([sq])
                        
                        if ys:
                            if sp_name_new in species_year:
                                species_year[sp_name_new].add(sq)
                            else:
                                species_year[sp_name_new] = set([sq])
                            
                    else:
                        if year >= START_YEAR and year <= END_YEAR:
                            years_sps[year] = dict({sp_name : set([sq])})
                            if ys:
                                years_sps[year] = dict({sp_name_new : set([sq])})
                                
                except:
                    continue
                
        return years_sps
    
    def read_interactions(self):
        """
        This method reads the interactions between the species from the corresponding source file.
        """
        inters = csv.DictReader(open('../input/xarxa3.csv', 'rb'), delimiter=';', quotechar='"')
        
        self.interactions = dict()
        self.predators = set()
        self.preys = set()
        self.youngs = set()
        for row in inters:
            predator = row["Depredador"]
            prey = row["Presa"]
            group = int(row["EspecieGrup"])
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
            
            if SP_SUFFIX in predator:
                self.youngs.add(predator)
            if SP_SUFFIX in prey:
                self.youngs.add(prey)
        
        return
    
    ##seasonality - we obtain the group of species that occurr together based on their seasonality
    def read_seasonality(self):
        """
        Returns a dictionary with the name of the species and their seasonality. The seasonality takes a string value from:
        'ss' (spring-summer), 'aw' (autumn-winter), 'ay' (all-year).
        """
        seasonality = csv.DictReader(open('../input/SeasonalityComplete.csv', 'rb'), delimiter=';', quotechar='"')
        
        seasons = dict()
        for row in seasonality:
            seasons[row["Species"]] = row["Seasonality"]
            
        if self.interactions == None:
            self.read_interactions()
            
        for sp in seasons.keys():
            sp_new = sp + ' ' + SP_SUFFIX
            if sp_new in self.youngs:
                seasons[sp_new] = 'ss'
        
        return seasons
    
    def read_pyrenees_quads(self):
        """
        This method reads exclusively the grids or locations that belong to the study to be performed.
        It should be a subset of the total locations for the whole set of species.
        """
        quads = csv.DictReader(open('../input/quadr_comarc.csv', 'rb'), delimiter=',', quotechar='"')
        quad_set = set()
        
        for row in quads:
            quad_set.add(row['ATRIBUT'])
        
        return quad_set
    
    def read_habitats(self):
        """
        This method obtains the sets of habitats for each species from the corresponding file.
        """
        habs = csv.DictReader(open('../input/EspeciesHabitatBONA.csv', 'rb'), delimiter=';', quotechar='"')
        habitats = dict()
        for row in habs:
            species = row['Species']
            habitat_code = row['HabitatCode']
            if habitats.has_key(species):
                habitats[species].add(habitat_code)
            else:
                habitats[species] = set([habitat_code])
                
        if self.interactions == None:
            self.read_interactions()
            
        if len(self.youngs) > 0:
            for sp in habitats.keys():
                new_sp = sp + ' ' + SP_SUFFIX
                if new_sp in self.youngs:
                    habitats[new_sp] = habitats[sp].copy()
            
        return habitats
    
    def read_habitats_min_max_elevations(self):
        """
        This method obtains the minimum and maximum altitudes for each habitat over the entire set of the study grids.
        It also obtains the overall min and max altitudes.
        """
        altitudes = csv.DictReader(open('../input/CaractQuadricules.csv', 'rb'), delimiter=';', quotechar='"')
        hab_elevations = dict()
        
        min_alt = None
        max_alt = None
        
        for row in altitudes:
            tmp_min = float(row['minAlt'])
            tmp_max = float(row['maxAlt'])
            
            if tmp_min != 0.0 and (min_alt == None or tmp_min < min_alt):
                min_alt = tmp_min
            if max_alt == None or tmp_max > max_alt:
                max_alt = tmp_max
            
            hab_code = row['CodiHabitat']
            if hab_elevations.has_key(hab_code):
                if hab_elevations[hab_code][0] == 0.0 or (tmp_min != 0.0 and tmp_min < hab_elevations[hab_code][0]):
                    hab_elevations[hab_code][0] = tmp_min
                if tmp_max > hab_elevations[hab_code][1]:
                    hab_elevations[hab_code][1] = tmp_max
            else:
                hab_elevations[hab_code] = [tmp_min,tmp_max]
            
        return hab_elevations, min_alt, max_alt
    
    def read_habitat_clusters(self):
        """
        Retrieves the clusters of grids, specifying which grids belong to which cluster
        """
        cls = csv.DictReader(open('../input/Clusters.csv', 'rb'), delimiter=';', quotechar='"')
        clusters = dict()
        for row in cls:
            cluster = int(row['Cluster6'])
            quad = row['Quadricula']
            if clusters.has_key(cluster):
                clusters[cluster].add(quad)
            else:
                clusters[cluster] = set([quad])
        
        return clusters
    
    def read_species_mass(self):
        """
        Returns a dictionary with the species names as keys and two values contained in another recursive dictionary, the 'biomass'
        and the taxonomic 'group' the species belongs to. 
        """
        masses = csv.DictReader(open('../input/SpeciesMass.csv', 'rb'), delimiter='\t', quotechar='"')
        
        mass = dict()
        for row in masses:
            sp_name = row["SpeciesName"].strip() 
            mass[sp_name] = dict()
            biomass = row["Mass_g"].strip()
            group = row["SpeciesBigGroup"].strip()
            if biomass == '':
                mass[sp_name]['biomass'] = '0.0'
            else:
                mass[sp_name]['biomass'] = biomass
            if group == '':
                mass[sp_name]['group'] = 'N/A'
            else:
                mass[sp_name]['group'] = group 
        
        return mass
    
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
    
    def process_data_observations(self, years_sps, obs_yr):
        """
        Returns two arrays and a dictionary containing: the number of species, the number of observations, and the number of observations
        per species respectively. Each value in the array is related to a year in the *obs_yr* dict.
        """
        #here we store the number of species per year
        species_nos = []
        #this structure keeps record of the number of observations made per year
        obs_nos = []
        #here we keep record of the number of observations per year per species
        obs_nos_ps = dict.fromkeys(years_sps.keys())
        
        for k in obs_nos_ps.keys():
            obs_nos_ps[k] = []
        
        years = obs_yr.keys()
        
        for y in years:
            species_nos.append(len(obs_yr[y]))
            composition = obs_yr[y]
            obs = 0
            
            for k in composition.keys():
                obs += composition[k]
            
            obs_nos.append(obs)
             
            for k1 in obs_nos_ps.keys():
                if composition.has_key(k1):    
                    obs_nos_ps[k1].append(composition[k1])
                else:
                    obs_nos_ps[k1].append(0)
        
        return species_nos, obs_nos, obs_nos_ps
    
    def differences_between_datasets(self):
        """
        This method is used to calculate differences between species data sets coming from different classifications
        in order determine what species are missing from which data sources.
        """
        self.read_interactions()
         
        sps = set(self.read_species().keys())
        
        net_sps = set(self.predators) | set(self.preys)
        my_sps = sps & net_sps
        
        print len(my_sps)
        
        #sps_habitats = set(read_habitats().keys())
        
        sps_seasons = set(self.read_seasonality().keys())
        
        diff = sorted(my_sps - sps_seasons)
        
        print 'species without seasonality:'
        for i in diff:
            print i
        
        mass = set(self.read_species_mass().keys())
        
        diff = sorted(my_sps - mass)
        print '\n \n \nspecies without mass:'
        for i in diff:
            print i
    
        habitats = set(self.read_habitats().keys())
        diff = sorted(my_sps - habitats)

        print '\n \n \nspecies without habitat:'
        for i in diff:
            print i