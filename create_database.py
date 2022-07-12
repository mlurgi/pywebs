
import string
import csv
import MySQLdb


def read_habitats():
    habs = csv.DictReader(open('../input/EspeciesHabitatBONA.csv', 'rb'), delimiter=';', quotechar='"')
    habitats = set()
    for row in habs:
        habitats.add(row['HabitatCode'])
        
    return habitats

def read_species_habitats():
    habs = csv.DictReader(open('../input/EspeciesHabitatBONA.csv', 'rb'), delimiter=';', quotechar='"')
    habitats = dict()
    
    habs_db = get_habitats_from_database()
    
    for row in habs:
        species = row['Species'].strip()
        habitat_code = row['HabitatCode'].strip()
        
        if habitat_code == '':
            continue
        
        if habitats.has_key(species):
            habitats[species].add(habs_db[habitat_code])
        else:
            habitats[species] = set([habs_db[habitat_code]])
                
    return habitats

def read_locations():
    qs = csv.DictReader(open('../input/vertebrats_2.txt', 'rb'), delimiter=',', quotechar='"')
    quads = set()
    for row in qs:
        if row['Quad'] == '9999':
            continue
        quads.add(row['Quad'])
    
    return sorted(quads)

def read_taxons():
    taxons = dict()
    
    taxs = csv.DictReader(open('../input/TAXONS_VERTE.DEL', 'rb'), delimiter=',', quotechar='"')
    for row in taxs:
        sp_name = row["Genre"].strip()+" "+row["Species"].strip()
        taxons[row["Code"]] = sp_name
        
    print 'Number of species in TAXONS_VERTE =', len(taxons)
        
    return taxons

def read_species():
    #here we read the species data from biocat database
    
    avail_species = csv.DictReader(open('../input/vertebrats_2.txt', 'rb'), delimiter=',', quotechar='"')
    taxons = read_taxons()
    
    taxons_set = set(taxons.values())
    
    
    #additionally, we also include the information from the ICO (Catalan Institute of Ornithology) dataset in order to obtain
    #a better resolver network 
    avail_birds = csv.DictReader(open('../input/ICO_data_in_plots.csv', 'rb'), delimiter=';', quotechar='"')
    
    
    species_set = set()
    #we read the data from the dataset and store it in memory
    for row in avail_species:    
        if not taxons.has_key(row["CodiSps"]):
            continue
        sp_name = taxons[row["CodiSps"]]
        
        if string.find(sp_name, '.') != -1 or string.find(sp_name, '-') != -1 or string.find(sp_name, '/') != -1:
            continue
        species_set.add(sp_name)
    
    ico_set = set()
    #to this dataset we add more information for the bird species (ICO's dataset)
    for row in avail_birds:
        sp_name = row["species"]
        
        if string.find(sp_name, '.') != -1 or string.find(sp_name, '-') != -1 or string.find(sp_name, '/') != -1:
            continue
        species_set.add(sp_name)
    
        ico_set.add(sp_name)
        
    species_dict = dict()
    mass_group = read_species_mass()
    seasons = read_seasonality()
    
    for sp in species_set:
        group = None
        mass = None
        if mass_group.has_key(sp):
            group = mass_group[sp]['group']
            mass = mass_group[sp]['biomass']
        
        season = None
        if seasons.has_key(sp):
            season = seasons[sp]
            
        species_dict[sp] = dict()
        species_dict[sp]['group'] = group
        species_dict[sp]['biomass'] = mass
        species_dict[sp]['season'] = season
        
#    species_sort = sorted(species_dict.keys())
#    for i in species_sort:
#        print i,' ',species_dict[i]['group'], ' ', species_dict[i]['biomass'],' ', species_dict[i]['season']
    
#    diff = ico_set - taxons_set
#    print '\n\n\n'
#    print 'The species that are in ico and not in biocat are: '
#    
#    for i in diff:
#        print i
    
    
    return species_dict

def read_observations():
    observations = []
    #here we read the species data from biocat database
    avail_species = csv.DictReader(open('../input/vertebrats_2.txt', 'rb'), delimiter=',', quotechar='"')
    taxons = read_taxons()
    
    taxons_set = set(taxons.values())
    #additionally, we also include the information from the ICO (Catalan Institute of Ornithology) dataset
    avail_birds = csv.DictReader(open('../input/ICO_data_in_plots.csv', 'rb'), delimiter=';', quotechar='"')
    
    db_grids = get_grids_from_database()
    db_species = get_species_index_from_database()
    
    #we read the data from the dataset and store it in memory
    for row in avail_species:
        if not taxons.has_key(row["CodiSps"]):
            continue
        sp_name = taxons[row["CodiSps"]]
        
        if string.find(sp_name, '.') != -1 or string.find(sp_name, '-') != -1 or string.find(sp_name, '/') != -1:
            continue
        
        if row['Quad'] == '9999':
            continue
        
        try:
            year = int(row['Any'])
        except:
            continue
            
        if year < 1000:
            continue
        
        species_info = dict()
        species_info['species'] = sp_name
        species_info['species_id'] = db_species[sp_name]
        species_info['year'] = year
        species_info['grid_id'] = db_grids[row['Quad']]
        species_info['source_id'] = 2
        observations.append(species_info)
      
    #to this dataset we add more information for the bird species (ICO's dataset)
    for row in avail_birds:
        #if row["plot"] in squares:
        sp_name = row["species"]
        
        if string.find(sp_name, '.') != -1 or string.find(sp_name, '-') != -1 or string.find(sp_name, '/') != -1:
            continue
                
        try:
            year = int(row['year'])
        except:
            continue
            
        species_info = dict()
        species_info['species'] = sp_name
        species_info['species_id'] = db_species[sp_name]
        species_info['year'] = year
        species_info['grid_id'] = db_grids[row['plot']]
        species_info['source_id'] = 1
        observations.append(species_info)
      
    return observations


def read_pyrenees_quads():
    quads = csv.DictReader(open('../input/quadr_comarc.csv', 'rb'), delimiter=',', quotechar='"')
    quad_set = set()
    
    db_grids = get_grids_from_database()
    
    for row in quads:
        quad_set.add(db_grids[row['ATRIBUT']])
    
    return quad_set

def read_species_mass():
    masses = csv.DictReader(open('../input/SpeciesMass.csv', 'rb'), delimiter='\t', quotechar='"')
    
    groups_db = get_groups_from_database()
    
    mass = dict()
    for row in masses:
        sp_name = row["SpeciesName"].strip() 
        mass[sp_name] = dict()
        biomass = row["Mass_g"].strip()
        group = row["SpeciesBigGroup"].strip()
        if biomass == '':
            mass[sp_name]['biomass'] = None
        else:
            mass[sp_name]['biomass'] = float(biomass)
        if group == '':
            mass[sp_name]['group'] = None
        else:
            mass[sp_name]['group'] = groups_db[group] 
    
    return mass

def read_seasonality():
    seasonality = csv.DictReader(open('../input/SeasonalityComplete.csv', 'rb'), delimiter=';', quotechar='"')
    
    sea_db = get_seasonalities_from_database()
    
    seasons = dict()
    for row in seasonality:
        seasons[row["Species"]] = sea_db[row["Seasonality"]]
    
    return seasons


def read_interactions():
    inters = csv.DictReader(open('../input/xarxa3.csv', 'rb'), delimiter=';', quotechar='"')
    interactions = dict()
    
    species_db = get_species_index_from_database()
    
    missing_species = set()
    
    for row in inters:
        predator = row["Depredador"].strip()
        prey = row["Presa"].strip()
        
        try:
            predator_id = species_db[predator]
        except:
            missing_species.add(predator)
            continue
        
        try:
            prey_id = species_db[prey]
        except:
            missing_species.add(prey)
            continue
        

        group = int(row["EspecieGrup"])
        if  predator_id in interactions.keys():
            inter_sp = interactions[predator_id]      
            if group in inter_sp.keys():
                if not prey_id in inter_sp[group]:
                    inter_sp[group].add(prey_id)
            else:
                inter_sp[group] = set([prey_id])
        else:
            interactions[predator_id] = dict({group:set([prey_id])})
        
    
    print 'species not in db:'
    for sp in sorted(missing_species):
        print sp
    
    print interactions
    
    return interactions
    

def read_clusters():
    cls = csv.DictReader(open('../input/Clusters.csv', 'rb'), delimiter=';', quotechar='"')
    clusters = dict()
    
    db_locs = get_locations_from_database()
    
    for row in cls:
        grid_id = db_locs[row['Quadricula'].strip()]
        clusters[grid_id] = dict()
        clusters[grid_id]['cluster_6'] = int(row['Cluster6']) 
        clusters[grid_id]['cluster_10'] = int(row['Cluster10'])
        clusters[grid_id]['cluster_2'] = int(row['Cluster2'])
                    
    return clusters


def get_groups_from_database():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor(MySQLdb.cursors.DictCursor)
    cur.execute("SELECT * FROM species_group")

    groups = cur.fetchall()
    cur.close()
    db.close()
    
    dict_groups = dict()
    for row in groups:
        dict_groups[row['group_name']] = row['id']
    
    return dict_groups

def get_seasonalities_from_database():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor(MySQLdb.cursors.DictCursor)
    cur.execute("SELECT * FROM seasonality")

    seasons = cur.fetchall()
    cur.close()
    db.close()

    dict_seasons = dict()
    for row in seasons:
        dict_seasons[row['seasonality_name']] = row['id']
        
    return dict_seasons

def get_grids_from_database():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor(MySQLdb.cursors.DictCursor)
    cur.execute("SELECT * FROM locations")

    grids = cur.fetchall()
    cur.close()
    db.close()
    
    dict_grids = dict()
    for row in grids:
        dict_grids[row['location_name']] = row['id']
    
    return dict_grids

def get_habitats_from_database():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor(MySQLdb.cursors.DictCursor)
    cur.execute("SELECT * FROM habitats")
    
    habitats = cur.fetchall()
    cur.close()
    db.close()
    
    dict_habitats = dict()
    for row in habitats:
        dict_habitats[row['habitat_name']] = row['id']
    
    return dict_habitats

def get_species_index_from_database():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor(MySQLdb.cursors.DictCursor)
    cur.execute("SELECT * FROM species")
    
    species = cur.fetchall()
    cur.close()
    db.close()
    
    dict_sps_idx = dict()
    for row in species:
        dict_sps_idx[row['species_name']] = row['id']
    
    return dict_sps_idx

def get_locations_from_database():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor(MySQLdb.cursors.DictCursor)
    cur.execute("SELECT * FROM locations")
    
    db_locations = cur.fetchall()
    cur.close()
    db.close()
    
    dict_locations = dict()
    for row in db_locations:
        dict_locations[row['location_name']] = row['id']
    
    return dict_locations


def insert_species_data(species_data):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()
    
    species_sort = sorted(species_data.keys())
    for i in species_sort:
        group = species_data[i]['group'] 
        if group == None:
            group = 'NULL'
        else:
            group = str(group)
        mass = species_data[i]['biomass']
        if mass == None:
            mass = 'NULL'
        else:
            mass = str(mass)
        season = species_data[i]['season']
        if season == None:
            season = 'NULL'
        else:
            season = str(season)
    
        query = "INSERT INTO species (species_name, species_mass, seasonality_id, group_id) VALUES (\""  + i + "\", "+ mass +", "+ season +", "+ group +")"
        
        #print query
        cur.execute(query)

    cur.close()
    db.commit()
    db.close()


def read_grid_features():
    feats = csv.DictReader(open('../input/CaractQuadricules.csv', 'rb'), delimiter=';', quotechar='"')
    
    features = []
    
    grids_db = get_grids_from_database()
    habs_db = get_habitats_from_database()
    
    for row in feats:
        row['quadricula'] = grids_db[row['quadricula']]
        row['CodiHabitat'] = habs_db[row['CodiHabitat']]
        features.append(row)
    
    print features
    
    return features


def insert_grid_features_data(grids_data):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()
    
    for row in grids_data:
        query = "INSERT INTO locations_features (location_id, habitat_id, area, min_altitude, max_altitude) VALUES ("  + str(row['quadricula']) + ", "+ str(row['CodiHabitat']) +", "+ row['superficie'] +", "+ row['minAlt'] +", "+ row['maxAlt'] +")"
        #print query
        
        cur.execute(query)

    cur.close()
    db.commit()
    db.close()
    
def insert_species_habitats_data(sps_habs_data):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()
    
    species = sorted(set(sps_habs_data.keys()))
    spid_db = get_species_index_from_database()
    
    suffixes = ['young', 'nests', 'eggs']
    
    for sp in species:
        habs_array = sps_habs_data[sp]
        
        #there are species in habitats that are not in the biocat or ico databases
        if not spid_db.has_key(sp):
            print 'Species ',sp,' is not in database'
            continue
        
        sp_id = spid_db[sp]
#        young_id = None
#        young_sp = sp+' young'
    
        species_children = set()
        for sf in suffixes:
            species_children.add(sp+' '+sf)
            
        
        
#        if spid_db.has_key(young_sp):
#            young_id = spid_db[young_sp]
        
        sp_children_ids = set()
        for sp_ch in species_children:
            if spid_db.has_key(sp_ch):
                sp_children_ids.add(spid_db[sp_ch])
        
        for h in habs_array:
            query = "INSERT INTO species_habitats (species_id, habitat_id) VALUES ("+ str(sp_id) + ", "+ str(h) +")"
            #print query
            cur.execute(query)
            for sp_ch_id in sp_children_ids:
                 query = "INSERT INTO species_habitats (species_id, habitat_id) VALUES ("+ str(sp_ch_id) + ", "+ str(h) +")"
                 cur.execute(query)
            
#            if young_id != None:
#                query = "INSERT INTO species_habitats (species_id, habitat_id) VALUES ("+ str(young_id) + ", "+ str(h) +")"
#                cur.execute(query)

    cur.close()
    db.commit()
    db.close()

def insert_observations_data(observations):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()
    for ob in observations:
        query = "INSERT INTO observations (species_id, year, location_id, source_id) VALUES ("+ str(ob['species_id']) + ", "+ str(ob['year']) +", "+ str(ob['grid_id']) +", "+ str(ob['source_id']) +")"
        #print query
        cur.execute(query)

    cur.close()
    db.commit()
    db.close()

def insert_locations():
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    
    cur = db.cursor()
    locs = read_locations()
    
    for loc in locs:
        if loc != '':
            cur.execute("INSERT INTO locations (location_name) VALUES (\""  + loc + "\")")

    cur.close()  
    db.commit()
    db.close()

def insert_interactions_data(interactions):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()

    for predator_id in sorted(interactions.keys()):
        predator_relations = interactions[predator_id]
        for group_id in sorted(predator_relations.keys()):
            preys = predator_relations[group_id]
            for prey_id in sorted(preys):
                query = "INSERT INTO interactions (predator_id, prey_id, type_id) VALUES ("+ str(predator_id) + ", "+ str(prey_id) +", "+ str(group_id) +")"
                cur.execute(query)

    cur.close()
    db.commit()
    db.close()

def insert_study_grids_data(grids):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()

    for grid_id in grids:
        query = "INSERT INTO study_locations (location_id) VALUES ("+ str(grid_id) +")"
        cur.execute(query)

    cur.close()
    db.commit()
    db.close()

def insert_clusters_data(clusters):
    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
                         user="bernat", # your username
                         passwd="bernat", # your password
                         db="piriweb") # name of the data base

    cur = db.cursor()

    for loc_id in clusters.keys():
        query = "INSERT INTO clusters (location_id, cluster_6, cluster_10, cluster_2) VALUES ("+ str(loc_id) +", "+ str(clusters[loc_id]['cluster_6']) +", "+ str(clusters[loc_id]['cluster_10']) +", "+ str(clusters[loc_id]['cluster_2']) +")"
        cur.execute(query)

    cur.close()
    db.commit()
    db.close()


if __name__ == '__main__':
#    db = MySQLdb.connect(host="158.109.62.88", # your host, usually localhost
#                         user="bernat", # your username
#                         passwd="bernat", # your password
#                         db="piriweb") # name of the data base
#
#    # you must create a Cursor object. It will let
#    #  you execute all the query you need
#    cur = db.cursor(MySQLdb.cursors.DictCursor)
#    #cur = MySQLdb.cursors.DictCursor(db) 
##
##    locs = read_locations()
##    
##    for loc in locs:
##        if loc != '':
##            print loc
##            cur.execute("INSERT INTO locations (location_name) VALUES (\""  + loc + "\")")
##
#
#    cur.execute("SELECT * FROM species_group")
#
#    # print all the first cell of all the rows
#    for row in cur.fetchall() :
#        print row['id'], row['group_name']
#
#    cur.close()
#    
##    db.commit()
#    

#    insert_locations()
#    insert_species_data(read_species())
#    
#    insert_grid_features_data(read_grid_features())
    
    insert_species_habitats_data(read_species_habitats())   

#    insert_observations_data(read_observations())

#    insert_interactions_data(read_interactions())

#    insert_study_grids_data(read_pyrenees_quads())

#    insert_clusters_data(read_clusters())

    print 'nothing'

