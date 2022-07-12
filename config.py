

#connection details for the database
DB_HOST = 'localhost' #'158.109.62.88'
DB_USER = 'miguel' #'bernat'
DB_PASSWD = 'miguel' #'bernat'
DB_NAME = 'piriweb'

DATA_AGGREGATION = True
DATA_AGGREGATION_INTERVAL = 1

#networks to be generated each INTERVAL number of years. If INTERVAL = 1, one network per year will be generated
INTERVAL = 5 #20
START_YEAR = 1990 #1960
END_YEAR = 2004

ELEVATION_SLOTS = 1

#whether to calculate the jaccard similarity matrix
JACCARD = False

#possible values for the different clustering associations: 2, 6, 10
CLUSTERING = 1

CALCULATE_SEASONS = ['all', 'ss', 'aw', 'ay'] #['all']

#this parameter states the measurement to be used to calculate the size of the nodes when the network is drawn
#possible values: 
#    1.- 'degree' = degree centrality
#    2.- 'biomass' = max biomass of the species
#    3.- 'clustering' = clustering coefficient of the node
NODE_DIAMETER = 'biomass'

#this parameter indicates the species that is going to be highlighted in the network
HIGHLIGHTED = 'Ursus arctos'

#using this parameter we can define a suffix for a different class for a given species
#which will preserve the original characteristics of the parent species but different 
#trophic relationships, which means that the interactions file (currently xarxa3) needs to
#specify different interactions for this class. Also the 'biomasses' file must be updated to include
#a different biomass for this species' class.
#An example of the suffix is 'young'.
SP_SUFFIX = 'young'

SPS_SUFFIXES = ['young', 'eggs', 'nests']

#nodes to be removed for the network analyses and calculations
REMOVED_NODES = ['ProdPri', 'Carronya', 'Invertebrates', 'Brossa', 'Fish']


#the following parameters are considered for the robustness to extinction experiments:
#######
#whether or not to calculate the robustness for the output
CALCULATE_ROBUSTNESS = False

#the delta value for the percentage of species to be deleted from the network
ROBUSTNESS_INTERVAL = 10

#maximum percentage of species to be removed
ROBUSTNESS_MAX_REMOVED = 60

#whether the extinction process is cumulative (first eliminate all the species to be removed and then calculate secondary extinction)
#or not (calculate secondary extinction dynamically once after every species removal)
#NOTE: only for static robustness assessment (no weighted networks) (values: True / False)
ROBUSTNESS_CUMULATIVE = False

#when considering weighted networks there is an index for dynamic adaptation of weights after any given removal
#which is given by the following formula: 
#             alpha_ij' = alpha_ij + (1/number_of_available_prey)*beta*gamma_ij
# where beta is a parameter that determines the amount of resource previously taken from the extinct prey that the predator
# is going to be able to extract from its remaining preys, and gamma_ij is the weight of the link that predator had with the
# extinct prey. The update of the link weights is performed for every predator that preys upon the a removed prey and for all the links
# that predator has with other prey
ROBUSTNESS_BETA = 0.5

#this threshold determines the percentage of the original sum of weights a predator must keep before going to extinction
#e.g.: if a given predator starts the robustness experiments with a total weight sum (the sum of all the weights of his in_links) of 10,
# then it will survive until reaching a weight sum of five - if after the extinction of N of its preys the sum of link weights goes below
# this number (5) is will get extinct (secondary extinction)
ROBUST_EXT_THRESHOLD = 0.5
#######

#removal experiments

REMOVAL = False
REMOVAL_BETA = 0.7
REMOVAL_EXT_THRESHOLD = 0.7
REMOVED_SPECIES = ['Aquila chrysaetos']

#network to be saved to a file and the name of the corresponding file
#NOTE: the networks are saved in graphml format
SAVE_NETWORK = '1999-5-2-ss'

NET_TO_READ = '1980-5-1-ss.graphml'

#this flag allows to specify whether to perform invasion experiments on the obtained networks and report the output of
#such experiments on the output generated file (values: True / False)
INVASION = False

#flag that indicates the simulation of various invasion scenarios based on the INTRODUCED_GROUPS and the 
#PARAMETERS_OFFSET for each of the parameters governing the invasion dynamics
INVASION_SCENARIOS = False

#the following parameters apply to the process of species invasion of a given network
#any of the following values:
#Bird', 'BirdPrey', 'Reptile', 'Mammal', 'MammalCarn', 'Fish', 'Amphibian', 'Invertebrate'
INTRODUCED_GROUP = 'BirdPrey'

#array of introduced groups for the simulation of various invasion scenarios
INTRODUCED_GROUPS = ['Bird', 'BirdPrey', 'Reptile', 'Mammal', 'MammalCarn', 'Amphibian']

#name of the introduced species
INTRODUCED_SPECIES = 'Invader'

#biomass of the introduced species. If None calculated as the mean of the biomasses of the species
#in the same group present in the network
BIOMASS_INTRODUCED = None   #None/0.0

#fraction of generalism to be applied to the invasive species
FRACTION_OF_GENERALISM = 0.1

#vulnerability fraction to be applied to the invasive species
FRACTION_OF_PREDATORS = 0.1

#extinction threshold (fraction of the original weight sum) for prey species on which the introduced species feeds
#and might affect the pressure posed over that species by all of its predators
INVASION_EXT_THRESHOLD = 0.1

#offset by which the three parameters above are going to be incremented for each iteration of the run when
#generating the various invasion scenarios
PARAMETERS_OFFSET = 0.2

#the following parameters are used for the configuration of the netcarto program, which is used to compute
#the compartments in the network based on the algorithms implemented by Guimera in the rgraph library.
#for more information about these parameters and what they are used for consult the rgraph documentation

#whether to calculate the networks' modularity using the rgraph library
RGRAPH_MOD = False

#random seed for initialisation purposes (must be an integer)
RGRAPH_SEED = 345347

#The recommended range for the iteration factor is [0.1, 1]
RGRAPH_ITERATION_F = 1.0

#Recommended values for the cooling factor are [0.990, 0.999]
RGRAPH_COOLING_F = 0.993

#number of randomizations (random networks with the same connectivity as the original) that are going
#to be performed in order to test whether the modularity obtained from the original network is significant
RGRAPH_RANDOMIZATIONS = 0


