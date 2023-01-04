#Template Python code for running experiments for CSCSU conference paper on identifying patient zero

import numpy as np
import networkx as nx
from itertools import product
import pandas as pd
import multiprocessing as mp
import time #IMPORTANT - if you choose to compare algorithms with respect to run time, you will need to know how this interacts with multiprocessing
import glob
import pickle
import math
import lzma
import cosasi

######################
### SAVE/LOAD DATA ###
######################
#write a dictionary using the pickle function
#input: d = dictionary to write to a file, outFile = file to save the dictionary to
def writeDict(d, fileName):
  with open(fileName, 'wb') as o:
    pickle.dump(d, o, protocol=pickle.HIGHEST_PROTOCOL)

def writeDict_lzma(d, fileName):
  with open(fileName, 'wb') as o:
    pickle.dump(d, o, protocol=pickle.HIGHEST_PROTOCOL)

#read a dictionary from a pickled file
#input: inFile = file where the dictionary is saved
#return: the dictionary saved in the given file
def readDict(fileName):
  d = {}
  with open(fileName, 'rb') as f:
    d = pickle.load(f)
  return d

def readDict_lzma(fileName):
  d = {}
  with lzma.open(fileName, 'rb') as f:
    d = pickle.load(f)
  return d

###
### METHODS FOR AGGREGATED CENTRALITIES ###
###

# Normalize the values of a given dictionary D
def normalize_dict(D):
  norm_dict = dict()

  for i in D.keys():
    if (sum(D.values()) == 0):
      norm_dict[i] = 0
    else:
      norm_dict[i] = float(D[i]) / sum(D.values())

  return norm_dict

# Function to obtain probabilities for specific nodes given centrality of all nodes
def obtainProbs(G, centrality):
  temp_dict = dict() # Dictionary to store centrality

  # Looping throguh nodes to update temp_dict
  for v in G.nodes():
    if (G.nodes[v]['state'] == 'R'):
      temp_dict[v] = centrality[v]

  # Normalizing to obtain probability dictionary
  prob_dict = normalize_dict(temp_dict)

  return prob_dict

# Given 2 dictionaries with probabilities, aggregate them by calculating their average
def aggregate_avg(probs_1, probs_2):
  return_dict = dict() # Dictionary to return

  # Looping through all keys in given dictionary    
  for v in probs_1.keys():
      return_dict[v] = np.average([probs_1[v], probs_2[v]])

  return return_dict

# Function to determine the source using various centrality measures of infected subgraph
def findSource_aggregated_centralities(I):
  prob_dict = dict() # Dictionary with probabilities (value) for each node (key)

  # Degree centrality
  deg_cent = nx.degree_centrality(I)
  prob_deg_cent = obtainProbs(I, deg_cent)
  
  # Harmonic centrality
  harmonic_cent = nx.harmonic_centrality(I)
  prob_harmonic_cent = normalize_dict(harmonic_cent)

  # Updating probabilities
  prob_dict = aggregate_avg(prob_deg_cent, prob_harmonic_cent)

  # Eigenvector centrality
  eigen_cent = nx.eigenvector_centrality(I, tol=1.0e-3)
  prob_eigen_cent = obtainProbs(I, eigen_cent)

  # Updating probabilities
  prob_dict = aggregate_avg(prob_dict, prob_eigen_cent)

  # Closeness centrality
  close_cent = nx.closeness_centrality(I)
  prob_close_cent = obtainProbs(I, close_cent)

  # Updating probabilities
  prob_dict = aggregate_avg(prob_dict, prob_close_cent)
  
  #sorting probabilities dictionary
  sorted_prob_dict = {k: v for k,v in sorted(prob_dict.items(), key=lambda item: item[1], reverse=True)}

  #list of nodes by rank (rank[0] = predicted source)
  rank = [x for x in sorted_prob_dict.keys()]

  return rank

#######################
### RUN EXPERIMENTS ###
#######################
#
def parallelTests(jobs, func, procs=1): #NOTE: this assumes each job is a tuple
  print("started: parallelTests")
  pool = mp.Pool(processes=procs)
  results = pool.map_async(func, jobs)
  # print("set results, getting results")
  results = results.get()
  # print("got results, closing and joining pool")
  pool.close()
  pool.join()

  data = {} #dictionary to store experiment results
  for experiment in results: #reorganize and store the results of the experiments in a dictionary
                             #this should include graph id, random seed, and information about the parameters used to generate the graph and the spread
    pass
  return data


def rumorCentrality():
  pass

def jordanCentrality():
  pass

def netsleuth():
  pass

def lisn():
  pass

#function to run experiments on graphs 
def runTests(args):
  # print(args) #test output

  #unpack the arguments, this will depend on what the jobs look like
  G, infection_rate, t_max, graphID, graph_type, random_seed, size, avg_deg, graph_param_1, graph_param_2 = args 

  transmission_rate = infection_rate * avg_deg #transmission rate of current experiment
  spread_seed = np.random.randint(0,999999) #random seed for spread

  #NOTE: this can change
  #pick a random source
  #simulate spread from this source
  contagion = simSpread(G, infection_rate, t_max) #simulating spread
  true_source_node = contagion.get_infected_indices(step=0)[0] #source of spread sim

  true_source_cosasi = contagion.get_source() #cosasi 
  I = contagion.get_infected_subgraph(step=t_max) #storing infected subgraph at time-step t_max

  #apply algorithms to estimate/predict patient 0
  #collect metrics
  #return results


  

#######################
### SIMULATE SPREAD ###
#######################
#
def simSpread(G, infection_rate, random_seed, t_max): #this will require some parameters, it is worth thinking carefully about how the spread will work

  #implement methods from cosasi package below

  #defining spread model using cosasi package
  contagion = cosasi.StaticNetworkContagion(
    G=G,
    model="si", #simulating SI-model
    infection_rate=infection_rate,
    number_infected=1, #single source node
    seed=random_seed
  )

  #simulating spread for t_max time-steps
  contagion.forward(t_max)

  #returning spread model
  return contagion 

###############
### METRICS ###
###############
#

#function to generate random BA graphs with given params
def genGraphsBA(size, m, graph_type):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.barabasi_albert_graph(size, m, random_seed) #generating random BA graph

  #calculating avg degree
  num_edges = len(G.edges())
  avg_deg = (2 * num_edges) / size

  return_tuple = tuple((G, graph_type, random_seed, size, avg_deg, m, None)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed,...)
  return return_tuple

#function to generate random ER graphs with given params
def genGraphsER(size, p, graph_type):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.erdos_renyi_graph(size, p, random_seed) #generating random ER graph

  #calculating avg degree
  num_edges = len(G.edges())
  avg_deg = (2 * num_edges) / size

  return_tuple = tuple((G, graph_type, random_seed, size, avg_deg, p, None)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed,...)
  return return_tuple

#function to generate random WS graphs with given params
def genGraphsWS(size, k, p, graph_type):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.watts_strogatz_graph(size, k, p, random_seed) #generating random WS graph

  #calculating avg degree
  num_edges = len(G.edges())
  avg_deg = (2 * num_edges) / size

  return_tuple = tuple((G, graph_type, random_seed, size, avg_deg, p, k)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed, m,...)
  return return_tuple

#function to generate graphs for experiments
def genGraphs(sizes):
  print("generating graphs for experiments") #test output
  
  experiment_graphs = [] #list of tuples to return
  #format of tuples: (graph, graph id, graph type, random seed, size, avg_deg, graph param 1, graph param 2)
  #graph types: 0 = BA tree, 1 = BA dense, 2 = ER, 3 = WS, 4 = RW_1, 5 = RW_2, etc.

  #variables to track number of random graphs for experiments 
  num_sizes = len(sizes) #total number of random graph sizes
  num_rand_g_types = 4 #total number of random graph types (BA tree, BA dense, ER, WS)
  num_rand_g_per_size = 100 #total number of random graphs per size
  num_rand_g_per_type = num_sizes * num_rand_g_per_size #total number of graphs per randon graph type
  total_rand_graphs = num_rand_g_per_type * num_rand_g_types #total number of graphs to generate

  graph_id = 0 #stores graph_id which relates to graph type

  graph_type = 0 #current graph type being generated

  #while loop to generate random graphs
  while graph_id < total_rand_graphs:

    #loop to create graph of each size listed
    for size in sizes:

      # loop to create (num_rand_g_per_size) graphs for each size
      for iter in range(num_rand_g_per_size):

        #if-else statements to create graph type based on graph_id
        if (graph_id in range(0,num_rand_g_per_type)): #BA (tree)
          m = 1 #m-value for BA tree graph
          
          #updating experiments list with new tuple
          experiment_graphs.append(genGraphsBA(size, m, graph_id, graph_type))
        
        elif(graph_id in range(num_rand_g_per_type,2*num_rand_g_per_type)): #BA (dense)
          m = 3 #m-value for BA dense graph

          #updating experiments list with new tuple
          experiment_graphs.append(genGraphsBA(size, m, graph_id, graph_type))

        elif(graph_id in range(2*num_rand_g_per_type,3*num_rand_g_per_type)): #ER
          p = (math.log(size) + 1) / size #setting p-value to be in connected regime
          
          #updating experiments list with new tuple
          experiment_graphs.append(genGraphsER(size, p, graph_id, graph_type))

        elif(graph_id in range(3*num_rand_g_per_type,4*num_rand_g_per_type)): #WS
          #parameter values for WS model
          p = 0.4
          k = 4

          #updating experiments list with new tuple
          experiment_graphs.append(genGraphsWS(size, k, p, graph_id, graph_type))

        graph_id += 1 #incrementing graph_id
    
    #sanity check
    print("finished generating graph type:", graph_type)
    graph_type += 1 #incrementing graph type

  #generating real-world networks (FB wall posts, US airport, power grid)

  #Facebook wall posts network (largest connected component)
  #reading in FB network data using pandas
  df = pd.read_csv("facebook-wosn-wall/out.facebook-wosn-wall", header=None, delimiter='\s+')
  df.drop(df.columns[[2, 3]], axis = 1, inplace=True) #dropping last two columns
  #generating network from dataframe (edge list)
  FB_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
  #induced subgraph of largest connected component
  FB_lcc = FB_graph.subgraph(max(nx.connected_components(FB_graph), key=len))
  #calculating avg degree
  num_edges = len(FB_lcc.edges())
  avg_deg = (2 * num_edges) / size
  #adding FB graph to experiment graph list
  experiment_graphs.append(tuple((FB_lcc, graph_id, graph_type, None, len(FB_lcc), avg_deg, None, None))) 
  #sanity check
  print("finished generating graph type:", graph_type)
  graph_id += 1 #incrementing graph_id
  graph_type += 1 #incrementing graph type

  #US airport network
  #reading in airport network data using pandas
  df = pd.read_csv("opsahl-usairport/out.opsahl-usairport", header=None, delimiter='\s+')
  df.drop(df.columns[[2]], axis = 1, inplace=True) #dropping last two columns
  #generating network from dataframe (edge list)
  us_airport_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
  #calculating avg degree
  num_edges = len(us_airport_graph.edges())
  avg_deg = (2 * num_edges) / size
  #adding power grid graph to experiment graph list
  experiment_graphs.append(tuple((us_airport_graph, graph_id, graph_type, None, len(us_airport_graph), avg_deg, None, None)))
  #sanity check
  print("finished generating graph type:", graph_type)
  graph_id += 1 #incrementing graph_id
  graph_type += 1 #incrementing graph type

  #power grid network
  #reading in the power grid graph
  power_grid_graph = nx.read_gml("power/power.gml", label='id')
  #calculating avg degree
  num_edges = len(power_grid_graph.edges())
  avg_deg = (2 * num_edges) / size
  #adding power grid graph to experiment graph list
  experiment_graphs.append(tuple((power_grid_graph, graph_id, graph_type, None, len(power_grid_graph), avg_deg, None, None)))
  #sanity check
  print("finished generating graph type:", graph_type)
  graph_id += 1 #incrementing graph_id
  graph_type += 1 #incrementing graph type

  #returns list of tuples with relevant graph info
  return experiment_graphs

#function to generate a graph for experiments given a size and type
def genGraph(size, type):
  #format of tuples: (graph, graph type, random seed, size, avg_deg, graph param 1, graph param 2)
  #graph types: 0 = BA tree, 1 = BA dense, 2 = ER, 3 = WS, 4 = RW_1, 5 = RW_2, 6 = RW_3, etc.
  graph_info = tuple() #tuple to return containing relevant graph informationx

  #if-else statements to create graph type based on graph_id
  if (type == 0): #BA (tree)
    m = 1 #m-value for BA tree graph
    
    #updating experiments list with new tuple
    graph_info = genGraphsBA(size, m, type)
  elif(type == 1): #BA (dense)
    m = 3 #m-value for BA dense graph

    #updating experiments list with new tuple
    graph_info = genGraphsBA(size, m, type)
  elif(type == 2): #ER
    p = (math.log(size) + 1) / size #setting p-value to be in connected regime
    
    #updating experiments list with new tuple
    graph_info = genGraphsER(size, p, type)
  elif(type == 3): #WS
    #parameter values for WS model
    p = 0.4
    k = 4

    #updating experiments list with new tuple
    graph_info = genGraphsWS(size, k, p, type)
  #generating real-world networks (FB wall posts, US airport, power grid)
  elif(type == 4): #Facebook wall posts network (largest connected component)
    #reading in FB network data using pandas
    df = pd.read_csv("facebook-wosn-wall/out.facebook-wosn-wall", header=None, delimiter='\s+')
    df.drop(df.columns[[2, 3]], axis = 1, inplace=True) #dropping last two columns
    #generating network from dataframe (edge list)
    FB_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
    #induced subgraph of largest connected component
    FB_lcc = FB_graph.subgraph(max(nx.connected_components(FB_graph), key=len))
    #calculating avg degree
    num_edges = len(FB_lcc.edges())
    avg_deg = (2 * num_edges) / len(FB_lcc)
    #adding FB graph to experiment graph list
    graph_info = tuple((FB_lcc, type, None, len(FB_lcc), avg_deg, None, None))
  elif(type == 5): #US airport network
    #reading in airport network data using pandas
    df = pd.read_csv("opsahl-usairport/out.opsahl-usairport", header=None, delimiter='\s+')
    df.drop(df.columns[[2]], axis = 1, inplace=True) #dropping last two columns
    #generating network from dataframe (edge list)
    us_airport_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
    #calculating avg degree
    num_edges = len(us_airport_graph.edges())
    avg_deg = (2 * num_edges) / len(us_airport_graph)
    #adding power grid graph to experiment graph list
    graph_info = tuple((us_airport_graph, type, None, len(us_airport_graph), avg_deg, None, None))
  elif(type == 6): #power grid network
    #reading in the power grid graph
    power_grid_graph = nx.read_gml("power/power.gml", label='id')
    #calculating avg degree
    num_edges = len(power_grid_graph.edges())
    avg_deg = (2 * num_edges) / len(power_grid_graph)
    #adding power grid graph to experiment graph list
    graph_info = tuple((power_grid_graph, type, None, len(power_grid_graph), avg_deg, None, None))

  #returns tuple with relevant graph info
  return graph_info

#function to return list of infection rates to simulate
def getInfectionRates(max_infection_rate, inf_rate_increments):

  multiplier = 1 #multiplier to set upper bound for range
  range_upper_bound = max_infection_rate #upper bound for range

  #updating multiplier
  while(int(multiplier * inf_rate_increments) == 0):
    multiplier *= 10

  #updating upper bound for range
  while(type(range_upper_bound) != int):
    range_upper_bound *= multiplier
    if (int(range_upper_bound) != 0):
        range_upper_bound = int(range_upper_bound)
        break

  range_upper_bound += 1 #incremnting upper bound

  #setting infection rates
  infection_rates = [(x * inf_rate_increments) for x in range(1, range_upper_bound)]

  # print(infection_rates) #test output

  return infection_rates

#function to generate list of tuples for experiments
def setupExperiments(time_steps, infection_rates, sizes):

  experiments_list = [] #list of tuples with parameters for each experiments

  experiments_per_params = 10 #number of experiments for each set of parameters

  experiment_id = 0 #unique id for each experiment set of parameters
  num_rand_graph_types = 4 #number of random graph types
  total_graph_types = 7 #total number of graph types

  #nested for loops to store sets of parameters for random graphs
  for curr_type in range(num_rand_graph_types): #random graph type
    for size in sizes: #random graph sizes
      for beta in infection_rates: #infection rate for sim
        for t_max in time_steps: #t_max for sim
          for iter in range(experiments_per_params): #number of experiments for each parameter set
            #generating tuple with graph info (graph, type, size, seed, avg_deg, param_1, param_2)
            graph_tuple = genGraph(size, curr_type) 

            #tuple with experiment parameters (exp ID, graph info, t_max, infection rate)
            experiment_params = (experiment_id, graph_tuple, t_max, beta)
            experiment_id += 1 #incrementing experiments id
            experiments_list.append(experiment_params) #updating list with experiment tuples
            print(experiment_params) #test

  #nested for loops to store sets of parameters for real-world graphs
  for curr_type in range(num_rand_graph_types, total_graph_types): #random graph type
    for beta in infection_rates: #infection rate for sim
      for t_max in time_steps: #t_max for sim
        for iter in range(experiments_per_params): #number of experiments for each parameter set
          #generating tuple with graph info (graph, type, size, seed, avg_deg, param_1, param_2)
          graph_tuple = genGraph(None, curr_type) 

          #tuple with experiment parameters (exp ID, graph info, t_max, infection rate)
          experiment_params = (experiment_id, graph_tuple, t_max, beta)
          experiment_id += 1 #incrementing experiments id
          experiments_list.append(experiment_params) #updating list with experiment tuples
          print(experiment_params) #test
  
  return experiments_list

#main
if __name__=='__main__':
  outFile = 'patient_zero_data.dict'
  data = readDict(outFile) if glob.glob(outFile) else {}

  #parameters relevant to the experiments that you want to run
  sizes = [100, 500, 1000, 5000, 10000] #sizes of graphs
  t_max = 10 #max time step to sim spread
  time_steps = [t for t in range (0,t_max)] #range of time steps to sim spread
  max_infection_rate = 0.9 #max infection rate (less than 1)
  inf_rate_increments = 0.1 #increments of infection rates from zero to max-rate
  infection_rates = getInfectionRates(max_infection_rate, inf_rate_increments) #infection rates
  
  repeats = 10
  procs = 1

  #list of tuples with parameters for each experiment
  # experiments_list = setupExperiments(time_steps, infection_rates, sizes) 

  # experiment_dict = {"graphs": experiments_list} #dictionary of graphs to experiment on
  # writeDict(experiment_dict, "experiment_dict.pickle") #saving graph dictionary to pickle file to avoid regeneration of graphs

  #test output
  # for exp in experiments_list:
  #   print(exp)

  #a list of random graphs generated using the parameters above along with other relevant information
  #for each random graph, it might be a good idea to use a specific random seed so that experiments are reproducible
  #a different seed for each graph can be generated with something like np.random.randint(0,999999)
  #these can be tuples that looks something like (graph id, random seed, ...)
  # experiment_graphs = genGraphs(sizes) #calling function to generate graphs to experiment on
  # print("finished generating graphs for experiments") #test output

  #test output
  # for exp in experiment_graphs:
  #   print(exp)

  # graph_dict = {"graphs": experiment_graphs} #dictionary of graphs to experiment on
  # writeDict(graph_dict, "graph_dict.pickle") #saving graph dictionary to pickle file to avoid regeneration of graphs

  # writeDict_lzma(graph_dict, "graph_dictionary.xz") #saving graph dictionary to pickle file to avoid regeneration of graphs

  # print("storing graph dictionary from stored pickle file") #test output
  # graph_dict = readDict("graph_dict.pickle") #obtaining graph dicitonary in pickle file
  # experiment_graphs = graph_dict["graphs"] #storing graphs to experiment on from graph dictionary
  # print("stored graph dictionary/list") #test output
  
  #test output
  # for graph in experiment_graphs:
  #   print(graph)


  # newData = parallelTests(experiments_list, runTests, procs=procs)
  pass