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
import cosasi

######################
### SAVE/LOAD DATA ###
######################
#write a dictionary using the pickle function
#input: d = dictionary to write to a file, outFile = file to save the dictionary to
def writeDict(d, fileName):
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

#######################
### RUN EXPERIMENTS ###
#######################
#
def parallelTests(jobs, func, procs=1): #NOTE: this assumes each job is a tuple
  # print("started: parallelTests")
  pool = mp.Pool(processes=procs)
  # print("set pool, settting results")
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

#function to run experiments on graphs 
def runTests(args):
  G, graphID, random_seed, size, graph_param_1, graph_param_2 = args #unpack the arguments, this will depend on what the jobs look like

  #NOTE: this can change
  #pick a random source
  #simulate spread from this source
  #apply algorithms to estimate/predict patient 0
  #collect metrics
  #return results

  #spread parameters (TBD)
  # R0 = 1.5
  # recovery_rate = 0.4 #recovery rate
  # avg_deg = np.mean([G.degree(v) for v in G.nodes()])
  # infection_rate = (R0 * recovery_rate) / avg_deg #infection rate
  

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
    seed=random_seed
  )

  #simulating spread for t_max time-steps
  contagion.forward(t_max)

  pass

###############
### METRICS ###
###############
#

#function to generate random BA graphs with given params
def genGraphsBA(size, m, graph_id, graph_type):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.barabasi_albert_graph(size, m, random_seed) #generating random BA graph

  return_tuple = tuple((G, graph_id, graph_type, random_seed, size, m, None)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed,...)
  return return_tuple

#function to generate random ER graphs with given params
def genGraphsER(size, p, graph_id, graph_type):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.erdos_renyi_graph(size, p, random_seed) #generating random ER graph

  return_tuple = tuple((G, graph_id, graph_type, random_seed, size, p, None)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed,...)
  return return_tuple

#function to generate random WS graphs with given params
def genGraphsWS(size, k, p, graph_id, graph_type):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.watts_strogatz_graph(size, k, p, random_seed) #generating random WS graph

  return_tuple = tuple((G, graph_id, graph_type, random_seed, size, p, k)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed, m,...)
  return return_tuple

#function to generate graphs for simulation
def genGraphs(sizes):
  print("generating graphs for experiments") #test output
  
  experiment_graphs = [] #list of tuples to return
  #format of tuples: (graph, graph id, graph type, random seed, size, graph parameter 1, graph parameter 2)
  #graph types: 0 = BA tree, 1 = BA dense, 2 = ER, 3 = WS, 4 = RW_1, 5 = RW_2, etc.

  #variables to track number of random graphs for experiments 
  num_sizes = len(sizes) #total number of random graph sizes
  num_rand_g_types = 4 #total number of random graph types (BA tree, BA dense, ER, WS)
  num_rand_g_per_size = 100 #total number of random graphs per size (UPDATE VALUE FOR EXPERIMENTS)
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

  #Facebook wall posts network (largest connected component
  #reading in FB network data using pandas
  df = pd.read_csv("facebook-wosn-wall/out.facebook-wosn-wall", header=None, delimiter='\s+')
  df.drop(df.columns[[2, 3]], axis = 1, inplace=True) #dropping last two columns
  #generating network from dataframe (edge list)
  FB_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
  #induced subgraph of largest connected component
  FB_lcc = FB_graph.subgraph(max(nx.connected_components(FB_graph), key=len))
  #adding FB graph to experiment graph list
  experiment_graphs.append(tuple((FB_lcc, graph_id, graph_type, None, len(FB_lcc), None, None))) 
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
  #adding power grid graph to experiment graph list
  experiment_graphs.append(tuple((us_airport_graph, graph_id, graph_type, None, len(us_airport_graph), None, None)))
  #sanity check
  print("finished generating graph type:", graph_type)
  graph_id += 1 #incrementing graph_id
  graph_type += 1 #incrementing graph type

  #power grid network
  #reading in the power grid graph
  power_grid_graph = nx.read_gml("power/power.gml", label='id')
  #adding power grid graph to experiment graph list
  experiment_graphs.append(tuple((power_grid_graph, graph_id, graph_type, None, len(power_grid_graph), None, None)))
  #sanity check
  print("finished generating graph type:", graph_type)
  graph_id += 1 #incrementing graph_id
  graph_type += 1 #incrementing graph type

  #returns list of tuples with relevant graph info
  return experiment_graphs

#
if __name__=='__main__':
  outFile = 'patient_zero_data.dict'
  data = readDict(outFile) if glob.glob(outFile) else {}

  #parameters relevant to the experiments that you want to run
  sizes = [100, 500, 1000, 5000, 10000] #sizes of graphs
  # R0 = [1.5]
  t_max = 15 #max time step to sim spread
  time_steps = [t for t in range (0,t_max)] #range of time steps to sim spread
  infection_rates = [] #infection rates 
  repeats = 10 
  procs = 1

  #a list of random graphs generated using the parameters above along with other relevant information
  #for each random graph, it might be a good idea to use a specific random seed so that experiments are reproducible
  #a different seed for each graph can be generated with something like np.random.randint(0,999999)
  #these can be tuples that looks something like (graph id, random seed, ...)
  experiment_graphs = genGraphs(sizes) #calling function to generate graphs to experiment on

  #test output
  # for exp in experiment_graphs:
  #   print(exp)

  graph_dict = {"graphs": experiment_graphs} #dictionary of graphs to experiment on
  writeDict(graph_dict, "graph_dict.pickle") #saving graph dictionary to pickle file to avoid regeneration of graphs

  # graph_dict = readDict("graph_dict.pickle") #obtaining dicitonary in pickle file
  # experiment_graphs = graph_dict.values() #storing graphs to experiment on from pickle file
  
  #test output
  # for list in graph_dict_test.values():
  #   for val in list:
  #     print(val)

  # newData = parallelTests(experiment_graphs, runTests, procs=procs)
