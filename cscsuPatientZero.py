#Template Python code for running experiments for CSCSU conference paper on identifying patient zero

import numpy as np
import networkx as nx
from itertools import product
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

#function to run experiments 
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
def genGraphsBA(size, m, graph_id):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.barabasi_albert_graph(size, m, random_seed) #generating random BA graph

  return_tuple = tuple((G, graph_id, random_seed, size, m, None))

  #returns a tuple (Graph G, graph_id, random_seed,...)
  return return_tuple

#function to generate random ER graphs with given params
def genGraphsER(size, p, graph_id):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.erdos_renyi_graph(size, p, random_seed) #generating random ER graph

  return_tuple = tuple((G, graph_id, random_seed, size, p, None))

  #returns a tuple (Graph G, graph_id, random_seed,...)
  return return_tuple

#function to generate random WS graphs with given params
def genGraphsWS(size, k, p, graph_id):

  random_seed = np.random.randint(0,999999) #random seed to generate graph

  G = nx.watts_strogatz_graph(size, k, p, random_seed) #generating random WS graph

  #return k and p values in tuple?
  return_tuple = tuple((G, graph_id, random_seed, size, p, k)) #tuple to return for experiments list

  #returns a tuple (Graph G, graph_id, random_seed, m,...)
  return return_tuple

#function to generate graphs for simulation
def genGraphs(sizes):
  print("generating graphs for experiments")
  experiments = [] #list of tuples to return 
  #format of tuples: (graph, graph id, random seed, size, graph parameter 1, graph parameter 2)

  #variables to track number of random graphs for experiments 
  num_sizes = len(sizes) #total number of random graph sizes
  num_rand_g_types = 3 #total number of random graph types (random, real world, etc)
  num_rand_g_per_size = 1 #total number of random graphs per size (UPDATE VALUE FOR EXPERIMENTS)
  num_rand_g_per_type = num_sizes * num_rand_g_per_size #total number of graphs per randon graph type
  total_rand_graphs = num_rand_g_per_type * num_rand_g_types #total number of graphs to generate

  graph_id = 0 #stores graph_id which relates to graph type

  iter = 0 #iterator variable

  #while loop to generate random graphs
  while graph_id < total_rand_graphs:

    #loop to create graph of each size listed
    for size in sizes:

      for iter_2 in range(num_rand_g_per_size):

        #if-else statements to create graph type based on graph_id
        if (graph_id in range(0,num_rand_g_per_type)): #BA (tree)
          m = 1 #m-value for tree BA graph?
          
          #updating experiments list with new tuple
          experiments.append(genGraphsBA(size, m, graph_id))
        
        elif(graph_id in range(num_rand_g_per_type,2*num_rand_g_per_type)): #BA (dense)
          m = 3 #m-value for dense BA graph

          #updating experiments list with new tuple
          experiments.append(genGraphsBA(size, m, graph_id))

        elif(graph_id in range(2*num_rand_g_per_type,3*num_rand_g_per_type)): #ER
          p = (math.log(size) + 1) / size #setting p-value to be in connected regime
          
          #updating experiments list with new tuple
          experiments.append(genGraphsER(size, p, graph_id))

        graph_id += 1 #incrementing graph_id
        iter += 1
    
    

    #cutoff point for generating random graphs (TBD)
    # if (graph_id >= 3*num_sizes):
    #   break

  #generate real-world networks below?

  #returns list of tuples with relevant graph info
  return experiments

#
if __name__=='__main__':
  outFile = 'patient_zero_data.dict'
  data = readDict(outFile) if glob.glob(outFile) else {}

  #parameters relevant to the experiments that you want to run
  sizes = [100, 500, 1000, 5000, 10000] #sizes of graphs
  # R0 = [1.5]
  t_max = 15 #max time step to sim spread
  time_steps = [t for t in range (0,t_max)] #range of time steps to sim spread
  infection_rates = []
  repeats = 10 
  procs = 1

  #a list of random graphs generated using the parameters above along with other relevant information
  #for each random graph, it might be a good idea to use a specific random seed so that experiments are reproducible
  #a different seed for each graph can be generated with something like np.random.randint(0,999999)
  #these can be tuples that looks something like (graph id, random seed, ...)
  experiments = []

  experiments = genGraphs(sizes) #calling function to generate graphs of different sizes

  #test output
  for exp in experiments:
    print(exp)

  newData = parallelTests(experiments, runTests, procs=procs)
