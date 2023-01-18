import numpy as np
import networkx as nx
import multiprocessing as mp
import cosasi
import copy
import rumor_centrality_code
import pandas as pd
import pickle

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

def parallelTests(jobs, func, procs=1): #NOTE: this assumes each job is a tuple
  print("started: parallelTests")
  pool = mp.Pool(processes=procs)
  results = pool.map_async(func, [(i,job) for i,job in enumerate(jobs)])
  results = results.get()
  # results = [func((i,j)) for i,j in enumerate(jobs)] #comment out to use multiprocessing
  pool.close()
  pool.join()

  data = {} #dictionary to store experiment results
  for (experiment_id, info) in results:
    data[experiment_id] = info
  return data

#modified cosasi source code to generate infected
def get_infected_indices(G, history, step=0):
  """Retrieves the indices of all vertices in the infected compartment at the provided step.
  Parameters
  ----------
  step : int
      Iteration step
  Returns
  -------
  list
  """
  # nodes = list(G)

  def status_to_delta(status):
      """Converts the history's status to a vector representing movement in
      (+1) and out (-1) of the infected compartment
      Parameters
      ----------
      status : dict
          status dictionary from history, e.g. self.history[step]["status"]
      """
      # delta = np.zeros(len(self.G))
      delta = set()
      for idx in status:
          s = status[idx]
          if s == 1:
              # node became infected this step
              # delta[idx] = 1
              delta.add(idx)
      return delta

  if step >= len(history):
      raise ValueError(
          "Invalid step. Continue the simulation to reach this step."
      )
  # infected = np.zeros(len(self.G))
  infected = set()
  for s in range(step + 1):
      # infected += status_to_delta(self.history[s]["status"])
      infected = infected.union(status_to_delta(history[s]["status"]))
  return infected

def get_infected_subgraph(G, history, step=0):
  """Returns the subgraph of the contact network whose vertices are marked infected.
  Parameters
  ----------
  step : int
      Iteration step
  Returns
  -------
  NetworkX Graph
  Notes
  -----
  This is only guaranteed to be connected in the SI model.
  """
  infected_indices = get_infected_indices(G, history, step=step)
  not_infected_indices = set(G.nodes) - set(infected_indices)
  H = G.copy()
  H.remove_nodes_from(not_infected_indices)
  return H

def genExperimentParameters(graph_models, sizes, betas, num_rand_g_per_type, num_rand_per_beta):
  graph_types = []
  for d in graph_models:
    for s in sizes:
      c = copy.deepcopy(d)
      c['params']['n'] = s
      graph_types.append(c)
  experiments = []
  for d in graph_types:
    for gSeed in np.random.choice(100000, size=num_rand_g_per_type, replace=False):
      for beta in betas:
        for sSeed in np.random.choice(100000, size=num_rand_per_beta, replace=False):
          g = copy.deepcopy(d)
          g['seed'] = int(gSeed)
          s = {'source':np.random.choice(g['params']['n']), 'beta':beta, 'seed':int(sSeed)}
          info = {'graph info': g, 'spread info': s}
          experiments.append(info)
  return experiments

#function to run experiments on graphs 
def runTests(args):
  experiment_id, info = args
  # print(experiment_id, info['graph info'])
  G = genGraph(info)
  info['spread info']['transmission'] = info['spread info']['beta'] * np.mean(G.degree()) #transmission rate of current experiment
  # characteristic_time = 1 / info['spread info']['transmission'] #time required to reach a fraction of 1/e infected (~36%)
  top_N_val = 3 #value to use for top-N accuracy

  #pick a random source
  #simulate spread from this source
  contagion = cosasi.StaticNetworkContagion(
    G=G,
    model="si", #simulating SI-model
    infection_rate=info['spread info']['beta'],
    number_infected=1, #single source node
    seed=info['spread info']['seed'] #spread seed
  )
  #simulating spread for t_ma`x time-steps
  t_max = 1000 #can be max time_step of interest
  history = contagion.forward(t_max, verbose=True) ###set to some high number? 1000?
  # print("finished simulating spread", experiment_id)

  info['spread info']['source'] = contagion.get_source()[0] #true source of spread

  #time points at which to observe the graph
  time_steps = [10, 20, 30]
  results = {}
  for time_step in time_steps:
    I = contagion.get_infected_subgraph(step=time_step) #storing infected subgraph at time-step t_max
    # I = get_infected_subgraph(G, history, time_step)
  
    #method results
    rumor_result = rumor_centrality_code.rumor_centrality(I,G)
    jordan_result = cosasi.single_source.jordan_centrality(I, G)
    netsleuth_result = cosasi.single_source.netsleuth(I, G)
    lisn_result = cosasi.single_source.lisn(I, G, t=time_step)

    #inferred source for each method 
    rumor_source = rumor_result.rank()[0]
    jordan_source = jordan_result.rank()[0]
    netsleuth_source = netsleuth_result.rank()[0]
    lisn_source = lisn_result.rank()[0]

    #scores for each method
    rumor_scores = rumor_result.data['scores']
    jordan_scores = jordan_result.data['scores']
    netsleuth_scores = netsleuth_result.data['scores']
    lisn_scores = lisn_result.data['scores']

    #dictionary for each method with relevant metrics (result object, top-1 accuracy, top-N accuracy, distance, rank of true source)
    rumor_dict = {'result': rumor_result,
                  'estimated source': rumor_source,
                  'top-1': True if rumor_source == info['spread info']['source'] else False,
                  'top-N': True if rumor_source in rumor_result.topn(top_N_val) else False,
                  'distance': nx.shortest_path_length(I, rumor_source, info['spread info']['source']),
                  'rank': rumor_result.evaluate_solution_rank(info['spread info']['source']),
                  'scores': rumor_scores} 
    jordan_dict = {'result': jordan_result, 
                   'estimated source': jordan_source,
                   'top-1': True if jordan_source == info['spread info']['source'] else False,
                   'top-N': True if jordan_source in jordan_result.topn(top_N_val) else False,
                   'distance': nx.shortest_path_length(I, jordan_source, info['spread info']['source']),
                   'rank': jordan_result.evaluate_solution_rank(info['spread info']['source']),
                   'scores': jordan_scores}
    netsleuth_dict = {'result': netsleuth_result, 
                      'estimated source': netsleuth_source,
                      'top-1': True if netsleuth_source == info['spread info']['source'] else False,
                      'top-N': True if netsleuth_source in netsleuth_result.topn(top_N_val) else False,
                      'distance': nx.shortest_path_length(I, netsleuth_source, info['spread info']['source']),
                      'rank': netsleuth_result.evaluate_solution_rank(info['spread info']['source']),
                      'scores': netsleuth_scores}
    lisn_dict = {'result': lisn_result, 
                 'estimated source': lisn_source,
                 'top-1': True if lisn_source == info['spread info']['source'] else False,
                 'top-N': True if lisn_source in lisn_result.topn(top_N_val) else False,
                 'distance': nx.shortest_path_length(I, lisn_source, info['spread info']['source']),
                 'rank': lisn_result.evaluate_solution_rank(info['spread info']['source']),
                 'scores': lisn_scores}
    
    results[time_step] = {'rumor': rumor_dict,
                          'jordan': jordan_dict,
                          'netsleuth': netsleuth_dict,
                          'lisn': lisn_dict}

    # results[time_step] = {'rumor': rumor_centrality_code.rumor_centrality(I,G),
    #                       'jordan': cosasi.single_source.jordan_centrality(I, G),
    #                       'netsleuth': cosasi.single_source.netsleuth(I, G),
    #                       'lisn': cosasi.single_source.lisn(I, G, t=time_step)} #do we need to specify t here? (documentation says optional, but throws error if not specified)
  #for now, results for each method could be, for example, probability that each node is the source... we can do post processing later
      
  #apply algorithms to estimate/predict patient 0
  #collect metrics
  #return results

  info['results'] = results
  # print(experiment_id, info['results'])
  print(experiment_id)
  return (experiment_id, info)

#function to read in real world graphs
def readGraph(type):
  G = nx.Graph() #graph to return

  if (type == 'fb'): #Facebook wall posts network (largest connected component)
    #reading in FB network data using pandas
    df = pd.read_csv("facebook-wosn-wall/out.facebook-wosn-wall", header=None, delimiter='\s+')
    df.drop(df.columns[[2, 3]], axis = 1, inplace=True) #dropping last two columns
    #generating network from dataframe (edge list)
    FB_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
    #induced subgraph of largest connected component
    fb_LCC = FB_graph.subgraph(max(nx.connected_components(FB_graph), key=len))
    G = nx.convert_node_labels_to_integers(fb_LCC) #updating node labels
  elif(type == 'ap'): #US airport network
    #reading in airport network data using pandas
    df = pd.read_csv("opsahl-usairport/out.opsahl-usairport", header=None, delimiter='\s+')
    df.drop(df.columns[[2]], axis = 1, inplace=True) #dropping last two columns
    #generating network using pandas dataframe (edge list)
    airport_graph = nx.from_pandas_edgelist(df, source=0, target=1, edge_attr=None, create_using=nx.Graph)
    airport_LCC = airport_graph.subgraph(max(nx.connected_components(airport_graph), key=len))
    G = nx.convert_node_labels_to_integers(airport_LCC) #updating node labels
  elif(type == 'pg'): #power grid network
    #reading in the power grid graph
    G = nx.read_gml("power/power.gml", label='id')
  
  return G

def genGraph(info):
  if info['graph info']['type'][:2] == 'ba': 
    G = nx.barabasi_albert_graph(info['graph info']['params']['n'],
                                 info['graph info']['params']['m'],
                                 seed=info['graph info']['seed'])
  elif info['graph info']['type'][:2] == 'er': 
    er_G = nx.erdos_renyi_graph(info['graph info']['params']['n'],
                                info['graph info']['params']['p'](info['graph info']['params']['n']),
                                seed=info['graph info']['seed'])
    er_LCC = er_G.subgraph(max(nx.connected_components(er_G), key=len)) #using largest connected component of er graph
    G = nx.convert_node_labels_to_integers(er_LCC) #updating node labels
  elif info['graph info']['type'][:2] == 'ws': 
    G = nx.watts_strogatz_graph(info['graph info']['params']['n'],
                                info['graph info']['params']['k'],
                                info['graph info']['params']['p'],
                                seed=info['graph info']['seed'])
  else: G = readGraph(info['graph info']['type']) #function to read world graph in
  #update 'n' parameter in info for real world graphs?
  return G

def connectedProb(n):
  return (np.log(n) + 1) / n

if __name__=='__main__':
  graph_models = [{'type':'ba_1', 'params':{'m':1}},
                  {'type':'ba_3', 'params':{'m':3}},
                  {'type':'er', 'params':{'p':connectedProb}},
                  {'type':'ws_', 'params':{'p':0.01, 'k':4}}]
                  # {'type':'fb', 'params':{}},
                  # {'type':'ap', 'params':{}},
                  # {'type':'pg', 'params':{}}]
  sizes = [100, 250, 500]#, 5000, 10000]
  betas = [0.1, 0.2, 0.3, 0.4, 0.5] #or should this be a function of e.g. <k> or size?
  num_rand_g_per_type = 17
  num_rand_per_beta = 4
  experiments = genExperimentParameters(graph_models, sizes, betas, num_rand_g_per_type, num_rand_per_beta)
  graph_count = 0

  print(len(experiments))

  # for e in experiments: 
  #   print(e)


  # procs = 4
  # procs = mp.cpu_count() - 1
  # data = parallelTests(experiments, runTests, procs=procs)
  # for e in data: print(e, data[e])

  # writeDict(data, "synthetic_experiments_2.pickle")
  # d = readDict("synthetic_graph_experiments.pickle")
  # for item in d.items():
  #   print(item)
