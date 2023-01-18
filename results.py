import numpy as np
import networkx as nx
import multiprocessing as mp
import cosasi
import copy
import rumor_centrality_code
import pandas as pd
import pickle
import matplotlib.pyplot as plt

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

def connectedProb(n):
  return (np.log(n) + 1) / n

def avgDistanceRank(experiment_data):
  #storing graph types/sizes experimented on
  graph_types = set()
  graph_sizes = set()
  beta_values = set()
  for id in experiment_data: 
    graph_types.add(experiment_data[id]['graph info']['type'])
    graph_sizes.add(experiment_data[id]['graph info']['params']['n'])
    beta_values.add(experiment_data[id]['spread info']['beta'])

  #storing maximum values
  time_step = max(experiment_data[0]['results'].keys())
  graph_size = max(graph_sizes)
  beta = min(beta_values)
  
  graph_types = sorted(list(graph_types))
  methods = list(experiment_data[0]['results'][time_step].keys())

  #keys = distance, values = frequency
  distance_freq = dict()
  rank_freq = dict()

  ap_diameter = 0
  if len(graph_types) == 1 and graph_types[0] == 'ap':
    ap_diameter = 8 #diameter of airport network
  else:
    diameter_dict = readDict("diams.pickle") #graph diameters

  diameter_list = []

  for graph_type in graph_types:
    print(graph_type)
    for method in methods:
      # print(type, size, beta, time, method)
      for id in experiment_data:
        if (experiment_data[id]['graph info']['type'] == graph_type and
            experiment_data[id]['graph info']['params']['n'] == graph_size and
            experiment_data[id]['spread info']['beta'] == beta):
          
          #store distance data
          distance = experiment_data[id]['results'][time_step][method]['distance']
          distance_freq[distance] = distance_freq.setdefault(distance, 0) + 1

          #store rank data
          true_source = experiment_data[id]['spread info']['source'] #true source of spread
          cosasi_result = experiment_data[id]['results'][time_step][method]['result']
          rank = cosasi_result.evaluate_solution_rank(true_source) #rank of true source
          rank_freq[rank] = rank_freq.setdefault(rank, 0) + 1

          #storing current graph's diameter
          if (ap_diameter == 0):
            graph_seed = experiment_data[id]['graph info']['seed']
            diameter_list.append(diameter_dict[(graph_seed, graph_type, graph_size)])
          else:
            diameter_list.append(ap_diameter)

      # print(graph_type, size, beta, time, method, len(distance_freq))
      
      #sorting dictionary by keys
      distance_freq = dict(sorted(distance_freq.items()))
      rank_freq = dict(sorted(rank_freq.items()))

      avg_distance = np.mean(list(distance_freq.keys()))
      avg_rank = np.mean(list(rank_freq.keys()))
      avg_diameter = np.mean(diameter_list)
      avg_distance_diameter = avg_distance / avg_diameter
      print(method, ":", avg_rank / graph_size, avg_distance_diameter)

      # print(distance_freq, method, size, beta) #test
      
      distance_freq.clear()
      rank_freq.clear()
      diameter_list.clear()

def distanceFreqFigure(experiment_data):
  #storing graph types/sizes experimented on
  graph_types = set()
  graph_sizes = set()
  beta_values = set()
  for id in experiment_data: 
    graph_types.add(experiment_data[id]['graph info']['type'])
    graph_sizes.add(experiment_data[id]['graph info']['params']['n'])
    beta_values.add(experiment_data[id]['spread info']['beta'])

  graph_types = sorted(list(graph_types))
  graph_sizes = sorted(list(graph_sizes))
  beta_values = sorted(list(beta_values))
  time_steps = sorted(list(experiment_data[0]['results'].keys()))
  methods = list(experiment_data[0]['results'][time_steps[0]].keys())

  # graph_types = ['ba_1']
  # graph_sizes = [500]
  # beta_values = [0.1, 0.3, 0.5]

  # print(graph_types, graph_sizes, beta_values, time_steps)

  #keys = distance, values = frequency
  distance_freq = dict()

  for graph_type in graph_types:
    for size in graph_sizes:
      for beta in beta_values:
        for time in time_steps:
          for method in methods:
            # print(type, size, beta, time, method)
            for id in experiment_data:
              if (experiment_data[id]['graph info']['type'] == graph_type and
                  experiment_data[id]['graph info']['params']['n'] == size and
                  experiment_data[id]['spread info']['beta'] == beta):
                
                #store distance data
                distance = experiment_data[id]['results'][time][method]['distance']
                distance_freq[distance] = distance_freq.setdefault(distance, 0) + 1
            print(graph_type, size, beta, time, method, len(distance_freq))
            
            #sorting dictionary by keys
            distance_freq = dict(sorted(distance_freq.items()))
            #averaging values
            avg_value = sum(distance_freq.values())
            for distance in distance_freq:
              distance_freq[distance] = distance_freq[distance] / avg_value
            # print(distance_freq, method, size, beta)
            plt.figure()
            plt.bar(distance_freq.keys(), distance_freq.values())
            plt.xlabel("Distance from true source")
            plt.ylabel("Frequency (%)")
            plt.ylim(0,0.6)
            plt.xlim(-1,7)
            plt.xticks(list(distance_freq.keys()))
            # plt.show()
            plt.savefig("time_dependence/"+graph_type+"_"+str(size)+"_"+str(beta)+"_"+method+"_"+str(time)+".jpg")
            plt.close()
            distance_freq.clear()

def rankFreqFigure(experiment_data):
  #storing graph types/sizes experimented on
  graph_types = set()
  graph_sizes = set()
  beta_values = set()
  for id in experiment_data: 
    graph_types.add(experiment_data[id]['graph info']['type'])
    graph_sizes.add(experiment_data[id]['graph info']['params']['n'])
    beta_values.add(experiment_data[id]['spread info']['beta'])

  graph_types = sorted(list(graph_types))
  graph_sizes = sorted(list(graph_sizes))
  beta_values = sorted(list(beta_values))
  time_steps = sorted(list(experiment_data[0]['results'].keys()))
  methods = list(experiment_data[0]['results'][time_steps[0]].keys())
  
  # print(graph_types, graph_sizes, beta_values, time_steps)

  #keys = rank, values = frequency
  rank_freq = dict()

  for graph_type in graph_types:
    for size in graph_sizes:
      for beta in beta_values:
        for time in time_steps:
          for method in methods:
            # print(type, size, beta, time, method)
            for id in experiment_data:
              if (experiment_data[id]['graph info']['type'] == graph_type and
                  experiment_data[id]['graph info']['params']['n'] == size and
                  experiment_data[id]['spread info']['beta'] == beta):
                true_source = experiment_data[id]['spread info']['source'] #true source of spread
                cosasi_result = experiment_data[id]['results'][time][method]['result']
                rank = cosasi_result.evaluate_solution_rank(true_source) #rank of true source

                #store rank data
                rank_freq[rank] = rank_freq.setdefault(rank, 0) + 1

            # print(graph_type, size, beta, time, method, len(distance_freq))
            
            #sorting dictionary by keys
            rank_freq = dict(sorted(rank_freq.items()))
            #averaging values
            avg_value = sum(rank_freq.values())
            for rank in rank_freq:
              rank_freq[rank] = rank_freq[rank] / avg_value
            # print(rank_freq, method, size, beta)
            # print(sum(rank_freq.values()))
            plt.figure()
            plt.bar(rank_freq.keys(), rank_freq.values())
            plt.xlabel("Rank")
            plt.ylabel("Frequency (%)")
            # plt.ylim(0,75)
            plt.xlim(0,max(rank_freq.keys()) + 1)
            # plt.xticks(list(rank_freq.keys()))
            # plt.show()
            plt.savefig("rank_figs/"+graph_type+"_"+str(size)+"_"+str(beta)+"_"+method+"_"+str(time)+".jpg")
            plt.close()
            rank_freq.clear()

def transmissionRank(experiment_data):
  #storing graph types/sizes experimented on
  graph_types = set()
  graph_sizes = set()
  beta_values = set()
  for id in experiment_data: 
    graph_types.add(experiment_data[id]['graph info']['type'])
    graph_sizes.add(experiment_data[id]['graph info']['params']['n'])
    beta_values.add(experiment_data[id]['spread info']['beta'])

  graph_types = sorted(list(graph_types))
  graph_sizes = sorted(list(graph_sizes))
  beta_values = sorted(list(beta_values))
  time_steps = sorted(list(experiment_data[0]['results'].keys()))
  methods = list(experiment_data[0]['results'][time_steps[0]].keys())
  
  # print(graph_types, graph_sizes, beta_values, time_steps)

  #keys = transmission rates, values = frequency of top-N accuracy
  transmission_rank = dict()

  top_N = 10

  for graph_type in graph_types:
    for size in graph_sizes:
      # top_N = int(0.1 * size)
      for beta in beta_values:
        for time in time_steps:
          for method in methods:
            accuracy_ct = 0
            # print(type, size, beta, time, method)
            for id in experiment_data:
              if (experiment_data[id]['graph info']['type'] == graph_type and
                  experiment_data[id]['graph info']['params']['n'] == size and
                  experiment_data[id]['spread info']['beta'] == beta):

                true_source = experiment_data[id]['spread info']['source'] #true source of spread
                cosasi_result = experiment_data[id]['results'][time][method]['result']
                rank = cosasi_result.evaluate_solution_rank(true_source) #rank of true source

                #store transmission data
                transmission = experiment_data[id]['spread info']['transmission'] #transmission rate
                if (rank <= top_N):
                  transmission_rank[transmission] = transmission_rank.setdefault(transmission, 0) + 1
                else:
                  transmission_rank[transmission] = transmission_rank.setdefault(transmission, 0)
                
                accuracy_ct += 1
            
            #sorting dictionary by keys
            transmission_rank = dict(sorted(transmission_rank.items()))
            #averaging values
            # avg_value = sum(distance_freq.values())
            for transmission in transmission_rank:
              transmission_rank[transmission] = transmission_rank[transmission] / accuracy_ct

            print(graph_type, size, beta, method, time, top_N, transmission_rank)
            #plotting results
            plt.figure()
            plt.bar(transmission_rank.keys(), transmission_rank.values())
            plt.xlabel(r"Transmission rates $\beta \langle k \rangle$")
            plt.ylabel("Top-"+str(top_N)+" accuracy")
            # plt.ylim(0,75)
            # plt.xlim(0,10)
            # plt.xticks(list(transmission_rank.keys()))
            # plt.show()
            plt.savefig("transmission_figs/"+graph_type+"_"+str(size)+"_"+str(beta)+"_"+method+"_"+str(time)+".jpg")
            plt.close()
            transmission_rank.clear()

#main
if __name__=='__main__':
  #storing experimental data
  # experiment_data = readDict("airport_graph_experiments.pickle")
  experiment_data = readDict("synthetic_experiments_2.pickle")

  # avgDistanceRank(experiment_data)
  # distanceFreqFigure(experiment_data)
  # rankFreqFigure(experiment_data)
  # transmissionRank(experiment_data)

  # print(experiment_data)  