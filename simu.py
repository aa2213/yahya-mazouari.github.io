import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from dataclasses import dataclass
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from swarm_sim import *
import networkx as nx
import random

# Chemin vers le fichier source des donnes
PATH = './Traces.csv'
# Nombre de temps a analyser
MAXTEMPS = 300
# Constantes de distance
MAX_RANGE = 60000
MID_RANGE = 40000
MIN_RANGE = 20000

###
### Manipulation du fichier source.
###

print("### Importation des données ### ")
df = pd.read_csv(PATH, header=0)

print("### Reformatage des données : ajout d'un index sur le temps ### ")
#On genere un tableau qui servira d'index de temps pour chaque echantillon
names = ['1']
i=2
while len(names)<10000 :
  names.append(''+ str(i))
  i=i+1

#On modifie le df pour ajouter en tete de colone notre index de temps
save = df.columns
df.columns = names
df.loc[-1] = save
df.index = df.index + 1 
df = df.sort_index()
print("### Reformatage des données : ajout d'un index sur les satelites ### ")
# On genere les noms des satelite 
satnames = ['satx1']
i=1
while len(satnames)<300 :
  if (i !=1):
    satnames.append('satx'+ str(i))
  satnames.append('saty'+ str(i))
  satnames.append('satz'+ str(i))
  i=i+1
# On applique notre index
df['coords'] = satnames
df = df.set_index('coords', drop=True)

# On transpose notre df pour avoir le temps en ligne et les coordonees de satelites en colonnes
dft = df.transpose()


###
### Mise en forme de la donnee en memoire (Creation des structures de donnees necessaires)
###
def GetNodes(dft, time):
  nodes = {}
  node_id = 1
  while node_id < 101:
    node = Node(
      node_id - 1,
      dft.loc[str(time)]['satx' + str(node_id)],
      dft.loc[str(time)]['saty' + str(node_id)],
      dft.loc[str(time)]['satz' + str(node_id)]
    )
    nodes[node_id - 1] = node
    node_id += 1
  return nodes

def GetWeightedMatrix(SwarmObj):
  matrix20 = SwarmObj.neighbor_matrix(MIN_RANGE)
  matrix40 = SwarmObj.neighbor_matrix(MID_RANGE)
  matrix60 = SwarmObj.neighbor_matrix(MAX_RANGE)

  x = 0
  while x < len(matrix20):
    y = 0
    while y < len(matrix20):
      if matrix20[x][y] == 0 and matrix40[x][y] == 1:
        matrix20[x][y] = 2
      if matrix20[x][y] == 0 and matrix40[x][y] == 0 and matrix60[x][y] == 1:
        matrix20[x][y] = 3
      y += 1
    x += 1
  return matrix20

def generate_simulation_data(dft, max_temps=300):
  Positions = {}
  Swarms = {}
  Matrixes = {}

  for temps in range(max_temps):
    Positions[temps] = GetNodes(dft, temps + 1)
    Swarms[temps] = Swarm(MAX_RANGE, list(Positions[temps].values()))
    Matrixes[temps] = GetWeightedMatrix(Swarms[temps])

  return Positions, Swarms, Matrixes

Positions, Swarms, Matrixes = generate_simulation_data(dft, max_temps=MAXTEMPS)
def plot_swarm_3d(swarm, time_idx=0):
  nodes = swarm.nodes
  nx_graph = swarm.swarm_to_nxgraph()
  edges = list(nx_graph.edges())

  x_nodes = [node.x for node in nodes]
  y_nodes = [node.y for node in nodes]
  z_nodes = [node.z for node in nodes]

  edge_x, edge_y, edge_z = [], [], []
  for e0, e1 in edges:
    n0 = swarm.get_node_by_id(e0)
    n1 = swarm.get_node_by_id(e1)
    edge_x.extend([n0.x, n1.x, None])
    edge_y.extend([n0.y, n1.y, None])
    edge_z.extend([n0.z, n1.z, None])

  trace_edges = go.Scatter3d(
    x=edge_x, y=edge_y, z=edge_z,
    mode='lines',
    line=dict(color='gray', width=1),
    hoverinfo='none',
    name='Liens'
  )

  trace_nodes = go.Scatter3d(
    x=x_nodes, y=y_nodes, z=z_nodes,
    mode='markers',
    marker=dict(size=4, color='blue', opacity=0.8),
    text=[f'Sat {n.id}' for n in nodes],
    hoverinfo='text',
    name='Satellites'
  )

  fig = go.Figure(data=[trace_edges, trace_nodes])
  fig.update_layout(
    title=f'Topologie 3D de l’essaim à t={time_idx}',
    showlegend=False,
    template="plotly_white",
    scene=dict(
      xaxis_title='X',
      yaxis_title='Y',
      zaxis_title='Z'
    )
  )
  fig.show()


plot_swarm_3d(Swarms[0], time_idx=0)
# On itère sur les elements de notre df pour en extraire des Node comme implantes dans swarm_sim
def GetNodes(time):
  nodes = {}
  id = 1
  while id < 101:
    node = Node(id - 1, dft.loc[str(time)]['satx'+str(id)], dft.loc[str(time)]['saty'+str(id)], dft.loc[str(time)]['satz'+str(id)])
    nodes[id - 1] = node
    id = id +1
  return nodes  

# Cette fonction permet, a partir du swarm, d'extraire une matrice d'adjaceance qui prend en compte le cout des differents liens.
def GetWeightedMatrix(Swarm):
  # On prend comme base le reseau avec une portee minimale puis on vient ajouter les chemins necessitant plus de portee en multipliant par un coefficient.
  matrix20 = Swarm.neighbor_matrix(MIN_RANGE)
  matrix40 = Swarm.neighbor_matrix(MID_RANGE)
  matrix60 = Swarm.neighbor_matrix(MAX_RANGE)
  x=0
  y=0 
  while x < len(matrix20) :
    y=0
    while y < len(matrix20):
      if (matrix20[x][y] == 0 and matrix40[x][y] == 1) :
        matrix20[x][y]=2
      if (matrix20[x][y] == 0 and matrix40[x][y] == 0 and matrix60[x][y] == 1):
        matrix20[x][y]=3 
      y=y+1
    x=x+1
  return matrix20

print("### Génération des tableaux de nodes en fonction du temps ### ")
def GetPositions(): 
  Positions = {}
  temps = 0
  while temps < MAXTEMPS:
    Positions[temps] = GetNodes(temps+1)
    temps = temps+1
  return Positions
Positions = GetPositions()

print("### Génération des swarms correspondants a chaque instants ###")
def InitSwarms(Positions):
  Swarms = {}
  temps = 0
  while temps < MAXTEMPS:
    Swarms[temps] = Swarm(MAX_RANGE, list(Positions[temps].values()))
    temps = temps+1
  return Swarms
Swarms = InitSwarms(Positions)

print("### Génération des matrices d'adjacence ###")
def GenerateWeightedMatrix(Swarms):
  Matrixes = {}
  temps = 0
  while temps < MAXTEMPS:
    Matrixes[temps] = GetWeightedMatrix(Swarms[temps])
    temps = temps+1
  return Matrixes
Matrixes=GenerateWeightedMatrix(Swarms)

###
### Analyse de la donnee pour en extraire des metriques
###

@dataclass
class Metric:
  MeanDegree: float
  MeanClusterCoef : float
  Connexity : float #Nous prenon un float et non un boolen pour plus facilement calculer des statistiques comme le taux de connexite de l'essain au fil du temps
  Efficiency: float


# On itere sur tout les noeud d'un swarm et on moyenne leur coefficient de clustering
def GetMeanClusterCoef(Swarm):
  MeanClusterCoef = 0
  ClusterCoeffs = Swarm.cluster_coef()
  for f in ClusterCoeffs:
    MeanClusterCoef = MeanClusterCoef +f
  MeanClusterCoef = MeanClusterCoef/len(ClusterCoeffs)
  return MeanClusterCoef

# On itere sur tout les noeud d'un swarm et on moyenne leur degre
def GetMeanDegree(Swarm):
  MeanDegree = 0
  # Ici des resultats incoherents nous ont pousses a passer par l'intermediaire de networkx pour calculer le degre.
  Degrees = Swarm.swarm_to_nxgraph().degree()
  for f in Degrees:
    MeanDegree = MeanDegree + f[1]
  MeanDegree = MeanDegree/len(Degrees)
  return MeanDegree

# Les swarm etant initialises a partir de la plus grande portee, 
# un graphe disposant de plusieurs composante connectees est forcement non connexe
def GetConnexity(Swarm):
  if (len(Swarm.connected_components()) <= 1):
    return 1
  else:
    return 0

def GetCentrality(Swarm, Matrix):
  listCentrality = []
  degrees = Swarm.swarm_to_nxgraph().degree()
  for f in range(0,len(Matrix)):
    centrality = 0
    node = Swarm.get_node_by_id(f)
    neighbors = node.getNeighborsId();
    for d in degrees:
      if d[0] == f or neighbors.__contains__(d[0]):
        centrality+=d[1]
    listCentrality.append(centrality)

  return listCentrality

''' Renvoie id du noeud qui a la meilleure centralité '''
def topCentrality(Swarm, Matrix, NumberOfNode):
  listCentrality = GetCentrality(Swarm, Matrix)

  mapCentrality = {}
  for id in range(0,len(listCentrality)):
    mapCentrality[id] = listCentrality[id]
    
  sortedCentrality = {k: v for k, v in sorted(mapCentrality.items(), key=lambda item: item[1], reverse = True)}
  ordreSuppression = list(sortedCentrality.keys())
  idNoeud = ordreSuppression[0:NumberOfNode]

  return idNoeud

# Cette fonction permet de traduire une matrice d'adjacence coefficientee en un graph networkx.
# Cette forme est necessaire car swarm_sim ne prend pas en compte les poids dans ses fonctions.
def MatrixtoNxGraph(Matrix):
  G = nx.Graph()
  x=0
  y=0
  while x < len(Matrix) :
    y=0
    while y < len(Matrix):
      if x>y:
        if Matrix[x][y] != 0:
          G.add_edge(x, y, weight=Matrix[x][y])
      y=y+1
    x=x+1
  return G

# Cette fonction permet de calculer l'efficacité du réseau.
# L'efficacité est définie comme la moyenne des inverses des longueurs des plus courts chemins entre toutes les paires de nœuds.
def GetEfficiency(Swarm,Matrix):

  G = MatrixtoNxGraph(Matrix)
  N = len(Matrix)

  if N <= 1:
    return 0

  sum_theta = 0.0
  nodes = list(range(N))

  for source in nodes:
    # Calcul de toutes les distances depuis source en une seule passe
    try:
      distances = nx.single_source_dijkstra_path_length(G, source, weight='weight')
    except Exception:
      continue

    for target in nodes:
      if target <= source:
        continue
      if target in distances and distances[target] > 0:
        sum_theta += 1.0 / distances[target]

  number_of_pairs = N * (N - 1) / 2
  if number_of_pairs == 0:
    return 0

  return sum_theta / number_of_pairs


# Retourne l'importance des noeuds sous forme d'un dictionnaire avec l'id (id de 1 à 100 vu que l'on utilise la matrice)
def GetNumberOfShorthestPathByNode(Matrix):
  G = MatrixtoNxGraph(Matrix)
  #listConnectedComponents = list(nx.connected_components(G));
  listShortestPath = list(nx.all_pairs_dijkstra_path(G))
  #listShortestPath = list(Swarm.number_of_shortest_path_through_node()); 
  # Tableau du nombre de plus court chemin pour chaque satellite
  NumberOfSP = {}
  for i in listShortestPath:
    SourceNode = i[0]
    # Nombre de plus court chemin pour le satellite observé
    NumberOfSPForNode = 0
    #a activer si on ne veut pas calculer les shortest path dans les groupes non connectés
    #if(listConnectedComponents[0].__contains__(SourceNode)):
    for j in listShortestPath:
      TargetNode = j[0]
      # Graph non orienté, donc aller et retour = même chemin 
      # (ne se compte qu'une fois)
      for k in j[1]:
          if(k < TargetNode and not(TargetNode == SourceNode and k == TargetNode) and j[1][k].__contains__(SourceNode)):
            NumberOfSPForNode+=1
    NumberOfSP[SourceNode] = NumberOfSPForNode
  return NumberOfSP

# Retourne l'id du noeud le plus important
def GetTopImportanceNoeud(Matrix, NumberOfNode):
  ImportanceByNode = GetNumberOfShorthestPathByNode(Matrix)
  # tri de l'importance des noeuds
  SortedImportance = dict(sorted(ImportanceByNode.items(), key=lambda item: item[1]))
  # [::-1] : on met dans l'ordre décroissant
  # puis on retourne le nombre de noeuds voulus
  return list(SortedImportance.keys())[::-1][0:NumberOfNode]

# Cette fonction permet de faire tourner toutes les fonctions d'analyse de metrique definie plus haut sur un temps.
def AnalyzeSingleGraph(Position, Swarm, Matrix):
  MeanDegree = GetMeanDegree(Swarm)
  MeanClusterCoef = GetMeanClusterCoef(Swarm)
  Connexity = GetConnexity(Swarm)
  Efficiency = GetEfficiency(Swarm, Matrix)
  return Metric(MeanDegree,MeanClusterCoef,Connexity,Efficiency)

# Cette fonction permet de faire tourner toutes les fonctions d'analyse de metrique definie plus haut sur les differents temps.
def AnalyzeGraph(Positions, Swarms, Matrixes):
  Metrics = {}
  temps = 0
  while temps < MAXTEMPS:
    Metrics[temps] = AnalyzeSingleGraph(Positions[temps], Swarms[temps], Matrixes[temps])
    temps= temps+1
  return Metrics
Stats = AnalyzeGraph(Positions, Swarms, Matrixes)

###
### Selection des elements d'interet au sein de la donnee 
###

def GetWorstCase(Stats):
  counter=0
  worstEffi=1
  worstEffiIndex=0
  for Stat in Stats:
    if Stats[Stat].Connexity == 1:
      if Stats[Stat].Efficiency < worstEffi:
        worstEffi = Stats[Stat].Efficiency
        worstEffiIndex=counter
    counter += 1
  return worstEffiIndex

def GetBestCase(Stats):
  counter=0
  BestEffi=0
  BestEffiIndex=0
  for Stat in Stats:
    if Stats[Stat].Connexity == 1:
      if Stats[Stat].Efficiency >= BestEffi:
        BestEffi = Stats[Stat].Efficiency
        BestEffiIndex=counter
    counter += 1
  return BestEffiIndex

###
### Creation de scenarios en l'absence de noeuds en fonction des differentes strategies
###
# Cette fonction permet de retourner un tableau de node dans lequel ne figure pas les node d'index fournis
def CreateAlterWithout(Position, IndexList):
  UpdatedPosition = []
  for node in range(len(Position)):
    if not node in IndexList:
      UpdatedPosition.append(Position[node])
  return UpdatedPosition

def StatFromPosition(Position):
  NewSwarm = Swarm(MAX_RANGE, Position)
  NewMatrix = GetWeightedMatrix(NewSwarm)
  return AnalyzeSingleGraph(Position, NewSwarm, NewMatrix)

WorstIndex = GetWorstCase(Stats)
BestIndex = GetBestCase(Stats)

@dataclass
class ResilienceMetric:
    Redundancy: float
    Disparity: float
    Criticality: int
    CriticalNodes: list
    FlowRobustness: float
    RoutingCost: float
    Efficiency: float

def GetNodeIdsFromMatrix(Matrix):
  return list(range(len(Matrix)))

MAX_PATHS_FOR_DISPARITY = 20

def MatrixToNxGraphWithNodes(Matrix):
  G = MatrixtoNxGraph(Matrix)
  nodes = GetNodeIdsFromMatrix(Matrix)
  G.add_nodes_from(nodes)
  return G

def GetShortestPathsBetweenAllPairs(Matrix, max_paths=MAX_PATHS_FOR_DISPARITY):
  G = MatrixToNxGraphWithNodes(Matrix)
  nodes = GetNodeIdsFromMatrix(Matrix)
  shortest_paths = {}

  for source_index in range(len(nodes)):
    source = nodes[source_index]
    predecessors, distances = nx.dijkstra_predecessor_and_distance(G, source, weight='weight')
    for target in nodes[source_index + 1:]:
      shortest_paths[(source, target)] = GetLimitedShortestPaths(source, target, predecessors, max_paths)

  return shortest_paths

def GetLimitedShortestPaths(source, target, predecessors, max_paths):
  if target not in predecessors:
    return []

  paths = []

  def BuildPath(current_node, path):
    if len(paths) >= max_paths:
      return

    if current_node == source:
      paths.append(list(reversed(path)))
      return

    for predecessor in predecessors[current_node]:
      BuildPath(predecessor, path + [predecessor])

  BuildPath(target, [target])
  return paths

def GetShortestPathCountsFromGraph(G):
  nodes = list(G.nodes())
  number_of_pairs = len(nodes) * (len(nodes) - 1) / 2
  if number_of_pairs == 0:
    return {}

  path_counts = {}

  for source_index in range(len(nodes)):
    source = nodes[source_index]
    predecessors, distances = nx.dijkstra_predecessor_and_distance(G, source, weight='weight')
    counts = {}
    counts[source] = 1

    sorted_nodes = sorted(distances.keys(), key=lambda node: distances[node])
    for node in sorted_nodes:
      if node == source:
        continue
      counts[node] = 0
      for predecessor in predecessors[node]:
        counts[node] += counts[predecessor]

    for target in nodes[source_index + 1:]:
      if target in counts:
        path_counts[(source, target)] = counts[target]
      else:
        path_counts[(source, target)] = 0

  return path_counts

def GetRedundancyFromPathCounts(path_counts):
  number_of_pairs = len(path_counts)
  if number_of_pairs == 0:
    return 0

  number_of_alternative_paths = 0
  for pair in path_counts:
    number_of_alternative_paths += path_counts[pair]

  return number_of_alternative_paths / number_of_pairs

def GetPairDisparity(paths):
  if len(paths) <= 1:
    return 0

  disparity_sum = 0

  for i in range(len(paths)):
    for j in range(i + 1, len(paths)):
      internal_nodes_i = set(paths[i][1:-1])
      internal_nodes_j = set(paths[j][1:-1])
      max_internal_nodes = max(len(internal_nodes_i), len(internal_nodes_j))

      if max_internal_nodes > 0:
        different_nodes = internal_nodes_i.symmetric_difference(internal_nodes_j)
        disparity_sum += len(different_nodes) / max_internal_nodes


  return min(1,disparity_sum / (2*len(paths)))

def GetDisparityFromShortestPaths(shortest_paths):
  number_of_pairs = len(shortest_paths)
  if number_of_pairs == 0:
    return 0

  disparity_sum = 0
  for pair in shortest_paths:
    disparity_sum += GetPairDisparity(shortest_paths[pair])

  return disparity_sum / number_of_pairs

def GetCriticalityFromBetweenness(betweenness, threshold=0.02):
  critical_nodes = []

  for node in betweenness:
    if betweenness[node] >= threshold:
      critical_nodes.append(node)

  return critical_nodes

def GetFlowRobustness(Matrix):
  G = MatrixToNxGraphWithNodes(Matrix)
  nodes = GetNodeIdsFromMatrix(Matrix)
  number_of_pairs = len(nodes) * (len(nodes) - 1) / 2

  if number_of_pairs == 0:
    return 0

  connected_pairs = 0
  components = list(nx.connected_components(G))

  for component in components:
    component_size = len(component)
    connected_pairs += component_size * (component_size - 1) / 2

  return connected_pairs / number_of_pairs

def GetRoutingCost(Matrix):
  G = MatrixToNxGraphWithNodes(Matrix)
  nodes = GetNodeIdsFromMatrix(Matrix)
  total_cost = 0
  for source_index in range(len(nodes)):
    source = nodes[source_index]
    distances = nx.single_source_dijkstra_path_length(G, source, weight='weight')
    for target in nodes[source_index + 1:]:
      if target in distances:
        total_cost += distances[target]

  return total_cost

def AnalyzeSingleGraphResilience(Matrix, threshold=0.02):
  G = MatrixToNxGraphWithNodes(Matrix)
  shortest_paths = GetShortestPathsBetweenAllPairs(Matrix)
  path_counts = GetShortestPathCountsFromGraph(G)
  betweenness = nx.betweenness_centrality(G, normalized=True, weight='weight')
  critical_nodes = GetCriticalityFromBetweenness(betweenness, threshold)
  return ResilienceMetric(
    GetRedundancyFromPathCounts(path_counts),
    GetDisparityFromShortestPaths(shortest_paths),
    len(critical_nodes),
    critical_nodes,
    GetFlowRobustness(Matrix),
    GetRoutingCost(Matrix),
    GetEfficiency(None, Matrix)
  )

def AnalyzeGraphResilience(Matrixes, threshold=0.02):
  ResilienceStats = {}
  temps = 0
  while temps < MAXTEMPS:
    ResilienceStats[temps] = AnalyzeSingleGraphResilience(Matrixes[temps], threshold)
    temps = temps + 1
  return ResilienceStats

def PlotResilienceCurves(ResilienceStats):
  times = list(ResilienceStats.keys())
  redundancy_values = []
  disparity_values = []
  criticality_values = []
  flow_robustness_values = []
  routing_cost_values = []
  efficiency_values = []
  for temps in times:
    redundancy_values.append(ResilienceStats[temps].Redundancy)
    disparity_values.append(ResilienceStats[temps].Disparity * 100)
    criticality_values.append(ResilienceStats[temps].Criticality)
    flow_robustness_values.append(ResilienceStats[temps].FlowRobustness* 100)
    routing_cost_values.append(ResilienceStats[temps].RoutingCost)
    efficiency_values.append(ResilienceStats[temps].Efficiency * 100)

  plt.figure(figsize=(10, 5))
  plt.plot(times, efficiency_values, color='tab:purple')
  plt.xlabel('Temps')
  plt.ylabel('Efficacité (%)')
  plt.title('Efficacité en fonction du temps')
  plt.grid(True)
  plt.tight_layout()
  plt.savefig('courbe_efficacite.png')
  plt.close()

  plt.figure(figsize=(10, 5))
  plt.plot(times, routing_cost_values, color='tab:orange')
  plt.xlabel('Temps')
  plt.ylabel('Coût de routage (#)')
  plt.title('Coût de routage en fonction du temps')
  plt.grid(True)
  plt.tight_layout()
  plt.savefig('courbe_cout_routage.png')
  plt.close()

  plt.figure(figsize=(10, 5))
  plt.plot(times, flow_robustness_values, color='tab:purple')
  plt.xlabel('Temps')
  plt.ylabel('Robustesse de flux (%)')
  plt.title('Robustesse de flux en fonction du temps')
  plt.grid(True)
  plt.tight_layout()
  plt.savefig('courbe_robustesse_flux.png')
  plt.close()

  plt.figure(figsize=(10, 5))
  plt.plot(times, redundancy_values, color='tab:blue')
  plt.xlabel('Temps')
  plt.ylabel('Redondance (#)')
  plt.title('Redondance en fonction du temps')
  plt.grid(True)
  plt.tight_layout()
  plt.savefig('courbe_redondance.png')
  plt.close()

  plt.figure(figsize=(10, 5))
  plt.plot(times, disparity_values, color='tab:green')
  plt.xlabel('Temps')
  plt.ylabel('Disparite (%)')
  plt.title('Disparite en fonction du temps')
  plt.grid(True)
  plt.tight_layout()
  plt.savefig('courbe_disparite.png')
  plt.close()

  plt.figure(figsize=(10, 5))
  plt.plot(times, criticality_values, color='tab:red')
  plt.xlabel('Temps')
  plt.ylabel('Criticite (#)')
  plt.title('Criticite en fonction du temps')
  plt.grid(True)
  plt.tight_layout()
  plt.savefig('courbe_criticite.png')
  plt.close()

def GetStrategyDeletionOrders(ReferenceIndex, max_deleted):
  rng = random.Random(42)
  node_count = len(Positions[ReferenceIndex])
  importance_order = GetTopImportanceNoeud(Matrixes[ReferenceIndex], max_deleted)

  return {
    'metrique de base': [],
    'importance': list(importance_order),
    'centralite': topCentrality(Swarms[ReferenceIndex], Matrixes[ReferenceIndex], max_deleted),
    'random': rng.sample(range(node_count), max_deleted)
  }

def ComputeStrategyEvolution(ReferenceIndex, max_deleted=10, threshold=0.02):
  deletion_counts = list(range(max_deleted + 1))
  deletion_orders = GetStrategyDeletionOrders(ReferenceIndex, max_deleted)
  strategy_results = {}

  base_metric = Stats[ReferenceIndex]
  base_resilience = AnalyzeSingleGraphResilience(Matrixes[ReferenceIndex], threshold)

  for strategy in deletion_orders:
    strategy_results[strategy] = {
      'NombreSupprimes': deletion_counts,
      'Redondance': [],
      'Disparite': [],
      'Criticite': [],
      'DegreMoyen': [],
      'Efficacite': [],
      'RobustesseFlux': []
    }

    for nb_deleted in deletion_counts:
      if strategy == 'metrique de base':
        metric = base_metric
        resilience = base_resilience
        flow_robustness = GetFlowRobustness(Matrixes[ReferenceIndex])
      else:
        deleted_nodes = deletion_orders[strategy][0:nb_deleted]
        altered_position = CreateAlterWithout(Positions[ReferenceIndex], deleted_nodes)
        altered_swarm = Swarm(MAX_RANGE, altered_position)
        altered_matrix = GetWeightedMatrix(altered_swarm)
        metric = StatFromPosition(altered_position)
        resilience = AnalyzeSingleGraphResilience(altered_matrix, threshold)
        flow_robustness = GetFlowRobustness(altered_matrix)

      strategy_results[strategy]['Redondance'].append(resilience.Redundancy)
      strategy_results[strategy]['Disparite'].append(resilience.Disparity * 100)
      strategy_results[strategy]['Criticite'].append(resilience.Criticality)
      strategy_results[strategy]['DegreMoyen'].append(metric.MeanDegree)
      strategy_results[strategy]['Efficacite'].append(metric.Efficiency)
      strategy_results[strategy]['RobustesseFlux'].append(flow_robustness * 100)

  return strategy_results

def PlotStrategyMetric(strategy_results, metric_name, ylabel, title, filename):
  plt.figure(figsize=(10, 5))

  for strategy in strategy_results:
    plt.plot(
      strategy_results[strategy]['NombreSupprimes'],
      strategy_results[strategy][metric_name],
      marker='o',
      label=strategy
    )

  plt.xlabel('Nombre de satellites supprimes')
  plt.ylabel(ylabel)
  plt.title(title)
  plt.grid(True)
  plt.legend()
  plt.tight_layout()
  plt.savefig(filename)
  plt.close()

def PlotStrategyEvolutionCurves(ReferenceIndex, max_deleted=10, threshold=0.02):
  strategy_results = ComputeStrategyEvolution(ReferenceIndex, max_deleted, threshold)

  PlotStrategyMetric(
    strategy_results,
    'Redondance',
    'Redondance (#)',
    'Evolution de la redondance selon la strategie',
    'courbe_strategies_redondance.png'
  )

  PlotStrategyMetric(
    strategy_results,
    'Criticite',
    'Criticite (#)',
    'Evolution de la criticite selon la strategie',
    'courbe_strategies_criticite.png'
  )

  PlotStrategyMetric(
    strategy_results,
    'Disparite',
    'Disparite (%)',
    'Evolution de la disparite selon la strategie',
    'courbe_strategies_disparite.png'
  )

  PlotStrategyMetric(
    strategy_results,
    'DegreMoyen',
    'Degre moyen (#)',
    'Evolution du degre moyen selon la strategie',
    'courbe_strategies_degre_moyen.png'
  )

  PlotStrategyMetric(
    strategy_results,
    'Efficacite',
    'Efficacite',
    "Evolution de l'efficacite selon la strategie",
    'courbe_strategies_efficacite.png'
  )

  PlotStrategyMetric(
    strategy_results,
    'RobustesseFlux',
    'Robustesse de flux (%)',
    'Evolution de la robustesse de flux selon la strategie',
    'courbe_strategies_robustesse_flux.png'
  )

def StatResilienceFromPosition(Position, threshold=0.02):
  NewSwarm = Swarm(MAX_RANGE, Position)
  NewMatrix = GetWeightedMatrix(NewSwarm)
  return AnalyzeSingleGraphResilience(NewMatrix, threshold)

ResilienceStats = AnalyzeGraphResilience(Matrixes)
PlotResilienceCurves(ResilienceStats)
PlotStrategyEvolutionCurves(WorstIndex, 10)

