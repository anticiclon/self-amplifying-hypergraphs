# -*- coding: utf-8 -*-
"""
Created on Mon May 19 12:10:36 2025

@author: Trabajador
"""
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import graphviz_layout
import scenario_generator as generator
from matplotlib.patches import Circle
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import auxiliary_functions as aux

# =============================================================================
def replaceReactionNodesWithAdj2WithFlows(graph):
    """
    Given a directed graph, this function replaces nodes with adjacency 2 by a direct edge.
    A node with adjacency 2 is defined as having exactly one predecessor and one successor.

    Parameters:
    - graph (nx.DiGraph): A directed graph created with networkx.

    Returns:
    - nx.DiGraph: A new graph with nodes of adjacency 2 replaced by a direct edge.
    """
    # Make a copy of the graph to avoid modifying the original
    new_graph = graph.copy()
    
    # Find nodes with exactly one predecessor and one successor
    nodes_to_replace = []
    for node in new_graph.nodes:
        predecessors = list(new_graph.predecessors(node))
        successors = list(new_graph.successors(node))
        
        # Check if the node has exactly one predecessor and one successor
        if (new_graph.nodes[node]["type"] == 'reaction' and 
            len(predecessors) == 1 and len(successors) == 1):
            nodes_to_replace.append((predecessors[0], node, successors[0]))
    
    # Replace each identified node with a direct edge
    for pred, node, succ in nodes_to_replace:
        # Add a direct edge between predecessor and successor
        new_graph.add_edge(pred, succ, weight = 1, flow = graph.nodes[node]["flow"])
        # Remove the node with adjacency 2
        new_graph.remove_node(node)
    
    return new_graph
# =============================================================================



# =============================================================================
def createGraphFromStoichiometricMatrixWithFlowsSpecial(stoichiometric_matrix, dict_flows, species_names=None, reaction_names=None):
    n_species, n_reactions = stoichiometric_matrix.shape
    species = range(n_species)
    reactions = range(n_reactions)

    # Create a directed bipartite graph
    G = nx.DiGraph()
    
    # Default names for species and reactions if none provided
    if species_names is None:
        species_names = [f"S{i+1}" for i in species]
    if reaction_names is None:
        reaction_names = [f"R{i+1}" for i in reactions]
        
    # New Dict
    new_dict_flows = {}
    for key, value in dict_flows.items():
        new_dict_flows["R" + str(key + 1)] = value


    # Add species nodes
    for i in species_names:
        G.add_node(i, type='species')
    
    # Add reaction nodes
    for idx, j in enumerate(reaction_names):
        if j in set(new_dict_flows.keys()):
            G.add_node(j, type='reaction', flow = new_dict_flows[j])
        else:
            G.add_node(j, type='reaction')

    # Add edges based on the stoichiometric matrix
    for i in species:
        for j in reactions:
            stoich_coeff = stoichiometric_matrix[i, j]
            if stoich_coeff < 0:  # Reactant (consumed)
                if "flow" in G.nodes()[reaction_names[j]]:
                    G.add_edge(species_names[i], reaction_names[j], weight = abs(stoich_coeff), flow = G.nodes()[reaction_names[j]]["flow"])
                else:
                    G.add_edge(species_names[i], reaction_names[j], weight = abs(stoich_coeff))
            elif stoich_coeff > 0:  # Product (produced)
                if "flow" in G.nodes()[reaction_names[j]]:
                    G.add_edge(reaction_names[j], species_names[i], weight = stoich_coeff, flow = G.nodes()[reaction_names[j]]["flow"])
                else:
                    G.add_edge(reaction_names[j], species_names[i], weight = stoich_coeff)

    # nx.draw(G, with_labels = True, pos = graphviz_layout(G, prog="neato"))
    return G
# =============================================================================



# =============================================================================
def plotChemicalReactionNetworkWithFlows(graph, subgraph, a_sol, sizeSpecies, title="Chemical Reaction Network"):
    
    # Set up node positions using graphviz for better layout
    pos = graphviz_layout(graph, prog="neato")

    # Start the plot
    fig, ax = plt.subplots()
    fig.canvas.manager.full_screen_toggle()
        
    sub_species, sub_reactions, sub_food, sub_waste = subgraph["species"], subgraph["reaction"], subgraph["food"], subgraph["waste"]

    subset_auto = [
        (u, v) for u, v in graph.edges()
        if ((u in sub_species and v in sub_reactions) or 
                (u in sub_reactions and v in sub_species) or
                (u in sub_species and v in sub_species)) and "flow" in graph[u][v]
    ]
    
        
    # Step 2: Extract flow values and normalize
    flows = [graph[u][v]['flow'] for u, v in subset_auto]

    norm = mcolors.Normalize(vmin=min(flows), vmax=max(flows))
    # Step 3: Choose a colormap and map normalized values to colors
    cmap = cm.get_cmap('plasma')  # Use the 'hot' colormap
    # cmap = cm.get_cmap('magma')  # Use the 'hot' colormap
    edge_colors = [cmap(norm(flow)) for flow in flows]

    dictionary_flows_colors = dict(zip(flows, edge_colors))
    
    # Separate nodes by type and assign sizes
    species_nodes = [node for node in graph.nodes if graph.nodes[node]["type"] == "species"]
    reaction_nodes = [node for node in graph.nodes if graph.nodes[node]["type"] == "reaction"]
    sizeList = [sizeSpecies if graph.nodes[node]["type"] == "species" else 150 for node in graph.nodes]

    species_nodes_nodified = list(set(species_nodes) - set(['C1a']))    
    
    nx.draw_networkx_nodes(graph, {n: pos[n] for n in species_nodes},
                           nodelist=species_nodes_nodified,
                           node_shape='o', node_color='none', edgecolors='black',
                           node_size=sizeSpecies, label='Species')
    
    nx.draw_networkx_nodes(graph, {n: pos[n] for n in reaction_nodes},
                            nodelist=[i for i in reaction_nodes if "flow" not in graph.nodes[i]],
                            node_shape='s', node_color='none', 
                            edgecolors="black",
                            node_size=150, label='Reactions', linewidths=1)


    nx.draw_networkx_nodes(graph, {n: pos[n] for n in reaction_nodes},
                            nodelist=[i for i in reaction_nodes if "flow" in graph.nodes[i]],
                            node_shape='s', node_color='none', 
                            edgecolors=[dictionary_flows_colors[graph.nodes[i]["flow"]] for i in reaction_nodes if "flow" in graph.nodes[i]],
                            node_size=150, label='Reactions', linewidths=3)
    
    nx.draw_networkx_nodes(graph, {n: pos[n] for n in a_sol},
                            nodelist = a_sol,
                            node_shape = 'o', node_color='gold', 
                            edgecolors = "none",
                            node_size = sizeSpecies, linewidths=0, alpha = 0.2)

    # Draw nodes manually with Circle patches for dashed borders
    for node, (x, y) in pos.items():
        if node == 'C1a':
            # Add a Circle patch for each node
            circle = Circle(
                (x, y),  # Position
                radius = 8,  # Size of the node
                edgecolor = "black",  # Border color
                facecolor = "none",  # Fill color
                lw = 1.5,  # Line width of the border
                linestyle = (0, (2, 2)),  # Dashed style (dash length 5, gap 10)
            )
            ax.add_patch(circle)

    edgelist = [i for i in graph.edges if i[0] != 'C1a' and i[1] != 'C1a' if "flow" in graph.edges[i]]
    
    edge_colors= [cmap(norm(graph.edges[i]["flow"])) for i in edgelist]


    nx.draw_networkx_edges(graph, pos, edgelist, edge_color=edge_colors, node_size=sizeList,
                            arrowsize=15, connectionstyle="arc3", width = 3, alpha = 0.8)
    
    edgelist = [i for i in graph.edges if i[0] != 'C1a' and i[1] != 'C1a' if "flow" not in graph.edges[i]]
    
    nx.draw_networkx_edges(graph, pos, edgelist, edge_color="black", node_size=sizeList,
                            arrowsize=15, connectionstyle="arc3", width = 1)
    
    edgelist = [i for i in graph.edges if i[0] == 'C1a' or i[1] == 'C1a' if "flow" in graph.nodes[i[1]]]

    edge_colors= [cmap(norm(graph.nodes[i[1]]["flow"])) for i in edgelist]

    nx.draw_networkx_edges(graph, pos, edgelist, edge_color=edge_colors, node_size=sizeList,
                            arrowsize=15, connectionstyle="arc3", style = "--", width = 3)
    
    edgelist = [i for i in graph.edges if i[0] == 'C1a' or i[1] == 'C1a' if "flow" not in graph.nodes[i[1]]]

    nx.draw_networkx_edges(graph, pos, edgelist, edge_color="black", node_size=sizeList,
                                arrowsize=15, connectionstyle="arc3", width = 1)
    
    # Draw labels for species nodes only
    nx.draw_networkx_labels(graph, {n: pos[n] for n in species_nodes},
                            labels={node: node for node in species_nodes},
                            font_size=18)
    
    # Add a colorbar to represent flow intensity
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(flows)
    # plt.colorbar(sm, label='Flow Intensity')
    cbar = plt.colorbar(sm)
    cbar.set_label('Flow Intensity', fontsize=30)
    cbar.ax.tick_params(labelsize=30)
    
    # Set plot limits and equal aspect ratio
    x_values, y_values = zip(*pos.values())
    padding = 0  # Dynamic padding for better fit
    plt.xlim(min(x_values) - padding, max(x_values) + padding)
    plt.ylim(min(y_values) - padding, max(y_values) + padding)

    # Remove axis details for a cleaner look
    plt.axis('equal')
    plt.axis('off')
    
    # Adjust layout and show plot
    fig.tight_layout()
    plt.show()
# =============================================================================









# =============================================================================
def drawFormose():
    # Formose
    # -------------------------------------------------------------------------
    nameScenario = "formose"
    input_matrix, output_matrix = aux.readScenario(nameScenario)

    SM = output_matrix - input_matrix
    formose_species = ["C1a", "C2a", "C2b", "C3a", "C3b", 
                        "C3c", "C4a", "C4b", "C4c", "C5a", "C5b",
                        "C5c", "C5d", "C5e", "C6a", "C6b", "C6c", "C6d", "C6e",
                        "C7a", "C7b", "C7c", "C7d", "C7e", "C7f", "C8a", "C8b",
                        "C8c", "C8d"]
    
    flows = [ 588.62057179,  520.94924775,  460.95928687,   1.,         1000.,
      883.08574222,    4.32028878,   32.05029734,  749.93941776,    3.85930376,
       27.55519255,    3.4506451,    14.12527291,   10.33154373,    0.,
        2.2018831,    11.6648353,     8.30173182,    2.09486791,    7.38881373,
        7.38881373,    5.69303072,    5.69303072,    2.09486791,    2.09486791,
        2.09486791,    2.09486791,    1.,            1.,            1.,
        1.,            1.,            0.,            0.,            0.,
        1.09486791,  663.95681256]
    flows_normalizados = [float(i)/max(flows) for i in flows]
    reac_subgraph1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                      17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                      32, 33, 34, 35, 36, 37]  
    dict_flows = dict(zip(reac_subgraph1, flows_normalizados))
    G = createGraphFromStoichiometricMatrixWithFlowsSpecial(SM, dict_flows, formose_species)    
    
    new_G = replaceReactionNodesWithAdj2WithFlows(G)

    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
    # species_reac1 = []
    # species_reac1 = np.arange(SM.shape[0]) 
    species_reac1 =[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28]
    reac_subgraph1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                      17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                      32, 33, 34, 35, 36, 37]  
    food = [0]
    waste = []

    species_original = [i for i in G.nodes 
                          if G.nodes[i]["type"] == "species"]
    
    species_final1 = [species_original[i]
                      for i in species_reac1]
    
    reactions_original = [i for i in G.nodes 
                          if G.nodes[i]["type"] == "reaction"]
    reactions_modified = [i for i in new_G.nodes 
                          if new_G.nodes[i]["type"] == "reaction"]
    
    list_erased = []
    for idx, i in enumerate(reactions_original):
        if i not in set(reactions_modified):
            list_erased.append(idx)

    reac_subgraph_final1 = []
    for i in reac_subgraph1:
        if i not in set(list_erased):
            reac_subgraph_final1.append(i)
    reaction_final1 = [reactions_original[i]
                      for i in reac_subgraph_final1]

    new_food = [species_original[i]
                      for i in food]
    
    new_waste = [species_original[i]
                      for i in waste]        
    # VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

    
    subgraph = {"species": species_final1,
              "reaction": reaction_final1,
              "food": new_food,
              "waste": new_waste}

    # print(subgraph)
    a_sol = ["C2a", "C2b", "C3a", "C3b", "C3c", "C4a", "C4b", "C4c", "C5a",
             "C5b", "C5c", "C5d", "C5e", "C6a", "C6b", "C6c", "C6d", "C6e",
             "C7a", "C7b", "C7c", "C7e", "C7f", "C8a", "C8b", "C8c", "C8d"]

    plotChemicalReactionNetworkWithFlows(new_G, subgraph, a_sol, 1500, title="Chemical Reaction Network")
# =============================================================================










def main():
    drawFormose()
    
    
if __name__ == "__main__":
    main()













