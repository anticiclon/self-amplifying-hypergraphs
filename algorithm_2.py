# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 01:52:02 2024

@author: Trabajador
"""

# import networkx as nx
import numpy as np
import gurobipy as gb
import time
from typing import Any
import auxiliary_functions as aux


# =============================================================================
def computeSubhypergraphWithGreaterMAF(
        input_matrix: Any,
        output_matrix: Any,
        scenario_name: str,
        time_limit_iteration: int,
        max_steps = 1000,
        name_nodes = ""
    ) -> None:
    """
    Calculate the amplification factor of subgraphs and log the results to a file.

    Args:
        input_matrix: The input incidence matrix of the hypergraph.
        output_matrix: The output incidence matrix of the hypergraph.
        time_limit_iteration: Time limit for each algorithm iteration.
        max_steps: The maximum number of iterations.
        scenario_name: Name for the output file.
    """
    
    NULL_THRESHOLD = 0.1  # Threshold for determining null rows/columns
    
    # Preprocessing: Remove null rows and columns
    try:
        input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(
            input_matrix, output_matrix, NULL_THRESHOLD
        )
    except Exception as e:
        raise RuntimeError(f"Error in preprocessing matrices: {e}")
    
    # Mapping rows and columns
    row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
    
    # Decompose into independent components
    components = aux.giveMeMatrixByComponent(input_modified, output_modified)
    
    # Create a list of component dictionaries
    component_dicts = [
        {
            "input": comp["input"],
            "output": comp["output"],
            "nodes": [row_mapping[idx] for idx in comp["nodes"]],
            "arcs": [column_mapping[idx] for idx in comp["arcs"]],
        }
        for comp in components
    ]
    
    # Solve for each component
    solutions = []
    for component in component_dicts:
        try:
            solution = computeMAFinSubhypergraph(
                component["output"], component["input"], max_steps, time_limit_iteration
            )
        except Exception as e:
            raise RuntimeError(f"Error in solving subhypergraph: {e}")
        
        solutions.append({
            "x": solution[0],
            "alpha": solution[1],
            "step": solution[2],
            "alphaDict": solution[3],
            "a": solution[4],
            "z": solution[5],
            "time": solution[6],
            "nodes": component["nodes"],
            "arcs": component["arcs"],
            "dict_sol": solution[7],
        })
    
    # Find the solution with the MAF (alpha)
    best_solution = max(solutions, key=lambda sol: sol["alpha"])

    # Prepare output data
    output_lines = []
    output_lines.append(f"{solutions}\n")
    output_lines.append(f"Total time:\n{best_solution['time']}\n")
    output_lines.append(f"Q dimension:\n{input_matrix.shape}\n")
    output_lines.append(f"self-sufficiency Q:\n{aux.checkSelfSufficiently(input_matrix, output_matrix)[2]}\n")
    output_lines.append(f"Number of self-amplifying nodes (a):\n{len(best_solution['a'])}\n")
    output_lines.append(
        f"Self-amplifying nodes (a):\n{' '.join([str(name_nodes[best_solution['nodes'][i]]) for i in best_solution['a']])}\n"
        )
    # name_nodes
    output_lines.append(f"nodes (numbers):\n{[best_solution['nodes'][i] for i in best_solution['a']]}\n")
    output_lines.append(f"arcs (numbers):\n{[best_solution['arcs'][z] for z in best_solution['z']]}\n")
    output_lines.append(f"MAF:\n{best_solution['alpha']}\n")
    
    # Record the resulting arcs
    selected_arcs = [best_solution["arcs"][z] for z in best_solution["z"]]
    output_ssp = output_matrix[:, selected_arcs]
    output_ssm = input_matrix[:, selected_arcs]
    
    # Map the rows and columns for the final output
    new_column_mapping = {i: selected_arcs[i] for i in range(len(selected_arcs))}
    alter_column_mapping = {i: best_solution["z"][i] for i in range(len(best_solution["z"]))}
    row_mapping = {i: i for i in range(output_ssm.shape[0])}
        
    if name_nodes == "":
        arc_details = aux.recordArcsAltered(output_ssp, output_ssm, best_solution["x"], row_mapping, new_column_mapping, alter_column_mapping)
    else:
        arc_details = aux.recordArcsAlteredWithNames(output_ssp, output_ssm, best_solution["x"], row_mapping, new_column_mapping, alter_column_mapping, name_nodes)

    output_lines.append("Q_result arcs:\n")
    output_lines.extend(arc_details)
    
    # Write output to a file
    output_file = f"output\{scenario_name}_algorithm_2.txt"
    try:
        with open(output_file, 'w') as f:
            f.writelines(output_lines)
    except Exception as e:
        raise RuntimeError(f"Error writing to file {output_file}: {e}")
        
    return solutions[0], arc_details
# =============================================================================


# =============================================================================
def computeMAFinSubhypergraph(output_matrix, input_matrix, t_max, time_limit_iteration, accuracy = 1e-3):
    """
    Calculate the MAF for a subhypergraph.
    
    Parameters:
    - output_matrix (np.ndarray): The output incidence matrix of the hypergraph.
    - input_matrix (np.ndarray): The input incidence matrix of the hypergraph.
    - t_max (int): Maximum number of iterations for convergence.
    - time_limit_iteration (float): Time limit for each optimization step in seconds.
    - accuracy (float, optional): precision solution
    
    Returns:
    - Tuple[np.ndarray, float, int, dict, list, list, float, dict]: 
        (optimal_intensity, final_amplification_factor, steps, amplification_factors, active_nodes, active_arcs, total_time, iteration_data)
    """
    # Parameters
    # ---------------------------
    # Number Nodes (int) and Number Arcs (int)
    num_nodes, num_arcs = output_matrix.shape
    # Nodes (list)
    nodes = range(num_nodes)
    # Arcs (list)
    arcs = range(num_arcs)
    # Alpha_0 (float)
    x_0 = np.ones(num_arcs)
    #
    alpha_0 = np.min([sum(output_matrix[v, a] * x_0[a] 
                          for a in arcs)
                      /
                      sum(input_matrix[v, a] * x_0[a] 
                          for a in arcs) 
                      for v in nodes])
    # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    if alpha_0 < 0.001:
        x_0 = np.random.randint(1, 100, size = num_arcs)
        #
        alpha_0 = np.min([sum(output_matrix[v, a] * x_0[a] 
                          for a in arcs)
                      /
                      sum(input_matrix[v, a] * x_0[a] 
                          for a in arcs) 
                      for v in nodes])
    # --------------------------------------
    # -------------------------------------------------------------------------
    def computeMAFinSubhypergraphFixed(alpha_0, time_limit_iteration):

        # Parameters
        upper_bound = 1000
        bigM2 = upper_bound  # Big-M for constraint relaxation
    
        # Initialize model
        model = gb.Model("Amplification_Factor_Model_SN")
        model.Params.OutputFlag = 0  # Suppress console output
        model.Params.TimeLimit = time_limit_iteration
        model.Params.MIPGap = 0.00  # Set to zero for exact solutions, might be too strict for practical use
    
        # Define variables
        intensities = model.addVars(num_arcs, lb=0, ub=upper_bound, name="arc_intensities")
        amplification_factor = model.addVar(name = "amplification_factor")
        is_self_amplifying = model.addVars(num_nodes, vtype = gb.GRB.BINARY, name = "is_self_amplifying")
        is_active_arc = model.addVars(num_arcs, vtype = gb.GRB.BINARY, name = "is_active_arc")

        # Objective function
        model.setObjective(amplification_factor, gb.GRB.MAXIMIZE)
        # Constraints
        # --------------------------------------
        model.addConstrs(
            (
            amplification_factor <= gb.quicksum(output_matrix[v, a] * intensities[a] 
                                  for a in arcs)
                    - alpha_0 * gb.quicksum(input_matrix[v, a] * intensities[a]
                                            for a in arcs) 
                    + alpha_0 * sum(input_matrix[v, :])*upper_bound*(1 - is_self_amplifying[v])
                    for v in nodes),
            name = "name1")
        #
        model.addConstrs(
            (gb.quicksum(input_matrix[v, a] * intensities[a] 
                          for a in arcs) 
            >= is_active_arc[v]
            for v in nodes 
            if sum(input_matrix[v, :]) > 0),
            name = "name2")
        #
        model.addConstrs(
            (is_active_arc[a] <= gb.quicksum(is_self_amplifying[v] for v in nodes
                                  if output_matrix[v, a] > 0) 
              for a in arcs), 
            name = "name8")
        #
        model.addConstrs(
            (is_active_arc[a] <= gb.quicksum(is_self_amplifying[v] for v in nodes 
                                  if input_matrix[v, a] > 0) 
              for a in arcs), 
            name = "name9")
        #
        model.addConstrs(
            (intensities[a] <= bigM2 * is_active_arc[a]
              for a in arcs), 
            name = "name10")
        #
        model.addConstrs(
            (is_active_arc[a] <= intensities[a] 
              for a in arcs), 
            name = "name11")
        #   
        model.addConstr(gb.quicksum(is_self_amplifying[v] 
                                for v in nodes) >= 1, 
                    name = "name12")
        #
        model.addConstr(gb.quicksum(is_active_arc[a] 
                                for a in arcs) >= 1,
                    name = "name13")
        # --------------------------------------      
  
        # Solve the model
        model.optimize()
    
        if model.status != gb.GRB.OPTIMAL:
            model.computeIIS()
            model.write("infeasibility_report.ILP")
            print("Model infeasible. Check infeasibility report.")
            raise ValueError("Optimization failed: model is infeasible")
    
        # Extract solution
        optimal_intensities = np.array([intensities[r].X for r in range(num_arcs)])
        optimal_amplification_factor = amplification_factor.X
        active_nodes = [v for v in range(num_nodes) if is_self_amplifying[v].X > 0.5]
        active_arcs = [a for a in range(num_arcs) if is_active_arc[a].X > 0.5]
    
        return optimal_intensities, optimal_amplification_factor, active_nodes, active_arcs, model.MIPGap, model.NumVars, model.NumConstrs
    # -------------------------------------------------------------------------
        
    stop = False
    step = 0
    alpha_dict = {}
    # alphabar = 10000
    current_alpha = alpha_0
    previous_alpha = 0  # Starting with 0 for comparison
    start_time=time.time()
    
    iteration_data = {}
    
    counter = 1
    while not stop:
        print(counter, current_alpha)
        iteration_start_time = time.time()

        # Solve optimization model for this iteration
        optimal_intensities, alpha_bar, active_nodes, active_arcs, gap, num_vars, num_constrs = computeMAFinSubhypergraphFixed(current_alpha, time_limit_iteration)
        
        # Update current alpha considering only active nodes
        current_alpha = np.min([sum(output_matrix[v, a] * optimal_intensities[a]
                        for a in arcs)
                        /
                        sum(input_matrix[v, a] * optimal_intensities[a] 
                        for a in arcs) 
                    for v in active_nodes])
        
        counter = counter + 1
        
        iteration_end_time = time.time()
        
        iteration_data[step] = {
            "intensities": optimal_intensities,
            "active_nodes": active_nodes,
            "active_arcs": active_arcs,
            "alpha_bar": alpha_bar,
            "gap": gap,
            "variables": num_vars,
            "constraints": num_constrs,
            "step": step,
            "alpha": current_alpha,
            "time": iteration_end_time - iteration_start_time
        }
        
        
        if len(active_nodes) < 1 or np.abs(alpha_bar) < accuracy or step >= t_max or np.abs(current_alpha - previous_alpha) < accuracy:
            stop = True
            alpha_dict[step] = current_alpha
            return optimal_intensities, current_alpha, step, alpha_dict, active_nodes, active_arcs, time.time() - start_time, iteration_data
        else:
            alpha_dict[step] = current_alpha
            previous_alpha = current_alpha
            step += 1
# =============================================================================





















