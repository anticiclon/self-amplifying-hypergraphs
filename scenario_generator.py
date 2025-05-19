# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 12:29:30 2024

@author: Trabajador
"""

import numpy as np
import random as rd
import os
from typing import List, Tuple

# =============================================================================
def partitions(x: int) -> List[List[int]]:
    """
    https://cristianbastidas.com/my-blog/en/algorithms/python/challenge/guide/2021/06/06/parts-of-number.html
    function for getting all partitions of a number
    Example: 
    - input:3
    - output: [[3], [2, 1], [1, 1, 1]]
    """
    # -------------------------------------------------------------------------
    def _partitions(ones: List[int], x: int, origin: set = None) -> List[List[int]]:
        if origin is None:
            origin = set()
        
        total = []
        
        for i in range(x, 0, -1):
            ayuda = []
            current_sum = 0
            
            while current_sum < x:
                if current_sum + i <= x:
                    ayuda.append(i)
                    current_sum += i
                else:
                    i -= 1
                    if i == 0:
                        break
            
            if current_sum < x:
                ayuda.extend([1] * (x - current_sum))
            elif current_sum > x:
                continue

            sorted_ayuda = tuple(sorted(ayuda))  # Use tuple for hashability
            if sorted_ayuda not in origin:
                total.append(ayuda)
                origin.add(sorted_ayuda)
                total.extend(_partitions(ayuda, x, origin))

        return total
    # -------------------------------------------------------------------------
    # Start with all ones as a base case
    all_ones = [1] * x
    parts = _partitions(all_ones, x, set())
    
    # Ensure all ones partition is included
    if all_ones not in parts:
        parts.append(all_ones)

    return parts
# =============================================================================

# =============================================================================
def turnIntoVectors(parts, n_max, rxn_order):
    """
    Convert a list of partitions into stoichiometric vectors.

    :param parts: List of partitions where each partition is a list of integers
    :param n_max: Maximum cluster size considered
    :param rxn_order: Highest reaction order allowed
    :return: Numpy array of stoichiometric vectors
    """
    # -------------------------------------------------------------------------
    def partitionToVector(part, n_max):
        stoich_vec = np.zeros(n_max, dtype=int)
        for cluster in part:
            if 1 <= cluster <= n_max:  # Ensure cluster sizes are valid
                stoich_vec[cluster - 1] += 1
        return stoich_vec
    # -------------------------------------------------------------------------

    # Use list comprehension for efficiency and readability
    s_vec = [
        partitionToVector(part, n_max) 
        for part in parts 
        if len(part) <= rxn_order and all(1 <= j <= n_max for j in part)
    ]
    
    return np.array(s_vec)
# =============================================================================

# =============================================================================
def createReactionMatrices(number_species: int, number_reactions: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate input and output matrices for chemical reactions based on stoichiometric vectors.

    :param number_species: Number of species in the system
    :param number_reactions: Number of reactions to simulate
    :return: Tuple of (input_matrix, output_matrix) where each is a numpy array
    """
    # Generate all possible partitions for the number of species
    all_partitions = partitions(number_species)
    # Convert partitions to stoichiometric vectors
    stoichiometric_vectors = turnIntoVectors(all_partitions, number_species, number_species)
    num_vectors = len(stoichiometric_vectors)
    
    if number_reactions > num_vectors * (num_vectors - 1):
        raise ValueError("Requested number of reactions exceeds possible unique reactions.")

    # Generate all possible reaction pairs (excluding self-pairing)
    possible_reactions = [(i, j) for i in range(num_vectors) for j in range(num_vectors) if i != j]

    # Randomly select reactions
    selected_reactions = rd.sample(possible_reactions, number_reactions)

    # Create input and output matrices
    input_matrix = np.array([stoichiometric_vectors[i] for i, _ in selected_reactions]).T
    output_matrix = np.array([stoichiometric_vectors[j] for _, j in selected_reactions]).T

    return input_matrix, output_matrix
# =============================================================================

# =============================================================================        
def saveMatrices(mMinus, mPlus, name):
    np.savetxt("scenarios/" + name + "_minus.txt", mMinus, fmt = '%i')
    np.savetxt("scenarios/" + name + "_plus.txt", mPlus, fmt = '%i')
# =============================================================================

# =============================================================================
def scenarioGeneratorV1(number_species: int, number_reactions: int, version: int) -> None:
    """
    Generate a chemical reaction scenario including input and output matrices.

    :param number_species: The number of chemical species in the scenario.
    :param number_reactions: The number of chemical reactions to simulate.
    :param version: The version identifier for this scenario.
    :raises ValueError: If any input parameter is not a positive integer.
    """
    if not all(isinstance(param, int) and param > 0 for param in [number_species, number_reactions, version]):
        raise ValueError("All parameters must be positive integers.")

    input_matrix, output_matrix = createReactionMatrices(number_species, number_reactions)
    scenario_name = f"scenario_e{number_species}_r{number_reactions}_v{version}"
    saveMatrices(input_matrix, output_matrix, scenario_name)
# =============================================================================

# =============================================================================
def scenarioGeneratorV2(rows: int, cols: int, density_factor: float, version: int, max_value: int) -> None:
    """
    Generate input and output matrices for a scenario based on random selection.

    :param rows: Number of rows in matrices.
    :param cols: Number of columns in matrices.
    :param density_factor: Factor to control the density of selected elements, higher means fewer selections.
    :param version: Version of the scenario for file naming.
    :param max_value: Maximum value to fill in the initial matrix.
    :raises ValueError: If any input parameter is invalid.
    """
    if not (isinstance(rows, int) and isinstance(cols, int) and rows > 0 and cols > 0):
        raise ValueError("Rows and columns must be positive integers.")
    if not isinstance(density_factor, (int, float)) or density_factor <= 0:
        raise ValueError("Density factor must be a positive number.")
    if not isinstance(max_value, int) or max_value < 0:
        raise ValueError("Max value must be a non-negative integer.")
    if not isinstance(version, int) or version < 0:
        raise ValueError("Version must be a non-negative integer.")

    # Create random matrix
    matrix = np.random.randint(0, max_value + 1, size=(rows, cols))

    # Number of elements to select for each matrix
    n_elements = int(rows * cols / density_factor)

    # Generate unique pairs for input matrix
    input_pairs = set()
    while len(input_pairs) < n_elements:
        pair = tuple(int(np.random.randint(0, dim)) for dim in (rows, cols))
        # print("Type of pair:", type(pair))
        # print("Content of pair:", pair)
        input_pairs.add(pair)
        
    # Generate unique pairs for output matrix, not overlapping with input
    output_pairs = set()
    while len(output_pairs) < n_elements:
        pair = tuple(int(np.random.randint(0, dim)) for dim in (rows, cols))
        if pair not in input_pairs:
            # print("Type of pair:", type(pair))
            # print("Content of pair:", pair)
            output_pairs.add(pair)

    # Fill matrices based on selected pairs
    input_matrix = np.zeros((rows, cols), dtype=int)
    output_matrix = np.zeros((rows, cols), dtype=int)
    for i in range(rows):
        for j in range(cols):
            if (i, j) in input_pairs:
                input_matrix[i, j] = matrix[i, j]
            elif (i, j) in output_pairs:
                output_matrix[i, j] = matrix[i, j]

    # Save matrices
    os.makedirs('scenarios', exist_ok=True)
    np.savetxt(f"scenarios/n{rows}m{cols}d{int(density_factor)}max{max_value}v{version}_minus.txt", input_matrix, fmt="%d")
    np.savetxt(f"scenarios/n{rows}m{cols}d{int(density_factor)}max{max_value}v{version}_plus.txt", output_matrix, fmt="%d")
# =============================================================================

