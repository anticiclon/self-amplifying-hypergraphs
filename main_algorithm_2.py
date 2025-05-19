# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 01:03:40 2024

@author: Trabajador
"""

import algorithm_2 as algo
import auxiliary_functions as aux
import pandas as pd

def main():

    # Formose
    # -------------------------------------------------------------------------
    nameScenario = "formose"
    input_matrix, output_matrix = aux.readScenario(nameScenario)
    time_limit_iteration = 1000
    formose_species = ["C1a fomaldehyd", "C2a", "C2b", "C3a", "C3b", 
                        "C3c dihydroxy acetone", "C4a", "C4b", "C4c", "C5a", "C5b",
                        "C5c", "C5d", "C5e", "C6a", "C6b", "C6c", "C6d", "C6e",
                        "C7a", "C7b", "C7c", "C7d", "C7e", "C7f", "C8a", "C8b",
                        "C8c", "C8d"]
    algo.computeSubhypergraphWithGreaterMAF(input_matrix, output_matrix, nameScenario, time_limit_iteration, name_nodes = formose_species)
    # -------------------------------------------------------------------------

    # # E coli
    # # -------------------------------------------------------------------------
    # nameScenario = "e_coli"
    # input_matrix, output_matrix = aux.readScenario(nameScenario)
    # time_limit_iteration = 1000
    # file_path_species = "./scenarios/e_coli_species.csv"
    # ecoli_species = list(pd.read_csv(file_path_species, delimiter=';')["id"])
    # algo.computeSubhypergraphWithGreaterMAF(input_matrix, output_matrix, nameScenario, time_limit_iteration, ecoli_species)
    # # -------------------------------------------------------------------------






if __name__ == "__main__":
    main()





    




