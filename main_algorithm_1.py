# -*- coding: utf-8 -*-

import algorithm_1 as algo
import auxiliary_functions as aux
import pandas as pd
import numpy as np

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
    
    input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix)
    row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
    rows_kept = np.arange(input_modified.shape[0]) 
    name_species = [str(formose_species[row_mapping[i]]) for i in rows_kept]

    algo.recordAmplificationData(input_modified, output_modified, nameScenario, time_limit_iteration, name_species)
    # -------------------------------------------------------------------------

    # # e_coli
    # # -------------------------------------------------------------------------
    # nameScenario = "e_coli"
    # input_matrix, output_matrix = aux.readScenario(nameScenario)
    # time_limit_iteration = 1000
    # file_path_species = "./scenarios/e_coli_species.csv"
    # ecoli_species = list(pd.read_csv(file_path_species, delimiter=';')["id"])
    # # print(ecoli_species)
    # # print(len(ecoli_species))
    # input_modified, output_modified, null_rows, null_columns = aux.removeNullRowsAndColumns(input_matrix, output_matrix)
    # row_mapping, column_mapping = aux.mappingRowsAndColumns(input_matrix, null_rows, null_columns)
    # rows_kept = np.arange(input_modified.shape[0]) 
    # name_species = [str(ecoli_species[row_mapping[i]]) for i in rows_kept]
    # # print('---')
    # # print(name_species)
    # # print(len(name_species))
    # algo.recordAmplificationData(input_modified, output_modified, nameScenario, time_limit_iteration, ecoli_species)
    # # -------------------------------------------------------------------------

    # test
    # -------------------------------------------------------------------------
    input_matrix = np.array([[2., 0., 0., 0.],
                              [0., 1., 0., 0.],
                              [2., 0., 1., 0.],
                              [0., 0., 0., 1.]])
    
    output_matrix = np.array([[0., 1., 0., 1.],
                              [1., 0., 1., 0.],
                              [0., 1., 0., 1.],
                              [0., 0., 1., 0.]])
    time_limit_iteration = 1000
    nameScenario = "test_1"
    algo.recordAmplificationData(input_matrix, output_matrix, nameScenario, time_limit_iteration)
    # -------------------------------------------------------------------------



if __name__ == "__main__":
    main()
    

