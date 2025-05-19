# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:12:52 2024

@author: Trabajador
"""

import autocatalytic_cores_lib as code
import auxiliary_functions as aux
import pandas as pd

def main():

    
    # Formose
    # -------------------------------------------------------------------------
    nameScenario = "formose"
    input_matrix, output_matrix = aux.readScenario(nameScenario)
    
    excelfile = "output/" + nameScenario + ".xlsx"
    txtfile = "output/" + nameScenario + ".txt"
    histogramfile = "output/" + nameScenario + ".png"

    SM = output_matrix - input_matrix
    formose_species = ["C1a fomaldehyd", "C2a", "C2b", "C3a", "C3b", 
                        "C3c dihydroxy acetone", "C4a", "C4b", "C4c", "C5a", "C5b",
                        "C5c", "C5d", "C5e", "C6a", "C6b", "C6c", "C6d", "C6e",
                        "C7a", "C7b", "C7c", "C7d", "C7e", "C7f", "C8a", "C8b",
                        "C8c", "C8d"]
        
    df = code.computeAutocatalyticCores(SM = SM, excelfile = excelfile , names_sp = formose_species, txtfile=txtfile, histogramfile = histogramfile)
    # -------------------------------------------------------------------------




    # E coli 
    # # -------------------------------------------------------------------------
    # file = "scenarios"
    # nameScenario = "e_coli"
    # input_matrix, output_matrix = aux.readScenario(nameScenario)

    # excelfile = "output/" + nameScenario + ".xlsx"
    # txtfile = "output/" + nameScenario + ".txt"
    # histogramfile = "output/" + nameScenario + ".png"


    # file_path_species = "./scenarios/e_coli_species.csv"
    # ecoli_species = list(pd.read_csv(file_path_species, delimiter=';')["id"])

    # SM = output_matrix - input_matrix
    # df = code.computeAutocatalyticCores(SM = SM, excelfile = excelfile , names_sp = ecoli_species, txtfile=txtfile, histogramfile = histogramfile)
    # # # -------------------------------------------------------------------------








if __name__ == "__main__":
    main()
    
    