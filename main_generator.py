# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 14:56:18 2024

@author: Trabajador
"""

import scenario_generator as gene


def main():

    # -------------------------------------------------------------------------
    number_species = 100
    number_reactions = 100
    version = 32
    gene.scenarioGeneratorV1(number_species, number_reactions, version)
    # -------------------------------------------------------------------------

    # -------------------------------------------------------------------------
    rows = 20
    cols = 40
    density_factor = 5
    version = 38
    max_value = 100
    gene.scenarioGeneratorV2(rows, cols, density_factor, version, max_value)
    # -------------------------------------------------------------------------



if __name__ == "__main__":
    main()
    








