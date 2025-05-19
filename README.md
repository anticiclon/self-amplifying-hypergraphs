# Self-Amplifying Hypergraphs

This repository contains the code accompanying the paper  
**"Identifying Self-Amplifying Hypergraph Structures through Mathematical Optimization."**

## Overview

The repository provides implementations of three core algorithms discussed in the paper:

- **Algorithm 1:** Compute the maximal Amplification factor (MAF) in a hypergraphs.
- **Algorithm 2:** Extracts subhypergraphs with greater self-amplifying factor (MAF).
- **Core Enumeration:** Enumerates all autocatalytic cores and evaluates their MAFs.

---

## Repository Structure

### 🔁 `algorithm_1.py`

#### Function: `recordAmplificationData(input_matrix, output_matrix, nameScenario, time_limit_iteration, name_nodes, accuracy, max_steps)`

- Implements **Algorithm 1** from the paper.
- Expects self-sufficient input and output incidence matrices stored in the `scenarios/` folder.
- Writes results to the `output/` folder.
- Example usage: See [`main_algorithm_1.py`](./main_algorithm_1.py) for applications to the **Formose** and **E. coli** scenarios.

---

### 🧬 `algorithm_2.py`

#### Function: `computeSubhypergraphWithGreaterMAF(input_matrix, output_matrix, scenario_name, time_limit_iteration, max_steps, name_nodes, accuracy)`

- Implements **Algorithm 2**, including:
  - Self-sufficiency simplification.
  - Decomposition into independent components.
- Input matrices should be placed in the `scenarios/` folder.
- Output is saved in `.txt` files within the `output/` folder.
- Example usage: See [`main_algorithm_2.py`](./main_algorithm_2.py) for applications to the **Formose** and **E. coli** scenarios.

---

### 🧠 `autocatalytic_cores_lib.py`

#### Function: `computeMAF(output_matrix, input_matrix, max_steps, time_limit_iteration)`

- Enumerates autocatalytic cores using the method described in  
  [this paper](https://link.springer.com/article/10.1007/s10910-024-01576-x).
- Computes the Maximal Amplification Factor (MAF) of each core via Algorithm 1.
- Reads scenarios from `scenarios/`, and writes output to `output/`.
- Example usage: See [`main_cores.py`](./main_cores.py) for applications to the **Formose** and **E. coli** scenarios.

---

## 🧪 Experiments

The `experiments/` folder contains scenarios generated using the **SBbadger** tool, described in  
[this paper](https://academic.oup.com/bioinformatics/article/38/22/5064/6701964?login=true).  
These scenarios were used for the experimental results reported in our study.

---

## 📊 Visualization

The file `drawing.py` contains plotting functions for visualizing results.  
It includes an example that visualizes the Formose solution produced by Algorithm 1.

---

## 🧰 Scenario Generation

The `scenario_generator/` folder includes helper functions for:

- Generating random scenarios.
- Running tests on new instances.

---




















