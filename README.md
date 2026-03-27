# A Reproducible Python-Based Approach for Teaching the Equilibrium p–n Junction through analytical and numerical modelling

This repository contains the materials associated with the article:

**"A reproducible Python-based approach for teaching the equilibrium p–n junction through analytical and numerical modelling"**

---

## Overview

This project provides a classroom-ready framework for teaching the equilibrium p–n junction using:

- Electrostatic analogy  
- Analytical modelling  
- Numerical solution of the Poisson equation  
- Python-based simulations  

---

## Repository Structure

pn-junction-modeling/  
├── generate_all_figures.py   # main script  
├── README.md  
├── requirements.txt  
├── .gitignore  
├── LICENSE  
├── src/                     # core physics models  
├── notebooks/               # Jupyter notebooks  
├── plots_article/           # generated figures  
└── tables_article/          # generated tables  

---

## How to Reproduce Results

Run the main script:

python generate_all_figures.py

This will generate all figures and tables from the article.

---

## Notes

- Figures 1–2 are computed from physical models  
- Figures 3–6 are reproduced for visualization consistency  
- Output files are saved into:
  - plots_article/  
  - tables_article/  

---

## Requirements

- numpy  
- pandas  
- matplotlib  
- scipy  

---

## License

MIT License
