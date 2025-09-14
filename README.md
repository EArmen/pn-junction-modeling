# p-n Junction Modeling: PBL Module with Python and Electrostatic Analogies

This repository contains open-source teaching materials for the paper *"Enhancing Student Understanding of p-n Junctions: A Problem-Based Learning Module with Python Modeling and Electrostatic Analogies"* (JOCSE, 2025) by Alfred V. Petrosyan and Armen S. Yepiskoposyan.

## Overview
The module uses Problem-Based Learning (PBL) to teach p-n junctions to undergraduates, including non-STEM majors. It integrates:
- **Electrostatic analogies** (parallel-plate capacitor).
- **Numerical modeling** via Python (Poisson equation solver).
- **Interactive Jupyter notebooks** with Matplotlib visualizations.
- **Empirical validation**: Score improvement from 62% to 85% (p < 0.01).

Materials support a three-level trajectory: intuition (Level 1), modeling (Level 2), analysis (Level 3).

## Quick Start
1. Install dependencies: `pip install -r requirements.txt`
2. Generate figures: `jupyter run generate_figures.ipynb`
3. Run tasks: `jupyter notebook tasks/practical_tasks.ipynb`
4. Explore supplementary: `jupyter notebook supplementary/`

## Structure
- **code/**: Core solver and plotting scripts.
- **tasks/**: PBL exercises (e.g., depletion width, E(x) plots).
- **supplementary/**: Advanced analysis (errors, temperature effects).
- **plots/**: Generated figures (e.g., Fig. 2: depletion_width.png).
- **docs/**: Full/short PDFs of the paper.

## Key Features
- **Reproducibility**: All code is self-contained; no external data needed.
- **Accessibility**: Analogies for non-STEM; Matplotlib for visualizations.
- **Comparison Table** (from paper Table 1):

| Tool          | Cost | Complexity | Accuracy | Interactivity |
|---------------|------|------------|----------|---------------|
| Python (ours) | Free | Medium     | High     | Medium        |
| MATLAB        | Paid | High       | High     | Medium        |
| COMSOL        | Paid | Very High  | Very High| Low           |
| NanoHUB       | Free | Low        | Medium   | High          |

## Citation
Petrosyan, A.V., & Yepiskoposyan, A.S. (2025). Enhancing Student Understanding of p-n Junctions... *Journal of Computational Science Education*.

For questions: petrosyanalfred58@gmail.com or earmen04@gmail.com

---
All materials © 2025, licensed under MIT.
