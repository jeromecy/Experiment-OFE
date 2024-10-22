# Optimal Design for On-Farm Strip Trials—Systematic or Randomised?

## Overview

This repository contains the simulation work and manuscript for the paper titled "Optimal design for on-farm strip trials—systematic or randomised?" published in Field Crops Research. The paper can be accessed here.

## Abstract

**🌾 Context or Problem**  
Randomised designs are often preferred by agronomists and biometricians. However, for on-farm trials, the choice may depend on the experiment's objective. If the goal is to create a prescription map for each plot in a large strip trial, a systematic design might be better, even though it gets less attention.

**🎯 Objective**  
This study evaluates how well systematic designs with geographically weighted regression (GWR) models handle spatial variation and estimate continuous treatment effects in large strip trials through simulations.

**🔬 Methods**  
We used a hierarchical model with spatially correlated random parameters to generate simulated data for various large strip on-farm trial scenarios. GWR models analyzed the data assuming both linear and quadratic responses of yield to treatment effects.

**📊 Results**  
- **Quadratic Response**: Systematic design achieved lower mean squared errors (MSEs) with GWR.
- **Linear Response**: No significant difference in MSE between systematic and randomised designs, regardless of spatial variation.

**🏆 Conclusions**  
Systematic designs are superior for producing smooth spatial maps of optimal input levels for quadratic response models in large strip trials, even with significant spatial variation. Fixed bandwidths in GWR analysis should be selected based on plot configurations. For large strip trials, systematic designs provide better estimates of spatially-varying treatment effects than randomised designs.

**🌟 Implications**  
These findings offer practical recommendations for designing large strip trials, contributing valuable insights for improving their efficacy and planning.

## Repository Contents

- **Simulation Code**: Scripts used to generate and analyze the simulated data.
- **Manuscript**: The full manuscript of the published paper.
- **Data**: Simulated datasets used in the study.

## Citation

If you use this repository in your research, please cite the paper as follows:

```bibtex
@article{cao2024optimal,
  title={Optimal design for on-farm strip trials—systematic or randomised?},
  author={Cao, Zhanglong and Brown, Jordan and Gibberd, Mark and Easton, Julia and Rakshit, Suman},
  journal={Field Crops Research},
  volume={318},
  pages={109594},
  year={2024},
  doi ={https://doi.org/10.1016/j.fcr.2024.109594},
  publisher={Elsevier}
}
```
