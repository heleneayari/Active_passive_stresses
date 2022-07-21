# Active_passive_stresses

This repository is intended to run the analysis which has been used in
the following PRE paper:
https://arxiv.org/abs/2203.06475

It was run on Matlab 2017a.

+ Run the file: calcul_forces_fft2D_adhesion_3D_rect_2gaussians.m
to calculate the forces exerted by the cell on the substrate using the displacement field coming from the Finite Element Simulation (FEM).
+ Run the file : Calculate_stress_cell_rect_ad_2gaussians.m to calculate via ISM the deformation and elastic stress in the cell using again the displacement field coming from the FEM simulation.
+ Run BISM_cartesian_cell_rect_ad_2gaussians_BC.m to calculate the total Stress in the monolayer (Passive + Active) using the BISM algorithm developed by Nier et al.
+ Run compare_BISM_comsol.m to compare ISM BISM and FEM calculations.


You should finally retrieve the different figures published in the article.





  