# Transient_Heat_Transfer_Analysis
MATLAB-based transient heat transfer analysis project, completed in November 2022.

2D Transient Heat Transfer Analysis

This project presents the numerical solution of the 2D transient heat transfer equation using MATLAB. It implements various methods to solve the problem and investigates the stability conditions.

Overview:

	Key methods used:

	* Explicit Finite Difference Method: Solves the 2D heat equation and examines stability via the CFL condition.
	* PSOR Method: Compares performance and accuracy with the explicit method.
	* Laasonen (Implicit) Method: Provides a stable solution even for larger time steps.

 	Features:

	* Stability analysis of the explicit method.
	* Comparison between explicit, PSOR, and Laasonen methods.

	Contents:

	Report: Full explanation of methods and results.
	Code: MATLAB .m files for explicit, PSOR, and Laasonen methods.
	Results: Plots comparing the methods.

How to Run:
	* Run the transient_heat_transfer_2D.m script for the explicit method.
	* Execute Guass_seidel.m and lassonnen.m for PSOR and implicit solutions.
	* Generated plots will be saved in the Results/ folder.
