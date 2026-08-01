# RIT Brake System Analysis
This repository collects all of the most important brakes and driver controls systems models for ease of use and version control. At the time of initial commit, it contains models and other code from over 15 years of RIT Racing's history, the original intent and development of which are unknown to the initial committer. Questions can be directed to Oliver Owen (oliver.owen.pers@gmail.com). 

# Script objectives:
## F34BrakeDataPreprocessor.m
Takes as many sets of car data that the user wishes as input and processes the data to remove peaks and long periods where the vehicle is stationary. Optionally splits the input data set into individual files for each active driving stint. 

## BrakeCoeffOptimizer.m
Takes as many sets of car data that the user wishes as input, applies the same padFrac and front and rear cooling coefficient fit parameters to each dataset, and uses a nonlinear least-squares algorithm to optimize these parameters for the best fit to vehicle data. Fits four different padFrac functions: linear with temperature only, linear with temperature and pressure but no interaction, linear with temperature and pressure and interaction, and quadratic in temperature and linear in pressure, and selects the best function based on Akaike’s corrected Information Criterion. Works best with preprocessed data.
