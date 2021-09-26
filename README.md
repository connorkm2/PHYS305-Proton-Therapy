# PHYS305 Computational Physics
This ccomputational physics project was completed as part of a course at the University of Liverpool.

## Proton Therapy

This simulation code simulates the delivery of a proton beam to a possible patient. 

Simple particle-matter interactions can be observed including the production of Spread Out Bragg Peak, similar to that found in practice using physical materials but here simulated by appropriately weighting the beam intensity.
We also include 2-dimensional gaussian surface which provides flattening in the x-y plane, to ensure homogeneous delivery of the dose.

Various metrics are also output including the clincally used dose volume histograms which represent the amount of dose delivered to a specified area compared to that delivered outside of that area, to allow for the user to determine the safety/effectiveness of the delivery plan. 
