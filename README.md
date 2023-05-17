# Surface-to-Volume-Scaling
A MATLAB wokflow for analyzing cell line profiles in light microscopy images. The code is to be compiled in four modular parts.

## read_dv.m
Reads every DeltaVision image in specified folder and identifies possible cells using built-in 'TwoStage' method. 

## curatedcutouts.m
Saves each potential cell as a 'cutout' and 'key' image. These help eliminate invalid cells from the circle finder and adjust position and radii.

## cookiecutter.m
Creates eclipse images--for every curated cell, there is an image that has all cells other than the cell of interest blacked out to prevent line profile collisions.

## lineprofiler.m
Takes 360 line profiles on each eclipse image and plots analysis.
![alt text](https://imgur.com/a/zYwEhOV)
![alt text](https://imgur.com/a/wA9JRDj)
![alt text](https://imgur.com/a/4uX2Rwq)
![alt text](https://imgur.com/a/kEiSpW2)
