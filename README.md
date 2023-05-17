# Surface-to-Volume-Scaling
A MATLAB wokflow for analyzing cell line profiles in light microscopy images for quanitfying amount of cell membrane relative to size. The code is to be compiled in four modular parts.

## read_dv.m
Reads every DeltaVision image in specified folder and identifies possible cells using built-in 'TwoStage' method. 

## curatedcutouts.m
Saves each potential cell as a 'cutout' and 'key' image. These allow manual curation to help eliminate invalid cells from the circle finder, adjust position and radii.

## cookiecutter.m
Creates eclipse images--for every curated cell, there is an image that has all cells other than the cell of interest blacked out to prevent line profile collisions.

## lineprofiler.m
Takes 360 line profiles on each eclipse image and plots analysis.
![alt text](https://i.imgur.com/MTHVZ4R.jpeg)
![alt text](https://i.imgur.com/1bPc4aq.jpeg)
![alt text](https://i.imgur.com/cptcR2w.jpeg)
![alt_text](https://i.imgur.com/Vsa0hdf.jpeg)
