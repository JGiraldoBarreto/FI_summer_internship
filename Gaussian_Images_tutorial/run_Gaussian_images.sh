#!/bin/bash

##This program generates a Gaussian Images file, given a grid from where the images will be generated, and a Matrix whose components are the numner of images per model. 

## About the invers temperature value: It is possible that the initial matrix does not fulfill the low temperature limit, so in that case it is important to reach that limit by scaling the number of images. This can be done if the "temperature" is lowered by an inverse temperature value. The higher the inverse temperature value, the lower the temperature.

## About the grid: Is an array with shape Mx2, where M is the number of models that "generates" the images, and 2 is the dimension of the coordinate space (in this case (x,y)). The grid could cover all the space (e.g [ [1,1], [1,2],...[1,20],[2,1],...,[20,20] ] for a 20x20 grid), so the images are generate from a full landscape, or the grid can be a particular set of coordinates (for examples the coordinate points of a particular N-nodes path), so the images are generated specifically for that points "from the path".

## Syntax: python3 Gaussian_images.py "grid_file" "Matrix_with_number_of_images_per_model_file" "Inverse_Temp_value"

python3 Gaussian_images.py data/grid data/Num_Images_grid-3well 3.5

