#!/bin/bash

for i in `seq 0 12`
do
	python find_surface_contour.py $i
	python 3d_render.py $i
done
