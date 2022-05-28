# Nematic
Construct a dense liquid crystal of non-overlapping spherocylinders via Go and Mathematica

spherocyl_random_construct.nb is a notebook that once run, stuffs particles of a chosen aspect ratio into a box, trying to reach
a chosen volume fraction while always respecting that these particles should have a certain degree of orientational order. 
This is a starting configuration for a denser liquid crystal. The routine can only reach so high a volume fraction, and 
the point is to use the output in the next step, which is an implementation of MCM2 (https://pubs.acs.org/doi/10.1021/je500119r).
This next step (the .go file) reads the output from mathematica, and then slowly raises the volume of the particles to increase
the volume fraction of the whole collection. The particles are treate as inifintely hard, so any overlaps are "undone" by
translation and rotation, until a target final volume fraction is reached.

The algorithm in the notebook is adapted from Henrik Schumacher's answer in https://mathematica.stackexchange.com/questions/159484/creating-random-configurations-of-spherocylinders-or-cylinders

Process:

- Open "spherocyl_random_construct.nb", select your parameters (aspect ratio, initial volume fraction, number of particles)
- Run this so that it writes an out file with starting particle positions and vectors
- From a terminal, set your desired parameters in "nearN_ordered.go": target volume fraction, particle growth parameters, step limits
- Build the .go file with go build
- Run ./nearN_ordered.go [FILENAME FOM MATHEMATICA] [Y/N], where Y says to write out contact and positional data at the end



![alt text](https://imgur.com/v3hV4RK)
