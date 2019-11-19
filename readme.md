# Triangulation Solver for the RPC Model in Satellite Images

## Introduction
An RPC model maps a 3D point (lat, lon, alt) to a pixel (col, row) in a satellite image. It is a non-linear mapping. 

For a specific 3D point, if you have multiple observations in different images, then you can triangulate the 3D point's coordinates with RPC models. This repo implements this functionality.

The triangulation solver basically operates in three steps:

* linearize each RPC model in a local area, or you can say 'locally approxmiate each RPC model with an affine model';
* solve a linear system to get an initial guess of the 3D point's coordinates;
* then, use [ceres](http://ceres-solver.org/) to iteratively refine the initial guess by minimizing the average reprojection error.

## Installation





