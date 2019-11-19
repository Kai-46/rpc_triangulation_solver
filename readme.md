# Triangulation Solver for the RPC Model in Satellite Images

## Introduction
An RPC model maps a 3D point (lat, lon, alt) to a pixel (col, row) in a satellite image. It is a non-linear mapping. 

For a specific 3D point, if you have multiple observations in different images (called a feature track), then you can triangulate the 3D point's coordinates with RPC models. This repo implements this functionality.

The triangulation solver basically operates in three steps:

* linearize each RPC model in a local area, or you can say 'locally approxmiate each RPC model with an affine model';
* solve a linear system to get an initial guess of the 3D point's coordinates;
* then, use [Ceres](http://ceres-solver.org/) to iteratively refine the initial guess by minimizing the average reprojection error.

## Installation

* install [Ceres](http://ceres-solver.org/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) on your machine
* compile the C++ backbone "multi_rpc_triangulate/" via:
```{r, engine='bash'}
cd multi_rpc_triangulate && mkdir build && cd build && cmake .. && make
```
* the python interface "triangulate.py" uses python3 instead of python2; it requires the minimal dependency: numpy

## Quick Start
The folder "example/" contains two example input configuration files: "metas.json" for specifying RPC camera parameters and image sizes, "tracks.txt" for specifying feature tracks. 

The format of "metas.json" is,

```{r, engine='bash'}
{
"img1": {"rpc": {"colNum": [], "colDen": [], "rowNum": [], "rowDen": [], 
                 "latOff": , "latScale": , "lonOff": , "lonScale": , "latOff": , "latScale": , 
                 "colOff": , "colScale": , "rowOff": , "rowScale": }, 
         "width": , "height": },
"img2": {"rpc": {"colNum": [], "colDen": [], "rowNum": [], "rowDen": [], 
                 "latOff": , "latScale": , "lonOff": , "lonScale": , "latOff": , "latScale": , 
                 "colOff": , "colScale": , "rowOff": , "rowScale": }, 
         "width": , "height": }
}
```



