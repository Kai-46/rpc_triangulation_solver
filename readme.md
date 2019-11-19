# Triangulation Solver for the RPC Model in Satellite Images

**Project page**: [https://kai-46.github.io/VisSat/](https://kai-46.github.io/VisSat/)

## Introduction
An RPC model maps a 3D point (lat, lon, alt) to a pixel (col, row) in a satellite image. It is a non-linear mapping. 

For a specific 3D point, if you have multiple observations in different images (called a feature track), then you can triangulate the 3D point's coordinates with RPC models. This repo implements this functionality.

The triangulation solver basically operates in three steps:

1. linearize each RPC model in a local area, or you can say 'locally approxmiate each RPC model with an affine model';
2. solve a linear system to get an initial guess of the 3D point's coordinates;
3. then, use [Ceres](http://ceres-solver.org/) to iteratively refine the initial guess by minimizing the average reprojection error.

## Installation

* install [Ceres](http://ceres-solver.org/) and [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) on your machine
* compile the C++ backbone "multi_rpc_triangulate/" via:
```{r, engine='bash'}
cd multi_rpc_triangulate && mkdir build && cd build && cmake .. && make
```
* the python interface "triangulate.py" uses python3 instead of python2; it requires the minimal dependency: numpy

## Quick Start
The folder "example/" contains three example input configuration files: "metas.json" for specifying RPC camera parameters and image sizes, "tracks.txt" for specifying feature tracks, "bbx.json" for specifying the local area in which the RPC models will be linearized. You can triangulate this example by typing:
```{r, engine='bash'}
python3 triangulate.py
```
The solver will write the triangulation results to "example/results.txt". 

The format of the input file "metas.json" is,
```{r, engine='bash'}
{
"img1": {"rpc": {"colNum": [], "colDen": [], "rowNum": [], "rowDen": [], 
                 "latOff": , "latScale": , "lonOff": , "lonScale": , "altOff": , "altScale": , 
                 "colOff": , "colScale": , "rowOff": , "rowScale": }, 
         "width": , "height": },
"img2": {"rpc": {"colNum": [], "colDen": [], "rowNum": [], "rowDen": [], 
                 "latOff": , "latScale": , "lonOff": , "lonScale": , "altOff": , "altScale": , 
                 "colOff": , "colScale": , "rowOff": , "rowScale": }, 
         "width": , "height": }
...
}
```
The format of the input file "bbx.json" is,
```{r, engine='bash'}
{"lat_min": , "lat_max": , "lon_min": , "lon_max": , "alt_min": , "alt_max": }
```
The format of the input file "tracks.txt" is,
```{r, engine='bash'}
number_of_tracks
track_length img1 col row img2 col row img3 col row ...
track_length img1 col row img2 col row ...
...
```
The format of the output file "results.txt" is, 
```{r, engine='bash'}
number_of_3D_points
initial_lat initial_lon initial_alt inital_reproj_err final_lat final_lon final_alt final_reproj_err
initial_lat initial_lon initial_alt inital_reproj_err final_lat final_lon final_alt final_reproj_err
...
```

## License
This software uses the [3-clause BSD license](https://opensource.org/licenses/BSD-3-Clause).
