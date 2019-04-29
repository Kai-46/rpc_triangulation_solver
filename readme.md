## Triangulation solver for the RPC model used in satellite images

The triangulation solver requires the local linearization of RPC model as input. It basically operates in two steps.

* solve a linear system to get a initial guess of the 3D point coordinates
* then, use ceres to iteratively refine the initial guess by solving a non-linear least square problem
