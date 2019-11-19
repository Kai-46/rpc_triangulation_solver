#  ===============================================================================================================
#  Copyright (c) 2019, Cornell University. All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without modification, are permitted provided that
#  the following conditions are met:
#
#      * Redistributions of source code must retain the above copyright otice, this list of conditions and
#        the following disclaimer.
#
#      * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
#        the following disclaimer in the documentation and/or other materials provided with the distribution.
#
#      * Neither the name of Cornell University nor the names of its contributors may be used to endorse or
#        promote products derived from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
#  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
#  OF SUCH DAMAGE.
#
#  Author: Kai Zhang (kz298@cornell.edu)
#
#  The research is based upon work supported by the Office of the Director of National Intelligence (ODNI),
#  Intelligence Advanced Research Projects Activity (IARPA), via DOI/IBC Contract Number D17PC00287.
#  The U.S. Government is authorized to reproduce and distribute copies of this work for Governmental purposes.
#  ===============================================================================================================


# Copyright (C) 2015, Carlo de Franchis <carlo.de-franchis@cmla.ens-cachan.fr>
# Copyright (C) 2015, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>
# Copyright (C) 2015, Enric Meinhardt <enric.meinhardt@cmla.ens-cachan.fr>


import numpy as np

def apply_poly(poly, x, y, z):
    """
    Evaluates a 3-variables polynom of degree 3 on a triplet of numbers.
    Args:
        poly: list of the 20 coefficients of the 3-variate degree 3 polynom,
            ordered following the RPC convention.
        x, y, z: triplet of floats. They may be numpy arrays of same length.
    Returns:
        the value(s) of the polynom on the input point(s).
    """
    out = 0
    out += poly[0]
    out += poly[1]*y + poly[2]*x + poly[3]*z
    out += poly[4]*y*x + poly[5]*y*z +poly[6]*x*z
    out += poly[7]*y*y + poly[8]*x*x + poly[9]*z*z
    out += poly[10]*x*y*z
    out += poly[11]*y*y*y
    out += poly[12]*y*x*x + poly[13]*y*z*z + poly[14]*y*y*x
    out += poly[15]*x*x*x
    out += poly[16]*x*z*z + poly[17]*y*y*z + poly[18]*x*x*z
    out += poly[19]*z*z*z
    return out

def apply_rfm(num, den, x, y, z):
    """
    Evaluates a Rational Function Model (rfm), on a triplet of numbers.
    Args:
        num: list of the 20 coefficients of the numerator
        den: list of the 20 coefficients of the denominator
            All these coefficients are ordered following the RPC convention.
        x, y, z: triplet of floats. They may be numpy arrays of same length.
    Returns:
        the value(s) of the rfm on the input point(s).
    """
    return apply_poly(num, x, y, z) / apply_poly(den, x, y, z)

class RPCModel(object):
    def __init__(self, meta_dict):
        rpc_dict = meta_dict['rpc']
        # normalization constant
        self.colOff = rpc_dict['colOff']
        self.colScale = rpc_dict['colScale']

        self.rowOff = rpc_dict['rowOff']
        self.rowScale = rpc_dict['rowScale']

        self.latOff = rpc_dict['latOff']
        self.latScale = rpc_dict['latScale']

        self.lonOff = rpc_dict['lonOff']
        self.lonScale = rpc_dict['lonScale']

        self.altOff = rpc_dict['altOff']
        self.altScale = rpc_dict['altScale']

        # polynomial coefficients
        self.colNum = rpc_dict['colNum']
        self.colDen = rpc_dict['colDen']
        self.rowNum = rpc_dict['rowNum']
        self.rowDen = rpc_dict['rowDen']

        # img_width, img_height
        self.width = meta_dict['width']
        self.height = meta_dict['height']

    def projection(self, lat, lon, alt):
        cLon = (lon - self.lonOff) / self.lonScale
        cLat = (lat - self.latOff) / self.latScale
        cAlt = (alt - self.altOff) / self.altScale
        cCol = apply_rfm(self.colNum, self.colDen, cLat, cLon, cAlt)
        cRow = apply_rfm(self.rowNum, self.rowDen, cLat, cLon, cAlt)
        col = cCol*self.colScale + self.colOff
        row = cRow*self.rowScale + self.rowOff
        return col, row

    def to_string(self):
        string = ''

        tmp = [str(x) for x in self.colNum]
        string += ' '.join(tmp) + '\n'

        tmp = [str(x) for x in self.colDen]
        string += ' '.join(tmp) + '\n'

        tmp = [str(x) for x in self.rowNum]
        string += ' '.join(tmp) + '\n'

        tmp = [str(x) for x in self.rowDen]
        string += ' '.join(tmp) + '\n'

        string += '{} {} {} {} {} {} {} {} {} {}\n'.format(self.latOff, self.latScale, self.lonOff, self.lonScale, self.altOff, self.altScale,
                                                           self.colOff, self.colScale, self.rowOff, self.rowScale)
        return string


##################################################################################################
# local approximation of the rpc model with an affine model
##################################################################################################
def generate_grid(x_points, y_points, z_points):
    x_point_cnt = x_points.size
    y_point_cnt = y_points.size
    z_point_cnt = z_points.size
    point_cnt = x_point_cnt * y_point_cnt * z_point_cnt

    xx, yy = np.meshgrid(x_points, y_points, indexing='ij')
    xx = np.reshape(xx, (-1, 1))
    yy = np.reshape(yy, (-1, 1))
    xx = np.tile(xx, (z_point_cnt, 1))
    yy = np.tile(yy, (z_point_cnt, 1))

    zz = np.zeros((point_cnt, 1))
    for j in range(z_point_cnt):
        idx1 = j * x_point_cnt * y_point_cnt
        idx2 = (j + 1) * x_point_cnt * y_point_cnt
        zz[idx1:idx2, 0] = z_points[j]

    return xx, yy, zz

def solve_affine(xx, yy, zz, col, row, valid_mask=None):
    diff_size = np.array([yy.size - xx.size, zz.size - xx.size, col.size - xx.size, row.size - xx.size])
    assert (np.all(diff_size == 0))

    if valid_mask is not None:
        xx = xx[valid_mask].reshape((-1, 1))
        yy = yy[valid_mask].reshape((-1, 1))
        zz = zz[valid_mask].reshape((-1, 1))
        row = row[valid_mask].reshape((-1, 1))
        col = col[valid_mask].reshape((-1, 1))

    # construct a least square problem
    point_cnt = xx.size
    all_ones = np.ones((point_cnt, 1))
    all_zeros = np.zeros((point_cnt, 4))
    # construct the least square problem
    A1 = np.hstack((xx, yy, zz, all_ones, all_zeros))
    A2 = np.hstack((all_zeros, xx, yy, zz, all_ones))

    A = np.vstack((A1, A2))
    b = np.vstack((col, row))
    res = np.linalg.lstsq(A, b, rcond=-1)
    P = res[0].reshape((2, 4))

    return P

def affine_approx(rpc_model, bbx):
    lat_min, lat_max, lon_min, lon_max, alt_min, alt_max = bbx

    xy_axis_grid_points = 100
    z_axis_grid_points = 20

    lat_points = np.linspace(lat_min, lat_max, xy_axis_grid_points)
    lon_points = np.linspace(lon_min, lon_max, xy_axis_grid_points)
    alt_points = np.linspace(alt_min, alt_max, z_axis_grid_points)
    lat_points, lon_points, alt_points = generate_grid(lat_points, lon_points, alt_points)

    col, row = rpc_model.projection(lat_points, lon_points, alt_points)

    valid_mask = np.logical_and.reduce((col>=0, row>=0, col<rpc_model.width, row<rpc_model.height))

    P = solve_affine(lat_points, lon_points, alt_points, col, row, valid_mask)

    return P

if __name__ == '__main__':
    pass
