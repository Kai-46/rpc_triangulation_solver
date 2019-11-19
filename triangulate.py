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


import os
import json
import numpy as np
import json
from rpc_model import RPCModel, affine_approx
import tempfile

# meta_file is a json
# bbx_file is a json
# track_file is a text file
def triangulate(meta_file, bbx_file, track_file, out_file):
    with open(meta_file) as fp:
        metas = json.load(fp)

    img_names = sorted(metas.keys())
    img_name2id = dict([(img_names[i], i) for i in range(len(img_names))])

    with open(bbx_file) as fp:
        bbx = json.load(fp)
        bbx = (bbx['lat_min'], bbx['lat_max'], bbx['lon_min'], bbx['lon_max'], bbx['alt_min'], bbx['alt_max'])

    rpc_cameras = []
    affine_cameras = []
    for name in img_names:
        rpc_cam = RPCModel(metas[name])
        rpc_cameras.append(rpc_cam)

        # apprximate the rpc camera with an affine one in a local area
        affine_cam = affine_approx(rpc_cam, bbx)
        affine_cameras.append(list(affine_cam.flatten()))

    tmpdir = tempfile.TemporaryDirectory()
    # write rpc cameras
    camera_fname = os.path.join(tmpdir.name, 'rpc_cameras.txt')
    with open(camera_fname, 'w') as fp:
        fp.write('{}\n'.format(len(rpc_cameras)))
        for i in range(len(rpc_cameras)):
            fp.write('{}\n'.format(i))
            string = rpc_cameras[i].to_string()
            string += ' '.join([str(x) for x in affine_cameras[i]])
            fp.write(string + '\n')

    # write feature tracks
    lines_to_write = []
    with open(track_file) as fp_in:
        track_cnt = int(fp_in.readline().strip())
        lines_to_write.append('{}\n'.format(track_cnt))
        for i in range(track_cnt):
            line = fp_in.readline().strip().split(' ')
            if len(line) == 0:
                continue
            assert ((len(line) - 1) % 3 == 0)
            pixels_cnt = int(line[0])

            string = '{}'.format(pixels_cnt)
            line = line[1:]
            for j in range(pixels_cnt):
                name = line[3 * j].strip()
                col = line[3 * j + 1].strip()
                row = line[3 * j + 2].strip()
                string += ' {} {} {}'.format(img_name2id[name], col, row)
            lines_to_write.append(string + '\n')

    track_fname = os.path.join(tmpdir.name, 'feature_tracks.txt')
    with open(track_fname, 'w') as fp:
        fp.writelines(lines_to_write)

    # now launch the C++ program
    cmd = './multi_rpc_triangulate/build/multi_rpc_triangulate {} {} {}'.format(camera_fname, track_fname, out_file)
    os.system(cmd)

    tmpdir.cleanup()

if __name__ == '__main__':
    meta_file = './example/metas.json'
    bbx_file = './example/bbx.json'
    track_file = './example/tracks.txt'
    out_file = './example/results.txt'
    triangulate(meta_file, bbx_file, track_file, out_file)