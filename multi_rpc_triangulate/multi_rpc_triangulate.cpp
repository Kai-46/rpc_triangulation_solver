// ===============================================================================================================
// Copyright (c) 2019, Cornell University. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that
// the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright otice, this list of conditions and
//       the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
//       the following disclaimer in the documentation and/or other materials provided with the distribution.
//
//     * Neither the name of Cornell University nor the names of its contributors may be used to endorse or
//       promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.
//
// Author: Kai Zhang (kz298@cornell.edu)
//
// The research is based upon work supported by the Office of the Director of National Intelligence (ODNI),
// Intelligence Advanced Research Projects Activity (IARPA), via DOI/IBC Contract Number D17PC00287.
// The U.S. Government is authorized to reproduce and distribute copies of this work for Governmental purposes.
// ===============================================================================================================


#include "ceres/ceres.h"
#include "glog/logging.h"
#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <limits>

using Eigen::MatrixXd;

using ceres::AutoDiffCostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

using std::cout;
using std::vector;
using std::string;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::stod;
using std::runtime_error;
using std::setprecision;

typedef std::numeric_limits<double> dbl;
#define INIT_VAL -1e10


struct RPCCamera {
    RPCCamera() {}
    
    // RPC camera parameters
    double col_numera[20] = {INIT_VAL};
    double col_denomi[20] = {INIT_VAL};
    double row_numera[20] = {INIT_VAL};
    double row_denomi[20] = {INIT_VAL};
    double lat_off = INIT_VAL, lat_scale = INIT_VAL;
    double lon_off = INIT_VAL, lon_scale = INIT_VAL;
    double alt_off = INIT_VAL, alt_scale = INIT_VAL;
    double row_off = INIT_VAL, row_scale = INIT_VAL;
    double col_off = INIT_VAL, col_scale = INIT_VAL;

    // affine approximation of the RPC camera
    // [col, row]^T = M * [lat, lon, alt]^T
    double M11;
    double M12;
    double M13;
    double M14;

    double M21;
    double M22;
    double M23;
    double M24;
};

struct Observation {
    Observation(RPCCamera* cam, double col, double row): cam(cam), col(col), row(row) {}
    
    RPCCamera* cam = NULL;
    double col = INIT_VAL;
    double row = INIT_VAL;
};

struct ReprojResidual {
    ReprojResidual(Observation* pixel): pixel(pixel) {}
    
    template <typename T>
    bool operator() (const T* const lat, const T* const lon, const T* const alt,
                     T* residuals) const {
        RPCCamera& cam = *(this->pixel->cam);
        
        T lat_normed = (lat[0] - T(cam.lat_off)) / T(cam.lat_scale);
        T lon_normed = (lon[0] - T(cam.lon_off)) / T(cam.lon_scale);
        T alt_normed = (alt[0] - T(cam.alt_off)) / T(cam.alt_scale);
        
        T row_numera = this->apply_poly(cam.row_numera, lat_normed, lon_normed, alt_normed);
        T row_denomi = this->apply_poly(cam.row_denomi, lat_normed, lon_normed, alt_normed);
        
        T predict_row = row_numera / row_denomi * T(cam.row_scale) + T(cam.row_off);
        
        T col_numera = this->apply_poly(cam.col_numera, lat_normed, lon_normed, alt_normed);
        T col_denomi = this->apply_poly(cam.col_denomi, lat_normed, lon_normed, alt_normed);
        
        T predict_col = col_numera / col_denomi * T(cam.col_scale) + T(cam.col_off);
        
        residuals[0] = predict_row - T(this->pixel->row);
        residuals[1] = predict_col - T(this->pixel->col);
        return true;
    }
    
private:
    template <typename T>
    T apply_poly(const double* const poly, T x, T y, T z) const {
        T out = T(poly[0]);
        out += poly[1]*y + poly[2]*x + poly[3]*z;
        out += poly[4]*y*x + poly[5]*y*z +poly[6]*x*z;
        out += poly[7]*y*y + poly[8]*x*x + poly[9]*z*z;
        out += poly[10]*x*y*z;
        out += poly[11]*y*y*y;
        out += poly[12]*y*x*x + poly[13]*y*z*z + poly[14]*y*y*x;
        out += poly[15]*x*x*x;
        out += poly[16]*x*z*z + poly[17]*y*y*z + poly[18]*x*x*z;
        out += poly[19]*z*z*z;
        
        return out;
    }
    
private:
    Observation* pixel = NULL;
};

void solve_initial(const vector<Observation*>& pixels, vector<double>& initial) {
    assert (pixels.size() >= 2 && initial.size() == 3);

    MatrixXd A(2 * pixels.size(), 3);  // one observation contributes to two equations
    MatrixXd b(2 * pixels.size(), 1);  
    for (int i=0; i < pixels.size(); ++i) {
        A(i * 2, 0) = pixels[i]->cam->M11;
        A(i * 2, 1) = pixels[i]->cam->M12;
        A(i * 2, 2) = pixels[i]->cam->M13;
        b(i * 2, 0) = pixels[i]->col - pixels[i]->cam->M14;

        A(i * 2 + 1, 0) = pixels[i]->cam->M21;
        A(i * 2 + 1, 1) = pixels[i]->cam->M22;
        A(i * 2 + 1, 2) = pixels[i]->cam->M23;
        b(i * 2 + 1, 0) = pixels[i]->row - pixels[i]->cam->M24;
    }

    MatrixXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    for (int i=0; i < 3; ++i) {
        initial[i] = x(i, 0);
    }
}

void refine_initial(const vector<Observation*>& pixels, const vector<double>& initial, vector <double>& final, vector <double>& reproj_error) {
    assert (initial.size() == 3 && final.size() == 3 && reproj_error.size() == 2);

    double lat = initial[0];
    double lon = initial[1];
    double alt = initial[2];
    
    Problem problem;
    for (int i=0; i < pixels.size(); ++i) {
        problem.AddResidualBlock(
                                 new AutoDiffCostFunction<ReprojResidual, 2, 1, 1, 1>(new ReprojResidual(pixels[i])),
                                 NULL, &lat, &lon, &alt);
    }
    
    Solver::Options options;
    options.max_num_iterations = 100;
    // options.function_tolerance = 1e-10;
    options.linear_solver_type = ceres::DENSE_QR;
    // set the following options to true for debugging
    options.minimizer_progress_to_stdout = false;
    
    Solver::Summary summary;
    Solve(options, &problem, &summary);

    double init_error = sqrt(summary.initial_cost * 2 / pixels.size());
    double final_error = sqrt(summary.final_cost * 2 / pixels.size());

    // for debugging
    // cout << summary.BriefReport() << "\n";
    // cout << "\ninitial Point: (" << initial[0] << "," << initial[1] << "," << initial[2] << "), reproj_error: " << init_error << " pixels\n";
    // cout << "final Point:  (" << lat << "," << lon << "," << alt << "), reproj_error: " << final_error << " pixels\n";

    // output
    final[0] = lat;
    final[1] = lon;
    final[2] = alt;
    reproj_error[0] = init_error;
    reproj_error[1] = final_error;
}

void read_rpc_cameras(const string& fname, vector<RPCCamera*>& rpc_cameras) {
    // number of rpc cameras
    // camera_id
    // 20 column numerator coefficients
    // 20 column denominator coefficients
    // 20 row numerator coefficients
    // 20 row denominator coefficients
    // 10 normalization constants: lat off, lat scale, lon off, lon scale, alt off, alt scale, col off, col scale
    // 8 affine approximation coefficients: M11, M12, M13, M14, M21, M22, M23, M24
    // ...

    ifstream infile;
    infile.open(fname);
    if (!infile) {
        throw runtime_error("unable to open " + fname);
    }

    int cnt; 
    infile >> cnt; 
    for (int i=0; i<cnt; ++i) {
        int cam_id;
        // read camera id
        infile >> cam_id;
        assert(cam_id == i);

        rpc_cameras.push_back(new RPCCamera());
        RPCCamera *cam = rpc_cameras.back();

        // read rpc camera parameters
        for (int i = 0; i < 20; ++i) {
            infile >> cam->col_numera[i];
        }
        for (int i = 0; i < 20; ++i) {
            infile >> cam->col_denomi[i];
        }
        for (int i = 0; i < 20; ++i) {
            infile >> cam->row_numera[i];
        }
        for (int i = 0; i < 20; ++i) {
            infile >> cam->row_denomi[i];
        }
        infile >> cam->lat_off >> cam->lat_scale;
        infile >> cam->lon_off >> cam->lon_scale;
        infile >> cam->alt_off >> cam->alt_scale;
        infile >> cam->col_off >> cam->col_scale;
        infile >> cam->row_off >> cam->row_scale;

        // read affine approximation parameters
        infile >> cam->M11 >> cam->M12 >> cam->M13 >> cam->M14;
        infile >> cam->M21 >> cam->M22 >> cam->M23 >> cam->M24;
    }

    infile.close();
}

void triangulate_tracks(const string& cameras_fname, const string& tracks_fname, const string& results_fname) {
    vector<RPCCamera*> rpc_cameras;
    read_rpc_cameras(cameras_fname, rpc_cameras);

    ifstream infile;
    infile.open(tracks_fname);
    if (!infile) {
        throw runtime_error("unable to open " + tracks_fname);
    }

    ofstream outfile;
    outfile.open(results_fname);
    if (!outfile) {
        throw runtime_error("unable to open " + results_fname);
    }
    outfile << setprecision(dbl::max_digits10);

    int cnt; 
    infile >> cnt; 
    outfile << cnt << '\n';
    for (int i=0; i<cnt; ++i) {
        // read feature track length
        int len;
        infile >> len;
        // read all the observations for this track
        vector<Observation*> pixels;
        for (int j=0; j < len; ++j) {
            int cam_id;
            double col, row;
            infile >> cam_id >> col >> row;

            assert(cam_id < rpc_cameras.size());
            pixels.push_back(new Observation(rpc_cameras[cam_id], col, row));
        }

        // solve for (lat, lon, alt)
        vector<double> initial(3, INIT_VAL);
        vector<double> final(3, INIT_VAL);
        vector<double> reproj_error(2, INIT_VAL);
        solve_initial(pixels, initial);
        refine_initial(pixels, initial, final, reproj_error);

        // write results to file
        // each line is "intial lat, initial lon, initial alt, initial reproj err, final lat, final lon, final alt, final reproj err"
        outfile << initial[0] << " " << initial[1] << " " << initial[2] << " " << reproj_error[0] << " ";
        outfile << final[0] << " " << final[1] << " " << final[2] << " " << reproj_error[1] << "\n";

        // free memory
        for (int i = 0; i < pixels.size(); ++i) {
            delete pixels[i];
        }
    }

    // close file
    infile.close();
    outfile.close();

    // free memory
    for (int i = 0; i < rpc_cameras.size(); ++i) {
        delete rpc_cameras[i];
    }
}


int main(int argc, char** argv) {
    // program name, cameras_fname, tracks_fname, results_fname
    assert(argc == 4);
    google::InitGoogleLogging(argv[0]);
    string cameras_fname = string(argv[1]);
    string tracks_fname = string(argv[2]);
    string results_fname = string(argv[3]);

    triangulate_tracks(cameras_fname, tracks_fname, results_fname);
    return 0;
}
