// ================================================================================================
//  Copyright (C) 2024, SCHOOL OF GEODESY AND GEOMATICS, WUHAN UNIVERSITY, CHINA
//  All Rights Reserved.
//
//  Authors: Rongling Zhang, Pengcheng Wei et al.
//  Contact:hitzrl@foxmail.com
//  See LICENSE file for more detailed license information.
// ================================================================================================
// Thanks to the work of Liu, et al:
// https://github.com/hku-mars/BALM/tree/master/src/benchmark

#ifndef BAVOXEL_HPP
#define BAVOXEL_HPP

#include "tools.hpp"
#include <Eigen/Eigenvalues>
#include <thread>
#include <chrono>
#include <algorithm>
#include <pcl/common/transforms.h>
#include <iostream>

#include <unsupported/Eigen/MatrixFunctions>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation//sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/io/pcd_io.h>
#include<pcl/registration/correspondence_estimation.h>
#include <sophus/se3.hpp>
#include "AndersonAcceleration.h"
#define SAME_THRESHOLD 1e-6
static int layer_limit = 3;
static int layer_size[] = {8, 8, 8, 8};
static float eigen_value_array[] = {1.0/4.0, 1.0/4.0, 1.0/4.0,1.0/4.0};
static int min_ps = 15;
static double one_three = (1.0 / 3.0);

static double voxel_size;
static int life_span = 1000;
static int win_size;

static int merge_enable = 1;

class VOX_HESS
{
public:
    vector<const PointCluster*> sig_vecs;
    vector<const vector<PointCluster>*> plvec_voxels;
    vector<double> coeffs, coeffs_back;


    vector<pcl::PointCloud<PointType>::Ptr> plptrs;

    bool push_voxel(const vector<PointCluster> *vec_orig, const PointCluster *fix, double feat_eigen, int layer)
    {
        int process_size = 0;
        for(int i=0; i<win_size; i++)
            if((*vec_orig)[i].N >=15)
                process_size++;


        if(process_size < 2) {
            return false;
        };

        double coe = 1 - feat_eigen/eigen_value_array[layer];
        coe = coe * coe;
        coe = 0;
        for(int j=0; j<win_size; j++)
            coe += (*vec_orig)[j].N;

        plvec_voxels.push_back(vec_orig);
        sig_vecs.push_back(fix);
        coeffs.push_back(coe);
        pcl::PointCloud<PointType>::Ptr plptr(new pcl::PointCloud<PointType>());
        plptrs.push_back(plptr);
        return true;
    }

    void acc_evaluate2(const vector<IMUST> &xs, int head, int end, Eigen::MatrixXd &Hess, Eigen::VectorXd &JacT, double &residual)
    {
        Hess.setZero(); JacT.setZero(); residual = 0;
        vector<PointCluster> sig_tran(win_size);
        const int kk = 0;

        PLV(3) viRiTuk(win_size);
        PLM(3) viRiTukukT(win_size);

        vector<Eigen::Matrix<double, 3, 6>, Eigen::aligned_allocator<Eigen::Matrix<double, 3, 6>>> Auk(win_size);
        Eigen::Matrix3d umumT;

        for(int a=head; a<end; a++)
        {
            const vector<PointCluster> &sig_orig = *plvec_voxels[a];
            double coe = coeffs[a];

            PointCluster sig = *sig_vecs[a];
            for(int i=0; i<win_size; i++)
                if(sig_orig[i].N != 0)
                {
                    sig_tran[i].transform(sig_orig[i], xs[i]);
                    sig += sig_tran[i];
                }

            const Eigen::Vector3d &vBar = sig.v / sig.N;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(sig.P/sig.N - vBar * vBar.transpose());
            const Eigen::Vector3d &lmbd = saes.eigenvalues();
            const Eigen::Matrix3d &U = saes.eigenvectors();
            int NN = sig.N;

            Eigen::Vector3d u[3] = {U.col(0), U.col(1), U.col(2)};

            const Eigen::Vector3d &uk = u[kk];
            Eigen::Matrix3d ukukT = uk * uk.transpose();
            umumT.setZero();
            for(int i=0; i<3; i++)
                if(i !=  kk)
                    umumT += 2.0/(lmbd[kk] - lmbd[i]) * u[i] * u[i].transpose();

            for(int i=0; i<win_size; i++)
                if(sig_orig[i].N != 0)
                {
                    Eigen::Matrix3d Pi = sig_orig[i].P;
                    Eigen::Vector3d vi = sig_orig[i].v;
                    Eigen::Matrix3d Ri = xs[i].R;
                    double ni = sig_orig[i].N;

                    Eigen::Matrix3d vihat; vihat << SKEW_SYM_MATRX(vi);
                    Eigen::Vector3d RiTuk = Ri.transpose() * uk;
                    Eigen::Matrix3d RiTukhat; RiTukhat << SKEW_SYM_MATRX(RiTuk);

                    Eigen::Vector3d PiRiTuk = Pi * RiTuk;
                    viRiTuk[i] = vihat * RiTuk;
                    viRiTukukT[i] = viRiTuk[i] * uk.transpose();

                    Eigen::Vector3d ti_v = xs[i].p - vBar;
                    double ukTti_v = uk.dot(ti_v);

                    Eigen::Matrix3d combo1 = hat(PiRiTuk) + vihat * ukTti_v;
                    Eigen::Vector3d combo2 = Ri*vi + ni*ti_v;
                    Auk[i].block<3, 3>(0, 0) = (Ri*Pi + ti_v*vi.transpose()) * RiTukhat - Ri*combo1;
                    Auk[i].block<3, 3>(0, 3) = combo2 * uk.transpose() + combo2.dot(uk) * I33;
                    Auk[i] /= NN;

                    const Eigen::Matrix<double, 6, 1> &jjt = Auk[i].transpose() * uk;
                    JacT.block<6, 1>(6*i, 0) += coe * jjt;


                    const Eigen::Matrix3d &HRt = 2.0/NN * (1.0-ni/NN) * viRiTukukT[i];
                    Eigen::Matrix<double, 6, 6> Hb = Auk[i].transpose() * umumT * Auk[i];
                    Hb.block<3, 3>(0, 0) += 2.0/NN * (combo1 - RiTukhat*Pi) * RiTukhat - 2.0/NN/NN * viRiTuk[i] * viRiTuk[i].transpose() - 0.5*hat(jjt.block<3, 1>(0, 0));
                    Hb.block<3, 3>(0, 3) += HRt;
                    Hb.block<3, 3>(3, 0) += HRt.transpose();
                    Hb.block<3, 3>(3, 3) += 2.0/NN * (ni - ni*ni/NN) * ukukT;

                    Hess.block<6, 6>(6*i, 6*i) += coe * Hb;
                }

            for(int i=0; i<win_size-1; i++)
                if(sig_orig[i].N != 0)
                {
                    double ni = sig_orig[i].N;
                    for(int j=i+1; j<win_size; j++)
                        if(sig_orig[j].N != 0)
                        {
                            double nj = sig_orig[j].N;
                            Eigen::Matrix<double, 6, 6> Hb = Auk[i].transpose() * umumT * Auk[j];
                            Hb.block<3, 3>(0, 0) += -2.0/NN/NN * viRiTuk[i] * viRiTuk[j].transpose();
                            Hb.block<3, 3>(0, 3) += -2.0*nj/NN/NN * viRiTukukT[i];
                            Hb.block<3, 3>(3, 0) += -2.0*ni/NN/NN * viRiTukukT[j].transpose();
                            Hb.block<3, 3>(3, 3) += -2.0*ni*nj/NN/NN * ukukT;

                            Hess.block<6, 6>(6*i, 6*j) += coe * Hb;
                        }
                }

            residual += coe * lmbd[kk];
        }

        for(int i=1; i<win_size; i++)
            for(int j=0; j<i; j++)
                Hess.block<6, 6>(6*i, 6*j) = Hess.block<6, 6>(6*j, 6*i).transpose();

    }

    void left_evaluate(const vector<IMUST> &xs, int head, int end, Eigen::MatrixXd &Hess, Eigen::VectorXd &JacT, double &residual)
    {
        Hess.setZero(); JacT.setZero(); residual = 0;
        int l = 0;
        Eigen::Matrix<double, 3, 4> Sp;
        Sp.setZero(); Sp.block<3, 3>(0, 0).setIdentity();
        Eigen::Matrix4d F; F.setZero(); F(3, 3) = 1;

        PLM(4) T(win_size);
        for(int i=0; i<win_size; i++)
            T[i] << xs[i].R, xs[i].p, 0, 0, 0, 1;

        vector<PLM(4)*> Cs;
        for(int a=0; a<plvec_voxels.size(); a++)
        {
            const vector<PointCluster> &sig_orig = *plvec_voxels[a];
            PLM(4) *Co = new PLM(4)(win_size, Eigen::Matrix4d::Zero());
            for(int i=0; i<win_size; i++)
                Co->at(i) << sig_orig[i].P, sig_orig[i].v, sig_orig[i].v.transpose(), sig_orig[i].N;
            Cs.push_back(Co);
        }


        for(int a=head; a<end; a++)
        {

            double coe = coeffs[a];


            Eigen::Matrix4d C; C.setZero();


            PLM(4) &Co = *Cs[a];
            for(int i=0; i<win_size; i++)
                if((int)Co[i](3, 3) > 0)
                    C += T[i] * Co[i] * T[i].transpose();

            double NN = C(3, 3);
            C = C / NN;



            Eigen::Vector3d v_bar = C.block<3, 1>(0, 3);

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(C.block<3, 3>(0, 0) - v_bar * v_bar.transpose() );
            Eigen::Vector3d lmbd = saes.eigenvalues();
            Eigen::Matrix3d Uev = saes.eigenvectors();

            residual += coe * lmbd[l];

            Eigen::Vector3d u[3] = {Uev.col(0), Uev.col(1), Uev.col(2)};
            Eigen::Matrix<double, 4, 6> U[3];

            PLV(-1) g_kl(3);
            for(int k=0; k<3; k++)
            {
                g_kl[k].resize(6*win_size); g_kl[k].setZero();
                U[k].setZero();
                U[k].block<3, 3>(0, 0) = hat(u[k]);
                U[k].block<1, 3>(3, 3) = u[k];
            }

            for(int j=0; j<win_size; j++)
                for(int k=0; k<3; k++)
                    if(Co[j](3, 3) > 0.1)
                    {
                        Eigen::Matrix<double, 3, 4> SpTC = Sp * (T[j] - C*F) * Co[j] * T[j].transpose();
                        Eigen::Matrix<double, 1, 6> g1, g2;
                        g1 = u[l].transpose() * SpTC * U[k];
                        g2 = u[k].transpose() * SpTC * U[l];

                        g_kl[k].block<6, 1>(6*j, 0) = (g1 + g2).transpose() / NN;
                    }

            JacT += coe * g_kl[l];

            for(int i=0; i<win_size; i++)
                if(Co[i](3, 3) > 0.1)
                {
                    for(int j=0; j<win_size; j++)
                        if(Co[j](3, 3) > 0.1)
                        {
                            Eigen::Matrix4d Dij = Co[i] * F * Co[j];
                            Eigen::Matrix<double, 6, 6> Hs = -2.0/NN/NN * U[l].transpose() * T[i] * Dij * T[j].transpose() * U[l];

                            if(i == j)
                            {
                                Hs += 2/NN * U[l].transpose() * T[j] * Co[j] * T[j].transpose() * U[l];
                                Eigen::Vector3d SpTC = Sp * T[j] * Co[j] * (T[j] - C*F).transpose() * Sp.transpose() * u[l];
                                Eigen::Matrix3d h1 = hat(SpTC);
                                Eigen::Matrix3d h2 = hat(u[l]);

                                Hs.block<3, 3>(0, 0) += (h1*h2 + h2*h1) / NN;
                            }

                            Hess.block<6, 6>(6*i, 6*j) += coe * Hs;
                        }
                }

            for(int k=0; k<3; k++)
                if(k != l)
                    Hess += coe * 2.0/(lmbd[l] - lmbd[k]) * g_kl[k] * g_kl[k].transpose();

        }


    }

    void left_evaluate_acc2(const vector<IMUST> &xs, int head, int end, Eigen::MatrixXd &Hess, Eigen::VectorXd &JacT, double &residual)
    {
        bool isHessian= true;
        Hess.setZero(); JacT.setZero(); residual = 0;
        int l = 0;
        PLM(4) T(win_size);
        for(int i=0; i<win_size; i++)
            T[i] << xs[i].R, xs[i].p, 0, 0, 0, 1;

        vector<PLM(4)*> Cs;
        for(int a=0; a<plvec_voxels.size(); a++)
        {
            const vector<PointCluster> &sig_orig = *plvec_voxels[a];
            PLM(4) *Co = new PLM(4)(win_size, Eigen::Matrix4d::Zero());
            for(int i=0; i<win_size; i++)
                Co->at(i) << sig_orig[i].P, sig_orig[i].v, sig_orig[i].v.transpose(), sig_orig[i].N;
            Cs.push_back(Co);
        }

        for(int a=head; a<end; a++)
        {

            double coe = coeffs[a];
            Eigen::Matrix4d C; C.setZero();

            vector<int> Ns(win_size);

            PLM(4) &Co = *Cs[a];
            PLM(4) TC(win_size), TCT(win_size);
            for(int j=0; j<win_size; j++)
                if((int)Co[j](3, 3) > 0)
                {
                    TC[j] = T[j] * Co[j];
                    TCT[j] = TC[j] * T[j].transpose();
                    C += TCT[j];

                    Ns[j] = Co[j](3, 3);
                }

            double NN = C(3, 3);
            C = C / NN;
            Eigen::Vector3d v_bar = C.block<3, 1>(0, 3);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(C.block<3, 3>(0, 0) - v_bar * v_bar.transpose());
            Eigen::Vector3d lmbd = saes.eigenvalues();
            Eigen::Matrix3d Uev = saes.eigenvectors();


            residual += coe * lmbd[l];

            Eigen::Vector3d u[3] = {Uev.col(0), Uev.col(1), Uev.col(2)};
            Eigen::Matrix<double, 6, 4> U[3];
            PLV(6) g_kl[3];
            for(int k=0; k<3; k++)
            {
                g_kl[k].resize(win_size);
                U[k].setZero();
                U[k].block<3, 3>(0, 0) = hat(-u[k]);
                U[k].block<3, 1>(3, 3) = u[k];
            }

            PLV(6) UlTCF(win_size, Eigen::Matrix<double, 6, 1>::Zero());

            Eigen::VectorXd JacT_iter(6*win_size);
            for(int i=0; i<win_size; i++)
                if(Ns[i] != 0)
                {
                    Eigen::Matrix<double, 3, 4> temp = T[i].block<3, 4>(0, 0);
                    temp.block<3, 1>(0, 3) -= v_bar;
                    Eigen::Matrix<double, 4, 3> TC_TCFSp = TC[i] * temp.transpose();
                    for(int k=0; k<3; k++)
                    {
                        Eigen::Matrix<double, 6, 1> g1, g2;
                        g1 = U[k] * TC_TCFSp * u[l];
                        g2 = U[l] * TC_TCFSp * u[k];

                        g_kl[k][i] = (g1 + g2) / NN;
                    }

                    UlTCF[i] = (U[l] * TC[i]).block<6, 1>(0, 3);
                    JacT.block<6, 1>(6*i, 0) += coe * g_kl[l][i];


                    if(isHessian) {
                        Eigen::Matrix<double, 6, 6> Ha(-2.0 / NN / NN * UlTCF[i] * UlTCF[i].transpose());

                        Eigen::Matrix3d Ell = 1.0 / NN * hat(TC_TCFSp.block<3, 3>(0, 0) * u[l]) * hat(u[l]);
                        Ha.block<3, 3>(0, 0) += Ell + Ell.transpose();

                        for (int k = 0; k < 3; k++)
                            if (k != l)
                                Ha += 2.0 / (lmbd[l] - lmbd[k]) * g_kl[k][i] * g_kl[k][i].transpose();

                        Hess.block<6, 6>(6 * i, 6 * i) += coe * Ha;
                    }
                }
            if(isHessian) {
                for (int i = 0; i < win_size; i++)
                    if (Ns[i] != 0) {
                        Eigen::Matrix<double, 6, 6> Hb = U[l] * TCT[i] * U[l].transpose();
                        Hess.block<6, 6>(6 * i, 6 * i) += 2.0 / NN * coe * Hb;
                    }

                for (int i = 0; i < win_size - 1; i++)
                    if (Ns[i] != 0) {
                        for (int j = i + 1; j < win_size; j++)
                            if (Ns[j] != 0) {
                                Eigen::Matrix<double, 6, 6> Ha = -2.0 / NN / NN * UlTCF[i] * UlTCF[j].transpose();

                                for (int k = 0; k < 3; k++)
                                    if (k != l)
                                        Ha += 2.0 / (lmbd[l] - lmbd[k]) * g_kl[k][i] * g_kl[k][j].transpose();

                                Hess.block<6, 6>(6 * i, 6 * j) += coe * Ha;
                            }
                    }
            }
        }
        if(isHessian) {
            for (int i = 1; i < win_size; i++)
                for (int j = 0; j < i; j++)
                    Hess.block<6, 6>(6 * i, 6 * j) = Hess.block<6, 6>(6 * j, 6 * i).transpose();
        }

    }

    void evaluate_only_residual(const vector<IMUST> &xs, double &residual)
    {
        residual = 0;
        vector<PointCluster> sig_tran(win_size);
        int kk = 0; // The kk-th lambda value

        int gps_size = plvec_voxels.size();

        vector<double> ress(gps_size);

        for(int a=0; a<gps_size; a++)
        {
            const vector<PointCluster> &sig_orig = *plvec_voxels[a];
            PointCluster sig = *sig_vecs[a];

            for(int i=0; i<win_size; i++)
            {
                sig_tran[i].transform(sig_orig[i], xs[i]);
                sig += sig_tran[i];
            }

            Eigen::Vector3d vBar = sig.v / sig.N;
            Eigen::Matrix3d cmt = sig.P/sig.N - vBar * vBar.transpose();

            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(cmt);
            Eigen::Vector3d lmbd = saes.eigenvalues();

            residual += coeffs[a] * lmbd[kk];

            ress[a] = lmbd[kk];

        }

    }

    ~VOX_HESS()
    {
        int vsize = sig_vecs.size();
    }

};

class VOXEL_MERGE
{
public:
    vector<const PointCluster*> sig_vecs;
    vector<const vector<PointCluster>*> plvec_voxels;

    PLV(3) centers, directs, evalues;
    vector<pcl::PointCloud<PointType>::Ptr> plptrs;

    void push_voxel(const vector<PointCluster> *vec_orig, const PointCluster *fix, Eigen::Vector3d &center, Eigen::Vector3d &direct, Eigen::Vector3d &evalue, pcl::PointCloud<PointType>::Ptr plptr = nullptr)
    {
        int process_size = 0;
        for(int i=0; i<win_size; i++)
            if((*vec_orig)[i].N != 0)
                process_size++;

        if(process_size < 2) return;

        plvec_voxels.push_back(vec_orig);
        sig_vecs.push_back(fix);
        centers.push_back(center);
        directs.push_back(direct);
        evalues.push_back(evalue);
        plptrs.push_back(plptr);
    }

    void reorganize(VOX_HESS &voxhess, pcl::PointCloud<PointType> &pl_send, pcl::PointCloud<PointType> &pl_cent, vector<IMUST> &x_buf)
    {
        static double cos1 = cos(8/57.3);
        static double cos2 = cos(80/57.3);

        int vsize = centers.size();
        if(vsize <= 0) return;

        vector<vector<int>> groups;
        groups.push_back(vector<int>());
        groups[0].push_back(0);
        for(int i=1; i<vsize; i++)
        {
            Eigen::Vector3d c2 = centers[i];
            Eigen::Vector3d direct2 = directs[i];

            bool match = false;
            if(merge_enable)
            {
                int gsize = groups.size();
                for(int j=0; j<gsize; j++)
                {
                    int surf1 = groups[j][0];

                    Eigen::Vector3d c2c = c2 - centers[surf1];
                    double c2cd = c2c.norm();
                    c2c /= c2cd;
                    Eigen::Vector3d direct1 = directs[surf1];

                    double dot1 = fabs(direct1.dot(direct2));
                    double dot2 = fabs(c2c.dot(direct1));
                    double dot3 = fabs(c2c.dot(direct2));

                    bool c2flag = (dot2<cos2 && dot3<cos2) || (c2cd < 0.1);
                    if(dot1>cos1 && c2flag)
                    {
                        groups[j].push_back(i);
                        match = true;
                        break;
                    }
                }
            }

            if(!match)
            {
                groups.push_back(vector<int>());
                groups.back().push_back(i);
            }
        }

        int g1size = groups.size();

        for(int i=0; i<g1size; i++)
        {
            vector<int> &group = groups[i];
            int g2size = group.size();

            PointCluster *sig_vec = new PointCluster(*sig_vecs[group[0]]);
            vector<PointCluster> *plvec_voxel = new vector<PointCluster>(*plvec_voxels[group[0]]);
            pcl::PointCloud<PointType>::Ptr plptr = plptrs[group[0]];

            for(int j=1; j<g2size; j++)
            {
                *sig_vec += *sig_vecs[group[j]];
                const vector<PointCluster> &plvec_tem = *plvec_voxels[group[j]];

                for(int k=0; k<win_size; k++)
                    if(plvec_tem[k].N != 0)
                        (*plvec_voxel)[k] += plvec_tem[k];

                *plptr += *plptrs[group[j]];
            }

            int process_size = 0;
            for(int j=0; j<win_size; j++)
                if((*plvec_voxel)[j].N != 0)
                    process_size++;
            if(process_size < 2)
            {
                delete sig_vec; delete plvec_voxel;
                continue;
            }

            double coe = 0;
            for(int j=0; j<win_size; j++)
                coe += (*plvec_voxel)[j].N;

            voxhess.sig_vecs.push_back(sig_vec);
            voxhess.plvec_voxels.push_back(plvec_voxel);
            voxhess.coeffs.push_back(coe);
            voxhess.plptrs.push_back(plptr);
        }

    }

};

class OCTO_TREE_NODE
{
public:
    int octo_state; // 0(unknown), 1(mid node), 2(plane)
    int push_state;
    int layer;
//
    vector<PLV(3)> vec_orig, vec_tran;
    vector<PointCluster> sig_orig, sig_tran;
    PointCluster fix_point;
    PLV(3) vec_fix;

    OCTO_TREE_NODE *leaves[8];
    float voxel_center[3];
    float quater_length;

    Eigen::Vector3d center, direct, value_vector; // temporal
    double decision=-1;
    double ref;

    OCTO_TREE_NODE()
    {
        octo_state = 0; push_state = 0;   //todo force to non-cut   0(unknown), 1(mid node), 2(plane)
        vec_orig.resize(win_size); vec_tran.resize(win_size);
        sig_orig.resize(win_size); sig_tran.resize(win_size);
        for(int i=0; i<8; i++) leaves[i] = nullptr;
        ref = 255.0*rand()/(RAND_MAX + 1.0f);
        layer = 0;
    }

    bool judge_eigen(int win_count)
    {
        if(decision!=-1){
            bool cond = (decision < eigen_value_array[layer]);
            return (cond);
        }
        PointCluster covMat = fix_point;
        for(int i=0; i<win_count; i++)
            covMat += sig_tran[i];

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(covMat.cov());
        value_vector = saes.eigenvalues();
        center = covMat.v / covMat.N;
        direct = saes.eigenvectors().col(0);

        decision = saes.eigenvalues()[0] / saes.eigenvalues()[1];

        Eigen::Matrix3d cmt = covMat.P/covMat.N - center * center.transpose();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes_res(cmt);
        Eigen::Vector3d lmbd = saes_res.eigenvalues();
        double residual = covMat.N*lmbd[0];
        bool cond = (decision < eigen_value_array[layer])&&(residual<0.2);
        return (cond);
    }

    void cut_func(int ci)
    {
        PLV(3) &pvec_orig = vec_orig[ci];
        PLV(3) &pvec_tran = vec_tran[ci];

        uint a_size = pvec_tran.size();
        for(uint j=0; j<a_size; j++)
        {
            int xyz[3] = {0, 0, 0};
            for(uint k=0; k<3; k++)
                if(pvec_tran[j][k] > voxel_center[k])
                    xyz[k] = 1;
            int leafnum = 4*xyz[0] + 2*xyz[1] + xyz[2];
            if(leaves[leafnum] == nullptr)
            {
                leaves[leafnum] = new OCTO_TREE_NODE();
                leaves[leafnum]->voxel_center[0] = voxel_center[0] + (2*xyz[0]-1)*quater_length;
                leaves[leafnum]->voxel_center[1] = voxel_center[1] + (2*xyz[1]-1)*quater_length;
                leaves[leafnum]->voxel_center[2] = voxel_center[2] + (2*xyz[2]-1)*quater_length;
                leaves[leafnum]->quater_length = quater_length / 2;
                leaves[leafnum]->layer = layer + 1;
            }

            leaves[leafnum]->vec_orig[ci].push_back(pvec_orig[j]);
            leaves[leafnum]->vec_tran[ci].push_back(pvec_tran[j]);

            if(leaves[leafnum]->octo_state != 1)
            {
                leaves[leafnum]->sig_orig[ci].push(pvec_orig[j]);
                leaves[leafnum]->sig_tran[ci].push(pvec_tran[j]);
            }
        }

        PLV(3)().swap(pvec_orig); PLV(3)().swap(pvec_tran);
    }
    void recut_ransac(int win_count)
    {
        if(octo_state != 1)
        {
            int point_size = fix_point.N;
            for(int i=0; i<win_count; i++)
                point_size += sig_orig[i].N;

            push_state = 0;
            if(point_size <= min_ps){
                return;
            }

            if(judge_eigen(win_count))
            {
                push_state = 1;
            }else{
                cout<<"decision::"<<decision<<endl;
            }
        }
    }
    void recut(int win_count)
    {

        if(octo_state != 1)
        {
            int point_size = fix_point.N;
            for(int i=0; i<win_count; i++)
                point_size += sig_orig[i].N;

            push_state = 0;
            if(point_size <= min_ps){
                return;
            }

            if(judge_eigen(win_count)&&layer>0)
            {
                if(octo_state==0 && point_size>layer_size[layer])
                    octo_state = 2;

                point_size -= fix_point.N;
                if(point_size > min_ps)
                    push_state = 1;
                return;
            }
            else if(layer == layer_limit)
            {
                octo_state = 2;
                return;
            }

            octo_state = 1;
            vector<PointCluster>().swap(sig_orig);
            vector<PointCluster>().swap(sig_tran);
            for(int i=0; i<win_count; i++)
                cut_func(i);
        }
        else
            cut_func(win_count-1);

        for(int i=0; i<8; i++)
            if(leaves[i] != nullptr)
                leaves[i]->recut(win_count);
    }

    void to_margi(int mg_size, vector<IMUST> &x_poses, int win_count)
    {
        if(octo_state != 1)
        {
            if(!x_poses.empty())
                for(int i=0; i<win_count; i++)
                {
                    sig_tran[i].transform(sig_orig[i], x_poses[i]);
                    plvec_trans(vec_orig[i], vec_tran[i], x_poses[i]);
                }

            if(fix_point.N<50 && push_state==1)
                for(int i=0; i<mg_size; i++)
                {
                    fix_point += sig_tran[i];
                    vec_fix.insert(vec_fix.end(), vec_tran[i].begin(), vec_tran[i].end());
                }

            for(int i=mg_size; i<win_count; i++)
            {
                sig_orig[i-mg_size] = sig_orig[i];
                sig_tran[i-mg_size] = sig_tran[i];
                vec_orig[i-mg_size].swap(vec_orig[i]);
                vec_tran[i-mg_size].swap(vec_tran[i]);
            }

            for(int i=win_count-mg_size; i<win_count; i++)
            {
                sig_orig[i].clear(); sig_tran[i].clear();
                vec_orig[i].clear(); vec_tran[i].clear();
            }

        }
        else
            for(int i=0; i<8; i++)
                if(leaves[i] != nullptr)
                    leaves[i]->to_margi(mg_size, x_poses, win_count);

    }

    ~OCTO_TREE_NODE()
    {
        for(int i=0; i<8; i++)
            if(leaves[i] != nullptr)
                delete leaves[i];
    }

    void tras_display(pcl::PointCloud<PointType> &pl_feat, int win_count)
    {
        if(octo_state != 1)
        {

            PointType ap;
            ap.intensity = ref;

            int tsize = 0;
            for(int i=0; i<win_count; i++)
                tsize += vec_tran[i].size();
            if(tsize < 100) return;

            for(int i=0; i<win_count; i++)
                for(Eigen::Vector3d pvec : vec_tran[i])
                {
                    ap.x = pvec.x(); ap.y = pvec.y(); ap.z = pvec.z();
                    ap.normal_x = voxel_center[0];
                    ap.normal_y = voxel_center[1];
                    ap.normal_z = voxel_center[2];
                    ap.curvature = quater_length * 4;

                    pl_feat.push_back(ap);
                }

        }
        else
        {

            for(int i=0; i<8; i++)
                if(leaves[i] != nullptr)
                    leaves[i]->tras_display(pl_feat, win_count);
        }
    }

    void tras_merge(VOXEL_MERGE &vlmg, int win_count)
    {
        if(octo_state != 1)
        {
            if(push_state == 1)
            {
                pcl::PointCloud<PointType>::Ptr plptr(new pcl::PointCloud<PointType>());
                for(int i=0; i<win_count; i++)
                {
                    PointType ap; ap.intensity = i;
                    for(Eigen::Vector3d &pvec : vec_orig[i])
                    {
                        ap.x = pvec[0]; ap.y = pvec[1]; ap.z = pvec[2];
                        plptr->push_back(ap);
                    }
                }

                int psize = 0;
                for(int i=0; i<win_count; i++)
                    psize += vec_orig[i].size();

                if(psize > 100)
                    vlmg.push_voxel(&sig_orig, &fix_point, center, direct, value_vector, plptr);
            }
        }
        else
        {
            for(int i=0; i<8; i++)
                if(leaves[i] != nullptr)
                    leaves[i]->tras_merge(vlmg, win_count);
        }

    }

    void tras_opt(VOX_HESS &vox_opt, int win_count)
    {
        if(octo_state != 1)
        {
            int points_size = 0;
            for(int i=0; i<win_count; i++)
                points_size += sig_orig[i].N;

            if(points_size < min_ps){
                return;
            }

            if(push_state == 1) {
                vox_opt.push_voxel(&sig_orig, &fix_point, decision, layer);
            }

        }
        else
        {
            for(int i=0; i<8; i++)
                if (leaves[i] != nullptr)
                    leaves[i]->tras_opt(vox_opt, win_count);
        }
    }

};

class OCTO_TREE_ROOT: public OCTO_TREE_NODE {
public:
    bool is2opt;
    int life;
    vector<int> each_num;



    OCTO_TREE_ROOT() {
        is2opt = true;
        life = life_span;
        each_num.resize(win_size);
        for (int i = 0; i < win_size; i++) each_num[i] = 0;
    }
    void marginalize(int mg_size, vector<IMUST> &x_poses, int win_count)
    {
        to_margi(mg_size, x_poses, win_count);

        int left_size = 0;
        for(int i=mg_size; i<win_count; i++)
        {
            each_num[i-mg_size] = each_num[i];
            left_size += each_num[i-mg_size];
        }

        if(left_size == 0) is2opt = false;

        for(int i=win_count-mg_size; i<win_count; i++)
            each_num[i] = 0;
    }


};

inline bool iter_stop(Eigen::VectorXd &dx, double thre = 1e-7, int win_size = 0)
{
    if(win_size == 0)
        win_size = dx.rows() / 6;

    double angErr = 0, tranErr = 0;
    for(int i=0; i<win_size; i++)
    {
        angErr += dx.block<3, 1>(6*i, 0).norm();
        tranErr += dx.block<3, 1>(6*i+3, 0).norm();

    }
    angErr /= win_size; tranErr /= win_size;
    return (angErr < thre) && (tranErr < thre);
}

class BALM2
{
public:
    BALM2(){}

    double divide_thread(vector<IMUST> &x_stats, VOX_HESS &voxhess, vector<IMUST> &x_ab, Eigen::MatrixXd &Hess, Eigen::VectorXd &JacT)
    {
        int thd_num = 4;
        double residual = 0;
        Hess.setZero(); JacT.setZero();
        PLM(-1) hessians(thd_num);
        PLV(-1) jacobins(thd_num);

        for(int i=0; i<thd_num; i++)
        {
            hessians[i].resize(6*win_size, 6*win_size);
            jacobins[i].resize(6*win_size);
        }

        int tthd_num = thd_num;
        vector<double> resis(tthd_num, 0);
        int g_size = voxhess.plvec_voxels.size();
        if(g_size < tthd_num) tthd_num = 1;
        vector<thread*> mthreads(tthd_num);
        double part = 1.0 * g_size / tthd_num;

        for(int i=0; i<tthd_num; i++)
            mthreads[i] = new thread(&VOX_HESS::left_evaluate_acc2, &voxhess, x_stats, part*i, part*(i+1), ref(hessians[i]), ref(jacobins[i]), ref(resis[i]));

        for(int i=0; i<tthd_num; i++)
        {
            mthreads[i]->join();
            Hess += hessians[i];
            JacT += jacobins[i];
            residual += resis[i];
            delete mthreads[i];
        }

        return residual;
    }

    double only_residual(vector<IMUST> &x_stats, VOX_HESS &voxhess, vector<IMUST> &x_ab)
    {
        double residual1 = 0, residual2 = 0;

        voxhess.evaluate_only_residual(x_stats, residual2);
        return (residual1 + residual2);
    }
    Eigen::Matrix<double, 4, 4> LogMatrix(const Eigen::Matrix<double, 3, 3>& R_matrix,const Eigen::Matrix<double, 3, 1>& t_matrix)
    {
        Eigen::Matrix<double, 4, 4> T;
        T.block<3, 3>(0, 0)=R_matrix;
        T.block<3, 1>(0, 3)=t_matrix;
        T(3,3)=1.0;

        int N=3;
        Eigen::RealSchur<Eigen::Matrix<double, 4, 4>> schur(T);
        Eigen::Matrix<double, 4, 4> U = schur.matrixU();
        Eigen::Matrix<double, 4, 4> R = schur.matrixT();
        std::vector<bool> selected(N, true);
        Eigen::Matrix<double, 3, 3> mat_B = Eigen::Matrix<double, 3, 3>::Zero(N, N);
        Eigen::Matrix<double, 3, 3> mat_V = Eigen::Matrix<double, 3, 3>::Identity(N, N);

        for (int i = 0; i < N; i++)
        {
            if (selected[i] && fabs(R(i, i) - 1)> SAME_THRESHOLD)
            {
                int pair_second = -1;
                for (int j = i + 1; j <N; j++)
                {
                    if (fabs(R(j, j) - R(i, i)) < SAME_THRESHOLD)
                    {
                        pair_second = j;
                        selected[j] = false;
                        break;
                    }
                }
                if (pair_second > 0)
                {
                    selected[i] = false;
                    R(i, i) = R(i, i) < -1 ? -1 : R(i, i);
                    double theta = acos(R(i, i));
                    if (R(i, pair_second) < 0)
                    {
                        theta = -theta;
                    }
                    mat_B(i, pair_second) += theta;
                    mat_B(pair_second, i) += -theta;
                    mat_V(i, pair_second) += -theta / 2;
                    mat_V(pair_second, i) += theta / 2;
                    double coeff = 1 - (theta * R(i, pair_second)) / (2 * (1 - R(i, i)));
                    mat_V(i, i) += -coeff;
                    mat_V(pair_second, pair_second) += -coeff;
                }
            }
        }

        Eigen::Matrix<double, 4, 4> LogTrim = Eigen::Matrix<double, 4, 4>::Zero();
        LogTrim.block(0, 0, N, N) = mat_B;
        LogTrim.block(0, N, N, 1) = mat_V * R.block(0, N, N, 1);
        Eigen::Matrix<double, 4, 4> res = U * LogTrim * U.transpose();
        return res;
    }
    void damping_iter(vector<IMUST> &x_stats, VOX_HESS &voxhess)
    {

        double u =100, v =2;
        Eigen::MatrixXd D(6*win_size, 6*win_size), Hess(6*win_size, 6*win_size);
        Eigen::VectorXd JacT(6*win_size), dxi(6*win_size);

        D.setIdentity();
        double residual1, residual2,residual3 ,q;
        bool is_calc_hess = true;

        vector<IMUST> x_ab(win_size);
        x_ab[0] = x_stats[0];

        for(int i=1; i<win_size; i++)
        {
            x_ab[i].p = x_stats[i-1].R.transpose() * (x_stats[i].p - x_stats[i-1].p);
            x_ab[i].R = x_stats[i-1].R.transpose() * x_stats[i].R;

        }
        AndersonAcceleration accelerator_;
        Eigen::Matrix<double, 6, 1> se3,se3_AA;
        vector<IMUST> x_stats_AA = x_stats;
        bool AA_accept=0;
        Eigen::Transform<double, 3, Eigen::Affine> T;
        se3= Sophus::SE3d(x_ab[1].R, x_ab[1].p).log();
        accelerator_.init(5, 6, se3.data());


        for(int i=0; i<50; i++)
        {
            if(is_calc_hess)
                residual1 = divide_thread(x_stats, voxhess, x_ab, Hess, JacT);

            D.diagonal() = Hess.diagonal();


            dxi = (Hess + u*D).ldlt().solve(-JacT);

            for(int j=0; j<win_size; j++)
            {

                Eigen::Matrix3d dR = Exp(dxi.block<3, 1>(DVEL*j, 0));
                x_stats[j].R = dR * x_stats[j].R;
                x_stats[j].p = dR * x_stats[j].p + dxi.block<3, 1>(DVEL*j+3, 0);

            }
            residual2 = only_residual(x_stats, voxhess, x_ab);

            se3= Sophus::SE3d(x_ab[1].R, x_ab[1].p).log();
            se3_AA = accelerator_.compute(se3.data());

            x_stats_AA[1].R = Sophus::SE3d::exp(se3_AA).matrix().block<3,3>(0,0);
            x_stats_AA[1].p = Sophus::SE3d::exp(se3_AA).matrix().block<3,1>(0,3);

            residual3 = only_residual(x_stats_AA, voxhess, x_ab);

            if(residual3-residual2<0)
            {
                AA_accept=1;
                residual2=residual3;
            }
            else{
                AA_accept=0;
            }

            q = (residual1-residual2);
            double q1 = 0.5*dxi.dot(u*D*dxi-JacT);

            if(q >0)
            {
                if(AA_accept)
                {
                    x_stats = x_stats_AA;
                }
                else
                {
                    accelerator_.replace(se3.data());
                }
                q = q / q1;
                v = 2;
                q = 1 - pow(2*q-1, 3);
                u *= (q<one_three ? one_three:q);
                is_calc_hess = true;
            }
            else
            {
                accelerator_.reset(se3.data());
                u = u * v;
                v = 2 * v;
                is_calc_hess = false;

            }

//            printf("iter%d: (%lf %lf) u: %lf v: %.1lf q: %.3lf %lf %lf\n", i, residual1, residual2, u, v, q/q1, q1, q);

            if(fabs(residual1-residual2)/residual1 < 1e-8)
                break;

        }
        IMUST es0 = x_stats[0];
        for(uint i=0; i<x_stats.size(); i++)
        {
            x_stats[i].p = es0.R.transpose() * (x_stats[i].p - es0.p);
            x_stats[i].R = es0.R.transpose() * x_stats[i].R;

        }

    }
    inline Eigen::Matrix<double, 6, 1> LogToVec(const Eigen::Matrix4d& LogT)
    {
        Eigen::Matrix<double, 6, 1> res;
        res[0] = -LogT(1, 2);
        res[1] = LogT(0, 2);
        res[2] = -LogT(0, 1);
        res[3] = LogT(0, 3);
        res[4] = LogT(1, 3);
        res[5] = LogT(2, 3);
        return res;
    }
    inline Eigen::Matrix<double, 4, 4> VecToLog(const Eigen::Matrix<double, 6, 1>& v)
    {
        Eigen::Matrix<double, 4, 4>  m = Eigen::Matrix<double, 4, 4> ::Zero();
        m << 0, -v[2], v[1], v[3],
                v[2], 0, -v[0], v[4],
                -v[1], v[0], 0, v[5],
                0, 0, 0, 0;
        return m;
    }

};




#endif
