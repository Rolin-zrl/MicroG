
// ================================================================================================
//  Copyright (C) 2024, SCHOOL OF GEODESY AND GEOMATICS, WUHAN UNIVERSITY, CHINA
//  All Rights Reserved.
//
//  Authors: Rongling Zhang, Pengcheng Wei et al.
//  Contact:hitzrl@foxmail.com
//  See LICENSE file for more detailed license information.
// ================================================================================================

#ifndef MICROG_MICROG_H
#define MICROG_MICROG_H

#include <iostream>
#include <fstream>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/correspondence.h>
#include <pcl/keypoints/iss_3d.h>
#include <pcl/visualization/pcl_visualizer.h>


#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "prepro.h"
#include "gror.hpp"
#include "opt.hpp"
#include "finereg.h"
#include "evaluation.hpp"


namespace c2r {

    class MicroG {
    public:

        static int  DEBUG_DETAILS;

        std::string source_path;
        std::string target_path;
        std::string gt_path;
        std::string result_path;
        double resolution;
        int K;
        double L;

        void ReadData();
        void PreProcessing();
        void CoarseREG();
        void FineREG();

        void Evaluation();
        void RecordFiles();
        void Visualization();

    private:
        Eigen::Matrix4d CoarseTrans=Eigen::Matrix4d::Identity();
        Eigen::Matrix4d FineTrans=Eigen::Matrix4d::Identity();
        Eigen::Matrix4d GtTrans=Eigen::Matrix4d::Identity();
        pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr octree_S;
        pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr octree_T;
        pcl::PointCloud<pcl::PointXYZ>::Ptr origin_cloudS;
        pcl::PointCloud<pcl::PointXYZ>::Ptr origin_cloudT;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS;   //downsapmling
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT;
        pcl::PointCloud<pcl::PointXYZ>::Ptr corS;   //initial correspondences
        pcl::PointCloud<pcl::PointXYZ>::Ptr corT;
        pcl::CorrespondencesPtr corr;
        pcl::PointCloud<pcl::PointXYZ>::Ptr inliers_corS;   //Removal_based_hierarchical_strategy
        pcl::PointCloud<pcl::PointXYZ>::Ptr inliers_corT;
        pcl::CorrespondencesPtr corr_based_nodes;
        pcl::CorrespondencesPtr corr_based_edges;
        pcl::PointCloud<pcl::PointXYZ>::Ptr issS;  //ISS
        pcl::PointCloud<pcl::PointXYZ>::Ptr issT;

        double Crmse,Frmse,CRerror,FRerror,CTerror,FTerror;
        std::chrono::time_point<std::chrono::system_clock> t0,t1;
        double time_gror,time_opt,time_search,time_optAA;
    };

    int OutliersRate(pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_,Eigen::Matrix4d &GtTrans_,double &resolution_);
    void Corr2Point( pcl::CorrespondencesPtr& corr_,pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_);
    void ShowVGPointCloudCorr ( const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS,  const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT, const pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_,const pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_,const std::string& viewer_name);
    void ShowVGPointCloud ( const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS,  const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT,const std::string& viewer_name);
} // C2R

#endif //MICROG_MICROG_H
