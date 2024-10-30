#include <vector>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/octree/octree_search.h>
#include "tools.hpp"
#include "bavoxel.hpp"
#include <malloc.h>


#ifndef MICROG_FINEREG_H
#define MICROG_FINEREG_H

namespace finereg{

    class FineReg{
    public:
        VOX_HESS voxhess;
        void setOctree(pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree_S_,pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree_T_);
        void setCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT_);
        void setCorr(pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_);
        void setTrans(Eigen::Matrix4d& InitTrans_);
        void setResolution(double& res_);
        void setL(double & L_);
        void AdaptiveSearch();
        void Opt_AA(Eigen::Matrix4d& FineTrans_);

    private:
        pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr octree_S;
        pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr octree_T;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudS;
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT;
        pcl::PointCloud<pcl::PointXYZ>::Ptr corS;
        pcl::PointCloud<pcl::PointXYZ>::Ptr corT;
        double resolution;
        double L;
        Eigen::Matrix4d InitTrans;
        std::vector<IMUST> x_buf;
        unordered_map < VOXEL_LOC, OCTO_TREE_ROOT * > feat_map;
        };
}





#endif //MICROG_FINEREG_H
