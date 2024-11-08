#ifndef MICROG_PREPRO_H
#define MICROG_PREPRO_H


#include <pcl/octree/octree_search.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/fpfh_omp.h>
#include <pcl/keypoints/iss_3d.h>
#include <pcl/filters/voxel_grid.h>
namespace prepro {
    void OctreePointCloudVoxel(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr &cloudVG,
                               pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree, double &inlTh) ;

    void OutlierRemoval(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud,
                               pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree);

    void IssKeyPointExtration(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr &ISS,
                              pcl::PointIndicesPtr &ISS_Idx, double &resolution);

    void fpfhComputation(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud,
                         double& resolution, pcl::PointIndicesPtr& iss_Idx,
                         pcl::PointCloud<pcl::FPFHSignature33>::Ptr &fpfh_out);

    void CorrespondenceSearching(pcl::PointCloud<pcl::FPFHSignature33>::Ptr &fpfhs,
                                 pcl::PointCloud<pcl::FPFHSignature33>::Ptr &fpfht, pcl::Correspondences &corr,
                                 int& max_corr, std::vector<int> &corr_NOs, std::vector<int> &corr_NOt);
}
#endif //MICROG_PREPRO_H
