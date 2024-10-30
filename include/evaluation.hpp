#pragma once
#ifndef MICROG_EVALUATION_HPP
#define MICROG_EVALUATION_HPP


#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/common/distances.h>
#include <pcl/kdtree/kdtree_flann.h>


namespace evaluation {

    inline void CalculateRMSE(pcl::PointCloud<pcl::PointXYZ>::Ptr source_cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr target_cloud,double& rmse) {
        pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
        kdtree.setInputCloud(target_cloud);

        double sum_squared_error = 0.0;
        int valid_correspondences = 0;

        for (const auto& point : source_cloud->points) {
            std::vector<int> point_idx(1);
            std::vector<float> point_squared_distance(1);

            if (kdtree.nearestKSearch(point, 1, point_idx, point_squared_distance) > 0) {
                sum_squared_error += point_squared_distance[0];
                valid_correspondences++;
            }
        }

        rmse = std::sqrt(sum_squared_error / valid_correspondences)*100;
    };

    inline void  CalculateRTerror(const Eigen::Matrix4d& j, const Eigen::Matrix4d& i, double& e_t, double& e_r_f){

        Eigen::Matrix3d R_j = j.block(0, 0, 3, 3);
        Eigen::Vector3d t_j = j.block(0, 3, 3, 1);

        Eigen::Matrix3d R_i = i.block(0, 0, 3, 3);
        Eigen::Vector3d t_i = i.block(0, 3, 3, 1);
        Eigen::Matrix3d R = R_j * R_i.transpose();
        Eigen::Vector3d t = t_j - t_i;

        Eigen::Affine3d rot(R);
        e_r_f = Eigen::AngleAxisd(rot.rotation()).angle();
        e_r_f=e_r_f*180/M_PI;
        e_t = t.norm()*100;
    };



};
#endif //MICROG_EVALUATION_HPP