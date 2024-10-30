#include "finereg.h"

namespace finereg {
    void FineReg::AdaptiveSearch() {

        x_buf.resize(2);
        x_buf[0].R=Eigen::Matrix3d::Identity();
        x_buf[0].p.z();
        x_buf[0].t=1;
        x_buf[1].R=InitTrans.block<3, 3>(0, 0);
        x_buf[1].p=InitTrans.block<3, 1>(0, 3);
        x_buf[1].t=1;
        win_size = x_buf.size();

        pcl::PointCloud<pcl::PointXYZ>::Ptr voxel_cloudS(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointCloud<pcl::PointXYZ>::Ptr voxel_cloudT(new pcl::PointCloud<pcl::PointXYZ>);

        double VoxelSize = L * resolution;

        for (int i = 0; i < corS->size(); i++) {
            // build map
            Eigen::Vector3f center = 0.5 * (corS->points[i].getVector3fMap() +(InitTrans.block<3, 3>(0, 0).cast<float>() *corT->points[i].getVector3fMap() +InitTrans.block<3, 1>(0, 3).cast<float>()));

            float loc_xyz[3];
            for (int j = 0; j < 3; j++) {
                loc_xyz[j] = center[j] / VoxelSize;
                if (loc_xyz[j] < 0) loc_xyz[j] -= 1.0;
            }
            VOXEL_LOC position((int64_t) loc_xyz[0], (int64_t) loc_xyz[1], (int64_t) loc_xyz[2]);
            auto iter = feat_map.find(position);
            if (iter != feat_map.end()) {
                continue;
            }

            std::vector<int> pointIdxVec_S;
            std::vector<int> pointIdxVec_T;
            std::vector<float> dis_S;
            std::vector<float> dis_T;

            octree_S->radiusSearch(corS->points[i], L * resolution, pointIdxVec_S,dis_S);
            octree_T->radiusSearch(corT->points[i], L * resolution, pointIdxVec_T, dis_T);

            if (pointIdxVec_S.size() == 0 || pointIdxVec_T.size() == 0) {
                continue;
            }


            for (int j: pointIdxVec_S) {
                voxel_cloudS->push_back((*cloudS)[j]);
            }
            for (int k: pointIdxVec_T) {
                voxel_cloudT->push_back((*cloudT)[k]);
            }

            if (voxel_cloudS->size() + voxel_cloudT->size() < 8) {
                continue;
            }

            OCTO_TREE_ROOT *ot = new OCTO_TREE_ROOT();
            ot->voxel_center[0] = loc_xyz[0] * VoxelSize;
            ot->voxel_center[1] = loc_xyz[1] * VoxelSize;
            ot->voxel_center[2] = loc_xyz[2] * VoxelSize;
            ot->quater_length = VoxelSize / 2.0;

            for ( auto& point: *voxel_cloudS) {
                Eigen::Vector3d pvec_orig = Eigen::Vector3d(point.x, point.y, point.z);
                Eigen::Vector3d pvec_tran = x_buf[0].R * pvec_orig + x_buf[0].p;
                ot->vec_orig[0].push_back(pvec_orig);
                ot->sig_orig[0].push(pvec_orig);
                ot->vec_tran[0].push_back(pvec_tran);
                ot->sig_tran[0].push(pvec_tran);
                ot->each_num[0]++;
            }
            for ( auto& point: *voxel_cloudT) {
                Eigen::Vector3d pvec_orig = Eigen::Vector3d(point.x, point.y, point.z) ;
                Eigen::Vector3d pvec_tran = x_buf[1].R * pvec_orig + x_buf[1].p;

                ot->vec_orig[1].push_back(pvec_orig);
                ot->sig_orig[1].push(pvec_orig);
                ot->vec_tran[1].push_back(pvec_tran);
                ot->sig_tran[1].push(pvec_tran);
                ot->each_num[1]++;
            }

            feat_map[position] = ot;

            voxel_cloudS.reset(new pcl::PointCloud<pcl::PointXYZ>);
            voxel_cloudT.reset(new pcl::PointCloud<pcl::PointXYZ>);
        }

        for(auto iter=feat_map.begin(); iter!=feat_map.end() ; iter++)
        {

            iter->second->recut(win_size);
            iter->second->tras_opt(voxhess, win_size);
        }

        std::vector<int> planes(x_buf.size(), 0);
        for(int i=0; i<voxhess.plvec_voxels.size(); i++)
        {
            for(int j=0; j<voxhess.plvec_voxels[i]->size(); j++)
                if(voxhess.plvec_voxels[i]->at(j).N != 0)
                    planes[j]++;
        }
        sort(planes.begin(), planes.end());

        if(planes[0] < 8)
        {
            std::cout<<"*****WARNING: PLANAR SIZE TOO SMALL*****"<<std::endl;
        }

    };

    void FineReg::Opt_AA(Eigen::Matrix4d& FineTrans_){

        BALM2 opt_lsv;
        opt_lsv.damping_iter(x_buf, voxhess);

        FineTrans_.block<3, 3>(0, 0)=x_buf[1].R;
        FineTrans_.block<3, 1>(0, 3)=x_buf[1].p;
        FineTrans_(3,3)=1;
        for(auto iter=feat_map.begin(); iter!=feat_map.end();)
        {
            delete iter->second;
            feat_map.erase(iter++);
        }
        feat_map.clear();
        malloc_trim(0);
    };

    void FineReg::setOctree(pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree_S_,
                            pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree_T_) {octree_S=octree_S_;octree_T=octree_T_;};

    void FineReg::setCorr(pcl::PointCloud<pcl::PointXYZ>::Ptr &corS_,pcl::PointCloud<pcl::PointXYZ>::Ptr &corT_) {corS=corS_;corT=corT_;};

    void FineReg::setResolution(double &res_) {resolution=res_;};

    void FineReg::setL(double &L_) {L=L_;}

    void FineReg::setTrans(Eigen::Matrix4d &InitTrans_) {InitTrans=InitTrans_;}

    void FineReg::setCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloudS_, pcl::PointCloud<pcl::PointXYZ>::Ptr &cloudT_) {cloudS=cloudS_;cloudT=cloudT_;};

}