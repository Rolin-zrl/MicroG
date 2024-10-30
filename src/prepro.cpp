#include "prepro.h"


void prepro::OctreePointCloudVoxel(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr &cloudVG,
                                   pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree, double &inlTh) {

    octree->setResolution(inlTh);
    octree->setInputCloud(cloud);
    octree->addPointsFromInputCloud();
    octree->getOccupiedVoxelCenters(cloudVG->points);
};

void prepro::StaticsOutlierRemoval(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud,
                                   pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>::Ptr &octree) {
    pcl::StatisticalOutlierRemoval<pcl::PointXYZ> sor;
    pcl::PointCloud<pcl::PointXYZ>::Ptr remove_cloud(new pcl::PointCloud<pcl::PointXYZ>);
    sor.setInputCloud(cloud);
    sor.setMeanK(20);
    sor.setStddevMulThresh(1.5);
    sor.filter(*cloud);
    sor.setNegative(true);
    sor.filter(*remove_cloud);
    for (size_t i = 0; i < remove_cloud->size(); i++)
        octree->deleteVoxelAtPoint(remove_cloud->at(i));
}

void prepro::IssKeyPointExtration(pcl::PointCloud<pcl::PointXYZ>::Ptr &cloud, pcl::PointCloud<pcl::PointXYZ>::Ptr &ISS,
                                  pcl::PointIndicesPtr &ISS_Idx, double &resolution) {

    double iss_salient_radius_ = 6 * resolution; //8
    double iss_non_max_radius_ = 4 * resolution;  //3
    double iss_gamma_21_(0.975);
    double iss_gamma_32_(0.975);
    double iss_min_neighbors_(8);
    int iss_threads_(16); //switch to the number of threads in your cpu for acceleration

    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
    pcl::ISSKeypoint3D<pcl::PointXYZ, pcl::PointXYZ> iss_detector;
    iss_detector.setSearchMethod(tree);
    iss_detector.setSalientRadius(iss_salient_radius_);
    iss_detector.setNonMaxRadius(iss_non_max_radius_);
    iss_detector.setThreshold21(iss_gamma_21_);
    iss_detector.setThreshold32(iss_gamma_32_);
    iss_detector.setMinNeighbors(iss_min_neighbors_);
    iss_detector.setNumberOfThreads(iss_threads_);
    iss_detector.setInputCloud(cloud);
    iss_detector.compute(*ISS);
    ISS_Idx->indices = iss_detector.getKeypointsIndices()->indices;
    ISS_Idx->header = iss_detector.getKeypointsIndices()->header;
}

void prepro::fpfhComputation(pcl::PointCloud<pcl::PointXYZ>::Ptr& cloud, double& resolution, pcl::PointIndicesPtr& iss_Idx, pcl::PointCloud<pcl::FPFHSignature33>::Ptr &fpfh_out) {
    pcl::PointCloud<pcl::Normal>::Ptr normal_(new pcl::PointCloud<pcl::Normal>);
    //compute normal
    pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
//    ne.setNumberOfThreads(16);
    ne.setInputCloud(cloud);
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(3 * resolution);//4
    ne.compute(*normal_);


    //compute fpfh using normals
    pcl::FPFHEstimationOMP<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_est;
    fpfh_est.setInputCloud(cloud);
    fpfh_est.setInputNormals(normal_);
    fpfh_est.setSearchMethod(tree);
    fpfh_est.setRadiusSearch(8 * resolution);  //14
    fpfh_est.setNumberOfThreads(16);
    fpfh_est.setIndices(iss_Idx);
    fpfh_est.compute(*fpfh_out);
}


void prepro::CorrespondenceSearching(pcl::PointCloud<pcl::FPFHSignature33>::Ptr &fpfhs,
                                     pcl::PointCloud<pcl::FPFHSignature33>::Ptr &fpfht, pcl::Correspondences &corr,
                                     int &max_corr, std::vector<int> &corr_NOs, std::vector<int> &corr_NOt) {
    int n = std::min(max_corr,
                     (int) fpfht->size()); //maximum number of correspondences to find for each source point
    corr.clear();
    corr_NOs.assign(fpfhs->size(), 0);
    corr_NOt.assign(fpfht->size(), 0);
    // Use a KdTree to search for the nearest matches in feature space
    pcl::KdTreeFLANN<pcl::FPFHSignature33> treeS;
    treeS.setInputCloud(fpfhs);
    pcl::KdTreeFLANN<pcl::FPFHSignature33> treeT;
    treeT.setInputCloud(fpfht);
    for (size_t i = 0; i < fpfhs->size(); i++) {
        std::vector<int> corrIdxTmp(n);
        std::vector<float> corrDisTmp(n);
        //find the best n matches in target fpfh
        treeT.nearestKSearch(*fpfhs, i, n, corrIdxTmp, corrDisTmp);
        for (size_t j = 0; j < corrIdxTmp.size(); j++) {
            bool removeFlag = true;
            int searchIdx = corrIdxTmp[j];
            std::vector<int> corrIdxTmpT(n);
            std::vector<float> corrDisTmpT(n);
            treeS.nearestKSearch(*fpfht, searchIdx, n, corrIdxTmpT, corrDisTmpT);
            for (size_t k = 0; k < n; k++) {
                if (corrIdxTmpT.data()[k] == i) {
                    removeFlag = false;
                    break;
                }
            }
            if (removeFlag == false) {
                pcl::Correspondence corrTabTmp;
                corrTabTmp.index_query = i;
                corrTabTmp.index_match = corrIdxTmp[j];
                corrTabTmp.distance = corrDisTmp[j];
                corr.push_back(corrTabTmp);
                corr_NOs[i]++;
                corr_NOt[corrIdxTmp[j]]++;
            }
        }
    }
}