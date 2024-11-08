// ================================================================================================
//  Copyright (C) 2024, SCHOOL OF GEODESY AND GEOMATICS, WUHAN UNIVERSITY, CHINA
//  All Rights Reserved.
//
//  Authors: Rongling Zhang, Pengcheng Wei et al.
//  Contact:hitzrl@foxmail.com
//  See LICENSE file for more detailed license information.
// ================================================================================================

#include "../include/MicroG.h"


namespace c2r {

    int  MicroG::DEBUG_DETAILS;


    void MicroG::ReadData(){
        origin_cloudS.reset(new pcl::PointCloud<pcl::PointXYZ>());
        origin_cloudT.reset(new pcl::PointCloud<pcl::PointXYZ>());
        if (source_path.substr(source_path.find_last_of('.') + 1) == "pcd") {
            pcl::io::loadPCDFile(source_path, *origin_cloudS);
            pcl::io::loadPCDFile(target_path, *origin_cloudT);
        }else if (source_path.substr(source_path.find_last_of('.') + 1) == "ply") {
            pcl::io::loadPLYFile(source_path, *origin_cloudS);
            pcl::io::loadPLYFile(target_path, *origin_cloudT);
        }

        if(gt_path.empty()){
            std::cout << "**********SUCCESS: POINTCLOUD LOADED WITHOUT GROUND TRUTH**********" << std::endl;
            return;
        } else{
            std::cout << "**********SUCCESS: POINTCLOUD LOADED WITH GROUND TRUTH**********" << std::endl;
        }

        std::ifstream infile;
        infile.open(gt_path, std::ios::in);
        if (!infile.is_open()) {
            std::cout << "**********ERROR: FAILED TO READ GROUND TRUTH**********" << std::endl;
            exit(0);
        }
        size_t numRow = GtTrans.rows();
        size_t numCol = GtTrans.cols();
        for (int j = 0; j < numRow; j++)//numRow
        {
            for (int i = 0; i < numCol; i++)//numCol
            {
                infile >> GtTrans(j, i);
            }

        }
    }

    void MicroG::PreProcessing() {

        if(DEBUG_DETAILS) { t0 = std::chrono::system_clock::now();}

        int max_corr = 5;// neighbor number in descriptor searching
        std::cout << "**********STEP 1.0：PREPROCESSING**********" << std::endl;

        if(DEBUG_DETAILS) {
            t1 = std::chrono::system_clock::now();
            std::cout<<"-----NUM_ORIGIN_CLOUDS: "<<origin_cloudS->size()<<"-----"<<std::endl;
            std::cout<<"-----NUM_ORIGIN_CLOUDT: "<<origin_cloudT->size()<<"-----"<<std::endl;
        }
        //downsampling
        octree_S.reset(new pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>(resolution));
        octree_T.reset(new pcl::octree::OctreePointCloudSearch<pcl::PointXYZ>(resolution));
        cloudS.reset(new pcl::PointCloud<pcl::PointXYZ>);
        cloudT.reset(new pcl::PointCloud<pcl::PointXYZ>);
        prepro::OctreePointCloudVoxel(origin_cloudS, cloudS, octree_S, resolution);
        prepro::OctreePointCloudVoxel(origin_cloudT, cloudT, octree_T, resolution);
        //filter
        prepro::OutlierRemoval(cloudS, octree_S);
        prepro::OutlierRemoval(cloudT, octree_T);

        if(DEBUG_DETAILS) {
            t1 = std::chrono::system_clock::now();
            std::cout<<"-----TIME_DOWNSAMPLING: "<<double(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) / 1000.0<<"s-----"<<std::endl;
            std::cout<<"-----NUM_CLOUDS: "<<cloudS->size()<<"-----"<<std::endl;
            std::cout<<"-----NUM_CLOUDT: "<<cloudT->size()<<"-----"<<std::endl;
        }
        //extract ISS
        issS.reset(new pcl::PointCloud<pcl::PointXYZ>());
        issT.reset(new pcl::PointCloud<pcl::PointXYZ>());
        pcl::PointIndicesPtr issS_Idx(new pcl::PointIndices);
        pcl::PointIndicesPtr issT_Idx(new pcl::PointIndices);
        prepro::IssKeyPointExtration(cloudS, issS, issS_Idx, resolution);
        prepro::IssKeyPointExtration(cloudT, issT, issT_Idx, resolution);
        if(DEBUG_DETAILS) {
            t0 = std::chrono::system_clock::now();
            std::cout<<"-----TIME_EXTRACT_ISS: "<<double(std::chrono::duration_cast<std::chrono::milliseconds>(t0 - t1).count()) / 1000.0<<"s-----"<<std::endl;
            std::cout<<"-----NUM_ISS_CLOUDS: "<<issS->size()<<"-----"<<std::endl;
            std::cout<<"-----NUM_ISS_CLOUDT: "<<issT->size()<<"-----"<<std::endl;
        }

        //calculate FPFH
        pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhS(new pcl::PointCloud<pcl::FPFHSignature33>());
        pcl::PointCloud<pcl::FPFHSignature33>::Ptr fpfhT(new pcl::PointCloud<pcl::FPFHSignature33>());
        prepro::fpfhComputation(cloudS, resolution, issS_Idx, fpfhS);
        prepro::fpfhComputation(cloudT, resolution, issT_Idx, fpfhT);
        if(DEBUG_DETAILS) {
            t1 = std::chrono::system_clock::now();
            std::cout<<"-----TIME_CACULATE_FPFH: "<<double(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) / 1000.0<<"s-----"<<std::endl;
        }

        //match correspondences
        corr.reset(new pcl::Correspondences());
        std::vector<int> corr_NOS, corr_NOT;
        prepro::CorrespondenceSearching(fpfhS, fpfhT, *corr, max_corr, corr_NOS, corr_NOT);

        corS.reset(new pcl::PointCloud<pcl::PointXYZ>());
        corT.reset(new pcl::PointCloud<pcl::PointXYZ>());
        Corr2Point(corr,issS,issT,corS,corT);

        if(DEBUG_DETAILS) {
            t0 = std::chrono::system_clock::now();
            std::cout<<"-----TIME_MATCH_FPFH: "<<double(std::chrono::duration_cast<std::chrono::milliseconds>(t0 - t1).count()) / 1000.0<<"s-----"<<std::endl;
            std::cout<<"-----NUM_CORR: "<<corr->size()<<"-----"<<std::endl;
            if(!result_path.empty()){
                std::cout<<"-----RATE_INLIERS: "<<OutliersRate(corS,corT,GtTrans,resolution)<<"/"<<corS->size()<<"-----"<<std::endl;
            }
        }
    }

    void MicroG::CoarseREG(){
        std::cout << "**********STEP 2.0：COARSE REG**********"  << std::endl;
        std::cout << "**********STEP 2.1：REMOVAL OUTLIERS**********"  << std::endl;

        if(DEBUG_DETAILS) { t0 = std::chrono::system_clock::now();}
        //removal
        corr_based_nodes.reset(new pcl::Correspondences());
        pcl::PointCloud<pcl::PointXYZ>::Ptr pcs(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::registration::GRORInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, float> gror;
        gror.setInputSource(issS);
        gror.setInputTarget(issT);
        gror.setResolution(resolution);
        gror.setOptimalSelectionNumber(K);
        gror.setInputCorrespondences(corr);
        gror.align(*pcs);
        gror.get_edge_InlierCorrespondences(corr_based_edges);
        inliers_corS.reset(new pcl::PointCloud<pcl::PointXYZ>());
        inliers_corT.reset(new pcl::PointCloud<pcl::PointXYZ>());
        Corr2Point(corr_based_edges,issS,issT,inliers_corS,inliers_corT);
        if(DEBUG_DETAILS) {
            t1 = std::chrono::system_clock::now();
            time_gror=double(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) / 1000.0;
            std::cout<<"-----TIME_REMOVAL: "<<time_gror<<"s-----"<<std::endl;
            std::cout<<"-----NUM_CORR_INLIERS: "<<corr_based_edges->size()<<"-----"<<std::endl;
            if(!result_path.empty()){
                std::cout<<"-----RATE_INLIERS: "<<OutliersRate(inliers_corS,inliers_corT,GtTrans,resolution)<<"/"<<inliers_corT->size()<<"-----"<<std::endl;
            }
        }
        if(corr_based_edges->size()<5){
            std::cout<<"*****WARNING: CORRESPONDENCE SIZE TOO SMALL*****"<<std::endl;
        }
        //GNC-Welsch Optimization
        std::cout << "**********STEP 2.2：GNC-WELSCH OPTIMIZATION**********"  << std::endl;
        opt::Opt OPT;
        CoarseTrans=OPT.GNCWelsch_Opt(inliers_corS,inliers_corT);
        if(DEBUG_DETAILS) {
            t0 = std::chrono::system_clock::now();
            time_opt=double(std::chrono::duration_cast<std::chrono::milliseconds>(t0 - t1).count()) / 1000.0;
            std::cout<<"-----TIME_GNC-Welsch_Optimization: "<<time_opt<<"s-----"<<std::endl;
            std::cout<<"-----TIME_COARSE_REG: "<<time_gror+time_opt<<"s-----"<<std::endl;
        }
    };

    void MicroG::FineREG(){
        std::cout << "**********STEP 3.0：FINE REG**********"  << std::endl;
        if(DEBUG_DETAILS) {
            t0 = std::chrono::system_clock::now();
        }
        // setting input
        finereg::FineReg fr;
        fr.setOctree(octree_S,octree_T);
        fr.setCloud(origin_cloudS,origin_cloudT);
        fr.setCorr(inliers_corS,inliers_corT);
        fr.setTrans(CoarseTrans);
        fr.setResolution(resolution);
        fr.setL(L);

        std::cout << "**********STEP 3.1：ADAPTIVE SEARCH PLANAR**********"  << std::endl;

        fr.AdaptiveSearch();

        if(DEBUG_DETAILS) {
            t1 = std::chrono::system_clock::now();
            time_search=double(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) / 1000.0;
            std::cout<<"-----NUM_PLANARS: "<<fr.voxhess.plvec_voxels.size()<<"-----"<<std::endl;
            std::cout<<"-----TIME_ADAPTIVE_SEARCH: "<<time_search<<"s-----"<<std::endl;
        }

        std::cout << "**********STEP 3.2：JOINT OPTIMIZATION WITH AA**********"  << std::endl;

        fr.Opt_AA(FineTrans);

        if(DEBUG_DETAILS) {
            t0 = std::chrono::system_clock::now();
            time_optAA=double(std::chrono::duration_cast<std::chrono::milliseconds>(t0 - t1).count()) / 1000.0;
            std::cout<<"-----TIME_OPT_AA: "<<time_optAA<<"s-----"<<std::endl;
            std::cout<<"-----TIME_FINE_REG: "<<time_search+time_optAA<<"s-----"<<std::endl;
            std::cout<<"-----TIME_TOTAL_REG: "<<time_gror+time_opt+time_search+time_optAA<<"s-----"<<std::endl;
        }
    };

    void MicroG::Evaluation() {

        std::cout << "**********EVALUATION**********"  << std::endl;
        pcl::PointCloud<pcl::PointXYZ>::Ptr creg_cloudT(new  pcl::PointCloud<pcl::PointXYZ>);
        pcl::PointCloud<pcl::PointXYZ>::Ptr freg_cloudT(new  pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*cloudT, *creg_cloudT, CoarseTrans);
        pcl::transformPointCloud(*cloudT, *freg_cloudT, FineTrans);
        evaluation::CalculateRMSE(cloudS,creg_cloudT,Crmse);
        evaluation::CalculateRMSE(cloudS,freg_cloudT,Frmse);
        if(DEBUG_DETAILS) {
            std::cout<<"-----RMSE_COARSE(cm): "<<Crmse<<"-----"<<std::endl;
            std::cout<<"-----RMSE_FINE(cm): "<<Crmse<<"-----"<<std::endl;
        }
        if(!gt_path.empty()){
            evaluation::CalculateRTerror(CoarseTrans,GtTrans.inverse(),CTerror,CRerror);
            evaluation::CalculateRTerror(FineTrans,GtTrans.inverse(),FTerror,FRerror);
            if(DEBUG_DETAILS) {
                std::cout<<"-----TRANSLATION_ERROR_COARSE(cm): "<<CTerror<<"-----"<<std::endl;
                std::cout<<"-----ROTATION_ERROR_COARSE(deg): "<<CRerror<<"-----"<<std::endl;
                std::cout<<"-----TRANSLATION_ERROR_COARSE(cm): "<<FTerror<<"-----"<<std::endl;
                std::cout<<"-----ROTATION_ERROR_COARSE(deg): "<<FRerror<<"-----"<<std::endl;
            }
        }


    };

    void MicroG::RecordFiles() {

        size_t lastSlash = source_path.find_last_of('/');
        size_t lastDot = source_path.find_last_of('.');
        size_t lastSlash2 = target_path.find_last_of('/');
        size_t lastDot2 = target_path.find_last_of('.');
        std::string pair_name=source_path.substr(lastSlash + 1, lastDot - lastSlash - 1)+"-"+target_path.substr(lastSlash2 + 1, lastDot2 - lastSlash2 - 1);

        std::string result_csv = result_path + "/" +pair_name+ ".csv";
        std::ofstream outFile;
        outFile.open(result_csv.c_str(), ios::out);
        outFile.setf(ios::fixed, ios::floatfield);
        outFile <<"scene"<<","<< "pair_name" << ',' << "Coarse_R_error(deg)" << ',' << "Coarse_T_error(cm)"<< ',' <<"CR_time(s)"<<","<< "Fine_R_error(deg)" << ',' << "Fine_T_error(cm)" <<","<<"FR_time(s)"<<"," "Success" <<","<<"Total_time(s)"<<std::endl;
        outFile <<"pair" << "," <<"\""<<pair_name <<"\""<< ',' <<CRerror << "," << CTerror<<","<<time_gror+time_opt<<","<<FRerror<<","<<FTerror<<","<<time_search+time_optAA<<","<<1<<","<<time_gror+time_opt+time_search+time_optAA<<std::endl;
        outFile.close();

        std::string result_txt = result_path + "/" +pair_name+ ".txt";
        std::ofstream outTXT;
        outTXT.open(result_txt.c_str(), ios::out);
        outTXT.setf(ios::fixed, ios::floatfield);
        outTXT<<CoarseTrans<<std::endl<<std::endl<<FineTrans<<std::endl;

    }

    void MicroG::Visualization() {

         pcl::visualization::PCLVisualizer::Ptr viewerVGF(new pcl::visualization::PCLVisualizer);
        int v1(0), v2(0), v3(0), v4(0);

        viewerVGF->createViewPort(0.0,0.5,0.5,1.0,v1);
        viewerVGF->createViewPort(0.5,0.5,1.0,1.0,v2);
        viewerVGF->createViewPort(0.5,0.0,1.0,0.5,v3);
        viewerVGF->createViewPort(0.0,0.0,0.5,0.5,v4);
        pcl::PointCloud<pcl::PointXYZ>::Ptr creg_cloudT(new  pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*cloudT, *creg_cloudT, CoarseTrans);
        pcl::PointCloud<pcl::PointXYZ>::Ptr freg_cloudT(new  pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*cloudT, *freg_cloudT, FineTrans);

        ShowVGPointCloudCorr(viewerVGF,cloudS, cloudT,corS,corT,"Initial Correspondences",v1);
        ShowVGPointCloudCorr(viewerVGF,cloudS, cloudT,inliers_corS,inliers_corT,"After Removal Outliers",v2);
        ShowVGPointCloud(viewerVGF,cloudS, creg_cloudT,"Coarse Reg",v3);
        ShowVGPointCloud(viewerVGF,cloudS, freg_cloudT,"Fine Reg",v4);
        while (!viewerVGF->wasStopped()) {
            viewerVGF->spinOnce(100);
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    };

    int OutliersRate(pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_,Eigen::Matrix4d &GtTrans_,double &resolution_){
        int Num_inliers=0;
        pcl::PointCloud<pcl::PointXYZ>::Ptr corT_Reg(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*corT_,*corT_Reg,GtTrans_.cast<float>().inverse());
        for(int i=0;i<corS_->size();i++){
            if((corS_->points[i].getVector3fMap()-corT_Reg->points[i].getVector3fMap()).norm()<2*resolution_)
                Num_inliers++;
        }
        return Num_inliers;
    };
    void Corr2Point( pcl::CorrespondencesPtr& corr_,pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_,pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_) {
        for (int i = 0; i < corr_->size(); ++i) {
            int idx_S = (*corr_)[i].index_query;
            int idx_T = (*corr_)[i].index_match;
            corS_->push_back((*cloudS_)[idx_S]);
            corT_->push_back((*cloudT_)[idx_T]);
        }
    };
    void ShowVGPointCloudCorr (pcl::visualization::PCLVisualizer::Ptr& viewerVGF, const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS,  const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT, const pcl::PointCloud<pcl::PointXYZ>::Ptr& corS_, const pcl::PointCloud<pcl::PointXYZ>::Ptr& corT_,const std::string& viewer_name,int v)
    {
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorS(cloudS, 222.0, 185.0, 0.0); //blue
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorT(cloudT, 4.0, 151.0, 210.0); //green
        viewerVGF->setBackgroundColor(255, 255, 255);
        viewerVGF->addPointCloud<pcl::PointXYZ>(cloudS,colorS,viewer_name+"S",v);
        viewerVGF->addPointCloud<pcl::PointXYZ>(cloudT,colorT,viewer_name+"T",v);
        viewerVGF->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, viewer_name+"S");
        viewerVGF->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, viewer_name+"T");
        for(int i=0;i<corS_->size();i++)
        {
            std::stringstream ss;
            ss<<viewer_name<<i;
            viewerVGF->addLine(corS_->points[i],corT_->points[i],230/255.0,52/255.0,4/255.0,ss.str(),v);
        }
        viewerVGF->addText(viewer_name,10,10,20,0,0, 0,viewer_name,v);
    };
    void ShowVGPointCloud ( pcl::visualization::PCLVisualizer::Ptr& viewerVGF,const pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudS, const  pcl::PointCloud<pcl::PointXYZ>::Ptr& cloudT,const std::string& viewer_name,int v) {

        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorS(cloudS, 222.0, 185.0, 0.0); //blue
        pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> colorT(cloudT, 4.0, 151.0, 210.0); //green
        viewerVGF->setBackgroundColor(255, 255, 255);
        viewerVGF->addPointCloud<pcl::PointXYZ>(cloudS, colorS,viewer_name+"S",v);
        viewerVGF->addPointCloud<pcl::PointXYZ>(cloudT, colorT,viewer_name+"T",v);
        viewerVGF->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, viewer_name+"S");
        viewerVGF->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 1, viewer_name+"T");
        viewerVGF->addText(viewer_name,10,10,20,0,0, 0,viewer_name,v);
    }




} // C2R