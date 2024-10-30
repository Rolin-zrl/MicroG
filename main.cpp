// ================================================================================================
//  Copyright (C) 2024, SCHOOL OF GEODESY AND GEOMATICS, WUHAN UNIVERSITY, CHINA
//  All Rights Reserved.
//
//  Authors: Rongling Zhang, Pengcheng Wei et al.
//  Contact:hitzrl@foxmail.com
//  See LICENSE file for more detailed license information.
// ================================================================================================

#include "include/MicroG.h"

int main(int argc, char** argv) {

    c2r::MicroG microG;

    //INPUT
    microG.result_path="../data/";
    microG.source_path=argv[1];
    microG.target_path=argv[2];
    if(argc==6) {
        microG.resolution = atof(argv[3]);  //depend on pointcloud 0.1
        microG.K=atoi(argv[4]);      //default 800
        microG.L=atof(argv[5]);     //default 2
    }else if(argc==7) {
        microG.gt_path = argv[3];
        microG.resolution = atof(argv[4]);
        microG.K = atoi(argv[5]);
        microG.L = atof(argv[6]);
    }else{
        std::cout<<"**********ERROR: INSUFFICIENT INPUT PARAMETERS**********"<<std::endl;
    }



    c2r::MicroG::DEBUG_DETAILS=1;

    microG.ReadData();
    microG.PreProcessing();
    microG.CoarseREG();
    microG.FineREG();

//    microG.Evaluation(); //T*target
//    microG.RecordFiles();
    microG.Visualization();



    return 0;
}