#ifndef MICROG_OPT_HPP
#define MICROG_OPT_HPP
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/common/transforms.h>

namespace opt{
    #define DIV_FACTOR		    1.4 // Division factor used for graduated non-convexity
    #define MAX_CORR_DIST	0.075// Maximum correspondence distance (also see comment of USE_ABSOLUTE_SCALE)
    #define ITERATION_NUMBER	100

    class Opt {
    public:
        double sigma = 2;
        float StartScale = 6;
        double par;
        int iteration_number_ = ITERATION_NUMBER;
        double max_corr_dist_ = MAX_CORR_DIST;
        double div_factor_ = DIV_FACTOR;
        bool decrease_mu_ = true;

        Eigen::Matrix4d
        GNCWelsch_Opt(pcl::PointCloud<pcl::PointXYZ>::Ptr &cor_SS, pcl::PointCloud<pcl::PointXYZ>::Ptr &cor_TT) {
            Eigen::Matrix4d trans, TransOutput_;
            trans.setIdentity();
            TransOutput_.setIdentity();
            pcl::PointCloud<pcl::PointXYZ>::Ptr cor_TT_copy(new pcl::PointCloud<pcl::PointXYZ>);
            pcl::copyPointCloud(*cor_TT, *cor_TT_copy);

            std::vector<double> s(cor_SS->size(), 1.0);
            Eigen::Matrix4d delta_last;
            double last_medistance = Medistance(cor_SS, cor_TT);
            delta_last.setIdentity();

            StartScale = float(10 * last_medistance);
            par = StartScale;
            div_factor_ = pow(par / max_corr_dist_, 1.0 / 20);

            for (int itr = 0; itr < iteration_number_; itr++) {
                if (decrease_mu_) {
                    if (itr % 2 == 0 && par > max_corr_dist_ && itr != 0) {
                        par /= div_factor_;
                    }
                }
                const int nvariable = 6;    // 3 for rotation and 3 for translation
                Eigen::MatrixXd JTJ(nvariable, nvariable);
                Eigen::MatrixXd D(nvariable, nvariable);
                Eigen::MatrixXd JTr(nvariable, 1);
                Eigen::MatrixXd J(nvariable, 1);
                Eigen::MatrixXd JT(nvariable, 1);

                double r;
                double r2 = 0.0;
                float temp;
                JT.setZero();
                JTJ.setZero();
                JTr.setZero();
                D.setZero();
                for (int c = 0; c < cor_TT_copy->size(); c++) {
                    Eigen::Vector3f p, q;
                    p = Eigen::Vector3f(cor_SS->points[c].x, cor_SS->points[c].y, cor_SS->points[c].z);
                    q = Eigen::Vector3f(cor_TT_copy->points[c].x, cor_TT_copy->points[c].y, cor_TT_copy->points[c].z);
                    int c2 = c;
                    Eigen::Vector3f rpq = p - q;
                    temp = pow(par * sigma, -2) * std::exp(-pow(rpq.norm(), 2) / (pow(par * sigma, 2)));

                    s[c2] = temp;
                    J.setZero();
                    J(1) = -q(2);
                    J(2) = q(1);
                    J(3) = -1;
                    r = rpq(0);
                    JTJ += J * J.transpose() * s[c2];
                    JTr += J * r * s[c2];
                    r2 += r * r * s[c2];
                    JT += J * sqrt(s[c2]);

                    J.setZero();
                    J(2) = -q(0);
                    J(0) = q(2);
                    J(4) = -1;
                    r = rpq(1);
                    JTJ += J * J.transpose() * s[c2];
                    JTr += J * r * s[c2];
                    r2 += r * r * s[c2];
                    JT += J * sqrt(s[c2]);

                    J.setZero();
                    J(0) = -q(1);
                    J(1) = q(0);
                    J(5) = -1;
                    r = rpq(2);
                    JTJ += J * J.transpose() * s[c2];
                    JTr += J * r * s[c2];
                    r2 += r * r * s[c2];
                    JT += J * sqrt(s[c2]);
                }

                Eigen::VectorXd result(nvariable, 1);
                //GM
                Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
                cod.compute(JTJ);
                result = cod.solve(-JTr);


                Eigen::Affine3d aff_mat;
                aff_mat.linear() = (Eigen::Matrix3d) Eigen::AngleAxisd(result(2), Eigen::Vector3d::UnitZ())
                                   * Eigen::AngleAxisd(result(1), Eigen::Vector3d::UnitY())
                                   * Eigen::AngleAxisd(result(0), Eigen::Vector3d::UnitX());

                aff_mat.translation() = Eigen::Vector3d(result(3), result(4), result(5));
                Eigen::Matrix4d delta = aff_mat.matrix().cast<double>();


                TransformPoints(cor_TT_copy, delta);
                trans = delta * trans;

                if ((delta_last - delta).norm() < 1e-12) {
                    break;
                }
                delta_last = delta;

            }
            TransOutput_ = trans * TransOutput_;

            return TransOutput_;
        };

    private:
        double Medistance(const pcl::PointCloud<pcl::PointXYZ>::Ptr &cor_SS,
                          const pcl::PointCloud<pcl::PointXYZ>::Ptr &cor_TT) {
            double dis_sum = 0;
            double med;
            for (int i = 0; i < cor_SS->size(); i++) {
                dis_sum += (Eigen::Vector3d(cor_SS->points[i].x, cor_SS->points[i].y, cor_SS->points[i].z) -
                            Eigen::Vector3d(cor_TT->points[i].x, cor_TT->points[i].y, cor_TT->points[i].z)).norm();
            }
            med = dis_sum / (cor_SS->size());
            return med;
        };

        void TransformPoints(pcl::PointCloud<pcl::PointXYZ>::Ptr &cor, Eigen::Matrix4d &Trans) {
            Trans(3, 3) = 1;
            pcl::transformPointCloud(*cor, *cor, Trans);
        };
    };

}

#endif //MICROG_OPT_HPP
