#ifndef ANDERSONACCELERATION_H_
#define ANDERSONACCELERATION_H_


#include <cassert>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>

typedef double Scalar;
template < int Rows, int Cols, int Options = (Eigen::ColMajor | Eigen::AutoAlign) >
using MatrixT = Eigen::Matrix<Scalar, Rows, Cols, Options>; ///< A typedef of the dense matrix of Eigen.

typedef MatrixT<Eigen::Dynamic, 1> VectorX;					///< A nd column vector.
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> MatrixXX;	///< A n by m matrix.

class AndersonAcceleration
{
public:
    AndersonAcceleration()
            :m_(-1), dim_(-1), iter_(-1), col_idx_(-1) {}

    void replace(const Scalar *u)
    {
        current_u_ = Eigen::Map<const VectorX>(u, dim_);
    }

    const VectorX& compute(const Scalar* g)
    {

        assert(iter_ >= 0);
        Eigen::Map<const VectorX> G(g, dim_);  //映射
        current_F_ = G - current_u_;

        if (iter_ == 0)
        {
            prev_dF_.col(0) = -current_F_;
            prev_dG_.col(0) = -G;
            current_u_ = G;
        }
        else
        {
            prev_dF_.col(col_idx_) += current_F_;
            prev_dG_.col(col_idx_) += G;

            Scalar eps = 1e-14;  //1e-14


            int m_k = std::min(m_, iter_);


            if (m_k == 1)
            {
                theta_(0) = 0;
                Scalar dF_sqrnorm = prev_dF_.col(col_idx_).squaredNorm();
                M_(0, 0) = dF_sqrnorm;
                Scalar dF_norm = std::sqrt(dF_sqrnorm);

                if (dF_norm > eps) {
                    theta_(0) = (prev_dF_.col(col_idx_) / dF_norm).dot(current_F_ / dF_norm);
                }
            }
            else
            {

                VectorX new_inner_prod = (prev_dF_.col(col_idx_).transpose() * prev_dF_.block(0, 0, dim_, m_k)).transpose();
                M_.block(col_idx_, 0, 1, m_k) = new_inner_prod.transpose();
                M_.block(0, col_idx_, m_k, 1) = new_inner_prod;

                MatrixXX M_aug(m_k + 1, m_k + 1);
                M_aug.block(0, 0, m_k, m_k) = M_;
                M_aug.row(m_k).setOnes();
                M_aug.col(m_k).setZero();

                MatrixXX prev_dF_aug(dim_+1,m_k+1);
                prev_dF_aug.block(0,0,dim_,m_k)=prev_dF_.block(0, 0, dim_, m_k);
                prev_dF_aug.row(dim_).setZero();
                prev_dF_aug.col(m_k).setZero();
                prev_dF_aug(dim_,m_k)=1;

                VectorX current_F_aug = VectorX::Ones(dim_ + 1);
                current_F_aug.head(dim_) = current_F_;

                cod_.compute(M_aug.block(0, 0, m_k+1, m_k+1));
                theta_.head(m_k) = cod_.solve(prev_dF_aug.transpose() * current_F_aug);
                double sum=theta_.head(m_k).sum();
                if(sum!=0) {
                    for (int ik = 0; ik < m_k; ik++) {
                        theta_(ik) /= sum;
                    }
                }


            }

            current_u_ = (G -prev_dG_.block(0, 0, dim_, m_k) * ((theta_.head(m_k).array()).matrix()));
            col_idx_ = (col_idx_ + 1) % m_;
            prev_dF_.col(col_idx_) = -current_F_;
            prev_dG_.col(col_idx_) = -G;
        }

        iter_++;
        return current_u_;
    }
    void reset(const Scalar *u)
    {
        iter_ = 0;
        col_idx_ = 0;
        current_u_ = Eigen::Map<const VectorX>(u, dim_);
    }

    // m: number of previous iterations used
    // d: dimension of variables
    // u0: initial variable values
    void init(int m, int d, const Scalar* u0)
    {
        assert(m > 0);
        m_ = m;
        dim_ = d;
        current_u_.resize(d);
        current_F_.resize(d);
        prev_dG_.resize(d, m);
        prev_dF_.resize(d, m);
        M_.resize(m, m);
        theta_.resize(m);
        dF_scale_.resize(m);
        current_u_ = Eigen::Map<const VectorX>(u0, d);
        iter_ = 0;
        col_idx_ = 0;
    }

private:

    VectorX current_u_;
    VectorX current_F_;
    MatrixXX prev_dG_;
    MatrixXX prev_dF_;
    MatrixXX M_;		// Normal equations matrix for the computing theta
    VectorX	theta_;	// theta value computed from normal equations
    VectorX dF_scale_;		// The scaling factor for each column of prev_dF
    Eigen::CompleteOrthogonalDecomposition<MatrixXX> cod_;


    int m_;		// Number of previous iterates used for Andreson Acceleration
    int dim_;	// Dimension of variables
    int iter_;	// Iteration count since initialization
    int col_idx_;	// Index for history matrix column to store the next value
    int m_k_;

    Eigen::Matrix4d current_T_;
    Eigen::Matrix4d current_F_T_;

    MatrixXX T_prev_dF_;
    MatrixXX T_prev_dG_;
};


#endif /* ANDERSONACCELERATION_H_ */
