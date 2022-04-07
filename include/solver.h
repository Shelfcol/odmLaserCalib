#ifndef CALIB_SOLVER_H
#define CALIB_SOLVER_H

#include <synchronizer.h>
#include <ceres/ceres.h>
#include <csm/math_utils.h>

template<typename T>
void ominus_d_template(const T x[3], T res[3]) {
	T c = ceres::cos(x[2]);
	T s = ceres::sin(x[2]);
	res[0] = -c*x[0]-s*x[1];
	res[1] =  s*x[0]-c*x[1];
	res[2] = -x[2];
}


template<typename T>
void oplus_d_template(const T x1[3], const T x2[3], T res[3]) {
	T c = ceres::cos(x1[2]);
	T s = ceres::sin(x1[2]);
	T x = x1[0] + c*x2[0] - s*x2[1];
	T y = x1[1] + s*x2[0] + c*x2[1];
 	T theta = x1[2] + x2[2];
	res[0] = x;
	res[1] = y;
	res[2] = theta;
}



class WheelLaserErr{
public:
  WheelLaserErr(const cSynchronizer::sync_data& sync_d):sync_d_(sync_d){}
  template<typename T>
  bool operator()(const T* const param, T* res) const { // r_L,r_R,b, lx,ly,l_theta
    T r_L = param[0];
    T r_R = param[1];
    T b = param[2];
    T lx = param[3];
    T ly = param[4];
    T l_theta = param[5];

    // calc wheel odom
    T J11 = r_L / T(2);
    T J12 = r_R / T(2);
    T J21 = - r_L / b;
    T J22 = r_R / b;

    T speed = J11 * sync_d_.velocity_left + J12 * sync_d_.velocity_right; // 速度
    T omega = J21 * sync_d_.velocity_left + J22 * sync_d_.velocity_right; // 角速度

    T o_theta = sync_d_.T * omega; // 旋转的角度

    T t1, t2;
    if (ceres::abs(o_theta) > 1e-12) // 有旋转
    {
      t1 = ceres::sin(o_theta) / o_theta;
      t2 = (T(1) - ceres::cos(o_theta)) / o_theta;
    }
    else
    {
      t1 = T(1);
      t2 = T(0);
    }

    T dx = t1 * speed * sync_d_.T;
    T dy = t2 * speed * sync_d_.T;
    T dtheta = o_theta;

    T l[3]={lx,ly,l_theta};
    T r[3]={dx,dy,dtheta};
    T res_right[3];
    ominus_d_template<T>(l, res_right);
    oplus_d_template<T>(res_right,r,res_right);
    oplus_d_template<T>(res_right,l,res_right);
    res[0]=T(sync_d_.scan_match_results[0])-res_right[0];
    res[1]=T(sync_d_.scan_match_results[1])-res_right[1];
    res[2]=T(sync_d_.scan_match_results[2])-res_right[2];

    return true;
  }
  static ceres::CostFunction* Create(const cSynchronizer::sync_data& sync_d){
    return  new ceres::AutoDiffCostFunction<WheelLaserErr,3,6>(new WheelLaserErr(sync_d));
  }
private:
  cSynchronizer::sync_data sync_d_;
};


class cSolver{

public:

  cSolver();

  struct solver_params {
    int mode;

    double max_cond_number;

    int outliers_iterations;
    double outliers_percentage;
  };

  struct calib_result {
    double radius_l, radius_r;
    double axle;

    /** Laser pose */
    double l[3]; // x,y,yaw
  };

  bool solve(const std::vector<cSynchronizer::sync_data> &calib_data, int mode, double max_cond_number, struct calib_result& res);
  bool solve(const std::vector<cSynchronizer::sync_data> &calib_data,struct calib_result& res);

  void calib(std::vector<cSynchronizer::sync_data> &calib_data, int outliers_iterations);

private:

  Eigen::VectorXd full_calibration_min(const Eigen::MatrixXd &M);

  Eigen::VectorXd numeric_calibration(const Eigen::MatrixXd &H);

  double calculate_error(const Eigen::VectorXd &x, const Eigen::MatrixXd &M);

  Eigen::VectorXd x_given_lambda(const Eigen::MatrixXd &M, const double &lambda, const Eigen::MatrixXd &W);

  void compute_disagreement(cSynchronizer::sync_data &calib_data, const struct calib_result &res);

  void estimate_noise(std::vector<cSynchronizer::sync_data> &calib_data, const struct calib_result &res,
                      double &std_x, double &std_y, double &std_th);

  double calculate_sd(const double array[], const int s, const int e);

  Eigen::MatrixXd compute_fim(const std::vector<cSynchronizer::sync_data> &calib_data, const struct calib_result &res,
                   const Eigen::Matrix3d &inf_sm);

};

#endif
