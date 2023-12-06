
#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <vector>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <string>
#include "matplotlibcpp.h"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
using namespace std;

// Aligned trajectories


struct Statistics {

public:
/// Root mean squared for the given values
double rmse = 0.0;

/// Mean of the given values
double mean = 0.0;

/// Median of the given values
double median = 0.0;

/// Standard deviation of given values
double std = 0.0;

/// Max of the given values
double max = 0.0;

/// Min of the given values
double min = 0.0;

/// 99th percentile
double ninetynine = 0.0;

/// Timestamp when these values occured at
std::vector<double> timestamps;

/// Values (e.g. error or nees at a given time)
std::vector<double> values;

/// Bound of these values (e.g. our expected covariance bound)
std::vector<double> values_bound;

/// Will calculate all values from our vectors of information
void calculate() {

    // Sort the data for easy finding of values
    std::vector<double> values_sorted = values;
    std::sort(values_sorted.begin(), values_sorted.end());

    // If we don't have any data, just return :(
    if (values_sorted.empty())
    return;

    // Now that its been sorted, can easily grab min and max
    min = values_sorted.at(0);
    max = values_sorted.at(values_sorted.size() - 1);

    // Compute median
    // ODD:  grab middle from the sorted vector
    // EVEN: average the middle two numbers
    if (values_sorted.size() == 1) {
    median = values_sorted.at(values_sorted.size() - 1);
    } else if (values_sorted.size() % 2 == 1) {
    median = values_sorted.at(values_sorted.size() / 2);
    } else if (values_sorted.size() > 1) {
    median = 0.5 * (values_sorted.at(values_sorted.size() / 2 - 1) + values_sorted.at(values_sorted.size() / 2));
    } else {
    median = 0.0;
    }

    // Compute mean and rmse
    mean = 0;
    for (size_t i = 0; i < values_sorted.size(); i++) {
    assert(!std::isnan(values_sorted.at(i)));
    mean += values_sorted.at(i);
    rmse += values_sorted.at(i) * values_sorted.at(i);
    }
    mean /= values_sorted.size();
    rmse = std::sqrt(rmse / values_sorted.size());

    // Using mean, compute standard deviation
    std = 0;
    for (size_t i = 0; i < values_sorted.size(); i++) {
    std += std::pow(values_sorted.at(i) - mean, 2);
    }
    std = std::sqrt(std / (values_sorted.size() - 1));

    // 99th percentile
    // TODO: is this correct?
    // TODO: http://sphweb.bumc.bu.edu/otlt/MPH-Modules/BS/BS704_Probability/BS704_Probability10.html
    ninetynine = mean + 2.326 * std;
}

/// Will clear any old values
void clear() {
    timestamps.clear();
    values.clear();
    values_bound.clear();
}
};


class ResultTrajectory{
public:
std::vector<Eigen::Matrix<double, 7, 1>> est_poses_aignedtoGT;
std::vector<Eigen::Matrix<double, 7, 1>> gt_poses_aignedtoEST;
std::vector<double> est_times, gt_times;
std::vector<Eigen::Matrix<double, 7, 1>> est_poses, gt_poses;

ResultTrajectory(std::string path_est, std::string path_gt, std::string alignment_method) {

    // Load from file
    load_data(path_est, est_times, est_poses);
    load_data(path_gt, gt_times, gt_poses);

    // Debug print amount
    // std::string base_filename1 = path_est.substr(path_est.find_last_of("/\\") + 1);
    // std::string base_filename2 = path_gt.substr(path_gt.find_last_of("/\\") + 1);
    // PRINT_DEBUG("[TRAJ]: loaded %d poses from %s\n",(int)est_times.size(),base_filename1.c_str());
    // PRINT_DEBUG("[TRAJ]: loaded %d poses from %s\n",(int)gt_times.size(),base_filename2.c_str());
    double len_gt = get_total_length(gt_poses);
    double len_est = get_total_length(est_poses);
    double ratio = len_est / len_gt;
    if (ratio > 1.1 || ratio < 0.9) {
        printf( "[TRAJ]: Trajectory is a bad ratio of %.2f length (est %.2f, gt %.2f)\n", ratio, len_est, len_gt);
        printf("[TRAJ]: %s\n", path_est.c_str());
    }

    // Intersect timestamps
    // perform_association(0, 0.02, est_times, gt_times, est_poses, gt_poses);
    // // Return failure if we didn't have any common timestamps
    // if (est_poses.size() < 3) {
    //     printf("[TRAJ]: unable to get enough common timestamps between trajectories.\n" );
    //     printf("[TRAJ]: does the estimated trajectory publish the rosbag timestamps??\n" );
    //     std::exit(EXIT_FAILURE);
    // }

    // Perform alignment of the trajectories
    Eigen::Matrix3d R_ESTtoGT, R_GTtoEST;
    Eigen::Vector3d t_ESTinGT, t_GTinEST;
    double s_ESTtoGT, s_GTtoEST;
    align_trajectory(est_poses, gt_poses, R_ESTtoGT, t_ESTinGT, s_ESTtoGT, alignment_method, 1);
    align_trajectory(gt_poses, est_poses, R_GTtoEST, t_GTinEST, s_GTtoEST, alignment_method, 1);

    // // Debug print to the user
    Eigen::Vector4d q_ESTtoGT = rot_2_quat(R_ESTtoGT);
    Eigen::Vector4d q_GTtoEST = rot_2_quat(R_GTtoEST);
    printf("[TRAJ]: q_ESTtoGT = %.3f, %.3f, %.3f, %.3f | p_ESTinGT = %.3f, %.3f, %.3f | s = %.2f\n", q_ESTtoGT(0), q_ESTtoGT(1),
                q_ESTtoGT(2), q_ESTtoGT(3), t_ESTinGT(0), t_ESTinGT(1), t_ESTinGT(2), s_ESTtoGT);
    // PRINT_DEBUG("[TRAJ]: q_GTtoEST = %.3f, %.3f, %.3f, %.3f | p_GTinEST = %.3f, %.3f, %.3f | s =
    //  %.2f\n",q_GTtoEST(0),q_GTtoEST(1),q_GTtoEST(2),q_GTtoEST(3),t_GTinEST(0),t_GTinEST(1),t_GTinEST(2),s_GTtoEST);

    // Finally lets calculate the aligned trajectories
    for (size_t i = 0; i < gt_times.size(); i++) {
      
      Eigen::Matrix<double, 7, 1> pose_ESTinGT, pose_GTinEST;
      pose_ESTinGT.block(0, 0, 3, 1) = s_ESTtoGT * R_ESTtoGT * est_poses.at(i).block(0, 0, 3, 1) + t_ESTinGT;
      pose_ESTinGT.block(3, 0, 4, 1) = quat_multiply(est_poses.at(i).block(3, 0, 4, 1), Inv(q_ESTtoGT));
      pose_GTinEST.block(0, 0, 3, 1) = s_GTtoEST * R_GTtoEST * gt_poses.at(i).block(0, 0, 3, 1) + t_GTinEST;
      pose_GTinEST.block(3, 0, 4, 1) = quat_multiply(gt_poses.at(i).block(3, 0, 4, 1), Inv(q_GTtoEST));
      est_poses_aignedtoGT.push_back(pose_ESTinGT);
      gt_poses_aignedtoEST.push_back(pose_GTinEST);
    }
}

Eigen::Matrix<double, 4, 1> Inv(Eigen::Matrix<double, 4, 1> q) {
  Eigen::Matrix<double, 4, 1> qinv;
  qinv.block(0, 0, 3, 1) = -q.block(0, 0, 3, 1);
  qinv(3, 0) = q(3, 0);
  return qinv;
}

void perform_association(double offset, double max_difference, std::vector<double> &est_times, std::vector<double> &gt_times,
                                     std::vector<Eigen::Matrix<double, 7, 1>> &est_poses,
                                     std::vector<Eigen::Matrix<double, 7, 1>> &gt_poses) {
  std::vector<Eigen::Matrix3d> est_covori, est_covpos, gt_covori, gt_covpos;
  perform_association(offset, max_difference, est_times, gt_times, est_poses, gt_poses, est_covori, est_covpos, gt_covori,
                                  gt_covpos);
}

void perform_association(double offset, double max_difference, std::vector<double> &est_times, std::vector<double> &gt_times,
                                     std::vector<Eigen::Matrix<double, 7, 1>> &est_poses,
                                     std::vector<Eigen::Matrix<double, 7, 1>> &gt_poses, std::vector<Eigen::Matrix3d> &est_covori,
                                     std::vector<Eigen::Matrix3d> &est_covpos, std::vector<Eigen::Matrix3d> &gt_covori,
                                     std::vector<Eigen::Matrix3d> &gt_covpos) {

  // Temp results which keeps only the matches
  std::vector<double> est_times_temp, gt_times_temp;
  std::vector<Eigen::Matrix<double, 7, 1>> est_poses_temp, gt_poses_temp;
  std::vector<Eigen::Matrix3d> est_covori_temp, est_covpos_temp, gt_covori_temp, gt_covpos_temp;

  // Iterator over gt (only ever increases to enforce injectivity of matches)
  size_t gt_pointer = 0;

  // Try to find closest GT pose for each estimate
  for (size_t i = 0; i < est_times.size(); i++) {

    // Default params
    double best_diff = max_difference;
    int best_gt_idx = -1;

    // Increment while too small and is not within our difference threshold
    while (gt_pointer < gt_times.size() && gt_times.at(gt_pointer) < (est_times.at(i) + offset) &&
           std::abs(gt_times.at(gt_pointer) - (est_times.at(i) + offset)) > max_difference) {
      gt_pointer++;
    }

    // If we are closer than max difference, see if we can do any better
    while (gt_pointer < gt_times.size() && std::abs(gt_times.at(gt_pointer) - (est_times.at(i) + offset)) <= max_difference) {
      // Break if we found a good match but are getting worse, we are done
      if (std::abs(gt_times.at(gt_pointer) - (est_times.at(i) + offset)) >= best_diff) {
        break;
      }
      // We have a closer match, save it and move on
      best_diff = std::abs(gt_times.at(gt_pointer) - (est_times.at(i) + offset));
      best_gt_idx = gt_pointer;
      gt_pointer++;
    }

    // Did we get a valid match
    if (best_gt_idx != -1) {

      // Save estimate and gt states for the match
      est_times_temp.push_back(gt_times.at(best_gt_idx));
      est_poses_temp.push_back(est_poses.at(i));
      gt_times_temp.push_back(gt_times.at(best_gt_idx));
      gt_poses_temp.push_back(gt_poses.at(best_gt_idx));

      // If we have covariance then also append it
      // If the groundtruth doesn't have covariance say it is 100% certain
      if (!est_covori.empty()) {
        assert(est_covori.size() == est_covpos.size());
        est_covori_temp.push_back(est_covori.at(i));
        est_covpos_temp.push_back(est_covpos.at(i));
        if (gt_covori.empty()) {
          gt_covori_temp.push_back(Eigen::Matrix3d::Zero());
          gt_covpos_temp.push_back(Eigen::Matrix3d::Zero());
        } else {
          assert(gt_covori.size() == gt_covpos.size());
          gt_covori_temp.push_back(gt_covori.at(best_gt_idx));
          gt_covpos_temp.push_back(gt_covpos.at(best_gt_idx));
        }
      }
    }
  }

  // Ensure that we have enough associations
  if (est_times.size() < 3) {
    printf( "[ALIGN]: Was unable to associate groundtruth and estimate trajectories\n" );
    printf("[ALIGN]: %d total matches....\n", (int)est_times.size());
    printf("[ALIGN]: Do the time of the files match??\n");
    return;
  }
  assert(est_times_temp.size() == gt_times_temp.size());
  // PRINT_DEBUG("[TRAJ]: %d est poses and %d gt poses => %d
  // matches\n",(int)est_times.size(),(int)gt_times.size(),(int)est_times_temp.size());

  // Overwrite with intersected values
  est_times = est_times_temp;
  est_poses = est_poses_temp;
  est_covori = est_covori_temp;
  est_covpos = est_covpos_temp;
  gt_times = gt_times_temp;
  gt_poses = gt_poses_temp;
  gt_covori = gt_covori_temp;
  gt_covpos = gt_covpos_temp;
}

void align_trajectory(const std::vector<Eigen::Matrix<double, 7, 1>> &traj_es,
                                       const std::vector<Eigen::Matrix<double, 7, 1>> &traj_gt, Eigen::Matrix3d &R, Eigen::Vector3d &t,
                                       double &s, std::string method, int n_aligned) {

  // Use the correct method
  if (method == "posyaw") {
    s = 1;
    align_posyaw(traj_es, traj_gt, R, t, n_aligned);
  // } else if (method == "posyawsingle") {
  //   s = 1;
  //   align_posyaw_single(traj_es, traj_gt, R, t);
  // } else if (method == "se3") {
  //   s = 1;
  //   align_se3(traj_es, traj_gt, R, t, n_aligned);
  // } else if (method == "se3single") {
  //   s = 1;
  //   align_se3_single(traj_es, traj_gt, R, t);
  // } else if (method == "sim3") {
  //   assert(n_aligned >= 2 || n_aligned == -1);
  //   align_sim3(traj_es, traj_gt, R, t, s, n_aligned);
  // } else if (method == "none") {
  //   s = 1;
  //   R.setIdentity();
  //   t.setZero();
  } else {
    printf("[ALIGN]: Invalid alignment method!\n" );
    printf("[ALIGN]: Possible options: posyaw, posyawsingle, se3, se3single, sim3, none\n");
    std::exit(EXIT_FAILURE);
  }
}

void align_sim3(const std::vector<Eigen::Matrix<double, 7, 1>> &traj_es,
                                 const std::vector<Eigen::Matrix<double, 7, 1>> &traj_gt, Eigen::Matrix3d &R, Eigen::Vector3d &t, double &s,
                                 int n_aligned) {

  // Need to have more than two to get
  assert(n_aligned >= 2 || n_aligned == -1);

  // Get just position vectors
  assert(!traj_es.empty());
  std::vector<Eigen::Vector3d> pos_est, pos_gt;
  for (size_t i = 0; i < traj_es.size() && i < traj_gt.size(); i++) {
    pos_est.push_back(traj_es.at(i).block(0, 0, 3, 1));
    pos_gt.push_back(traj_gt.at(i).block(0, 0, 3, 1));
  }

  // Align using the method of Umeyama
  align_umeyama(pos_est, pos_gt, R, t, s, false, false);
}

void align_se3_single(const std::vector<Eigen::Matrix<double, 7, 1>> &traj_es,
                                       const std::vector<Eigen::Matrix<double, 7, 1>> &traj_gt, Eigen::Matrix3d &R, Eigen::Vector3d &t) {

  // Get first ever poses
  Eigen::Vector4d q_es_0 = traj_es.at(0).block(3, 0, 4, 1);
  Eigen::Vector3d p_es_0 = traj_es.at(0).block(0, 0, 3, 1);

  Eigen::Vector4d q_gt_0 = traj_gt.at(0).block(3, 0, 4, 1);
  Eigen::Vector3d p_gt_0 = traj_gt.at(0).block(0, 0, 3, 1);

  // Get rotations from IMU frame to World (note JPL!)
  Eigen::Matrix3d g_rot = quat_2_Rot(q_gt_0).transpose();
  Eigen::Matrix3d est_rot = quat_2_Rot(q_es_0).transpose();

  R.noalias() = g_rot * est_rot.transpose();
  t.noalias() = p_gt_0 - R * p_es_0;
}

void align_se3(const std::vector<Eigen::Matrix<double, 7, 1>> &traj_es,
                                const std::vector<Eigen::Matrix<double, 7, 1>> &traj_gt, Eigen::Matrix3d &R, Eigen::Vector3d &t,
                                int n_aligned) {

  // If we only have one, just use the single alignment
  if (n_aligned == 1) {
    align_se3_single(traj_es, traj_gt, R, t);
  } else {

    // Get just position vectors
    assert(!traj_es.empty());
    std::vector<Eigen::Vector3d> pos_est, pos_gt;
    for (size_t i = 0; i < traj_es.size() && i < traj_gt.size(); i++) {
      pos_est.push_back(traj_es.at(i).block(0, 0, 3, 1));
      pos_gt.push_back(traj_gt.at(i).block(0, 0, 3, 1));
    }

    // Align using the method of Umeyama
    double s;
    align_umeyama(pos_est, pos_gt, R, t, s, true, false);
  }
}

void align_posyaw(const std::vector<Eigen::Matrix<double, 7, 1>> &traj_es,
                                   const std::vector<Eigen::Matrix<double, 7, 1>> &traj_gt, Eigen::Matrix3d &R, Eigen::Vector3d &t,
                                   int n_aligned) {

  // If we only have one, just use the single alignment
  if (n_aligned == 1) {
    align_posyaw_single(traj_es, traj_gt, R, t);
  } else {

    // Get just position vectors
    assert(!traj_es.empty());
    std::vector<Eigen::Vector3d> pos_est, pos_gt;
    for (size_t i = 0; i < traj_es.size() && i < traj_gt.size(); i++) {
      pos_est.push_back(traj_es.at(i).block(0, 0, 3, 1));
      pos_gt.push_back(traj_gt.at(i).block(0, 0, 3, 1));
    }

    // Align using the method of Umeyama
    double s;
    align_umeyama(pos_est, pos_gt, R, t, s, true, true);
    assert(s == 1);
  }
}

void align_umeyama(const std::vector<Eigen::Matrix<double, 3, 1>> &data, const std::vector<Eigen::Matrix<double, 3, 1>> &model,
                               Eigen::Matrix<double, 3, 3> &R, Eigen::Matrix<double, 3, 1> &t, double &s, bool known_scale, bool yaw_only) {

  assert(model.size() == data.size());

  // Substract mean of each trajectory
  Eigen::Matrix<double, 3, 1> mu_M = get_mean(model);
  Eigen::Matrix<double, 3, 1> mu_D = get_mean(data);
  std::vector<Eigen::Matrix<double, 3, 1>> model_zerocentered, data_zerocentered;
  for (size_t i = 0; i < model.size(); i++) {
    model_zerocentered.push_back(model[i] - mu_M);
    data_zerocentered.push_back(data[i] - mu_D);
  }

  // Get correlation matrix
  double n = model.size();
  Eigen::Matrix<double, 3, 3> C = Eigen::Matrix<double, 3, 3>::Zero();
  for (size_t i = 0; i < model_zerocentered.size(); i++) {
    C.noalias() += model_zerocentered[i] * data_zerocentered[i].transpose();
  }
  C *= 1.0 / n;

  // Get data sigma
  double sigma2 = 0;
  for (size_t i = 0; i < data_zerocentered.size(); i++) {
    sigma2 += data_zerocentered[i].dot(data_zerocentered[i]);
  }
  sigma2 *= 1.0 / n;

  // SVD decomposition
  Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd(C, Eigen::ComputeFullV | Eigen::ComputeFullU);

  Eigen::Matrix<double, 3, 3> U_svd = svd.matrixU();
  Eigen::Matrix<double, 3, 1> D_svd = svd.singularValues();
  Eigen::Matrix<double, 3, 3> V_svd = svd.matrixV();

  Eigen::Matrix<double, 3, 3> S = Eigen::Matrix<double, 3, 3>::Identity();
  if (U_svd.determinant() * V_svd.determinant() < 0) {
    S(2, 2) = -1;
  }

  // If only yaw, use that specific solver (optimizes over yaw angle)
  // Else get best full 3 dof rotation
  if (yaw_only) {
    Eigen::Matrix<double, 3, 3> rot_C = n * C.transpose();
    double theta = get_best_yaw(rot_C);
    R = rot_z(theta);
  } else {
    R.noalias() = U_svd * S * V_svd.transpose();
  }

  // If known scale, fix it
  if (known_scale) {
    s = 1;
  } else {
    // Get best scale
    s = 1.0 / sigma2 * (D_svd.asDiagonal() * S).trace();
  }

  // Get best translation
  t.noalias() = mu_M - s * R * mu_D;
}

void align_posyaw_single(const std::vector<Eigen::Matrix<double, 7, 1>> &traj_es,
                                          const std::vector<Eigen::Matrix<double, 7, 1>> &traj_gt, Eigen::Matrix3d &R, Eigen::Vector3d &t) {

  // Get first ever poses
  Eigen::Vector4d q_es_0 = traj_es.at(0).block(3, 0, 4, 1);
  Eigen::Vector3d p_es_0 = traj_es.at(0).block(0, 0, 3, 1);

  Eigen::Vector4d q_gt_0 = traj_gt.at(0).block(3, 0, 4, 1);
  Eigen::Vector3d p_gt_0 = traj_gt.at(0).block(0, 0, 3, 1);

  // Get rotations from IMU frame to World (note JPL!)
  Eigen::Matrix3d g_rot = quat_2_Rot(q_gt_0).transpose();
  Eigen::Matrix3d est_rot = quat_2_Rot(q_es_0).transpose();

  // Data matrix for the Frobenius problem
  Eigen::Matrix3d C_R = est_rot * g_rot.transpose();

  // Recover yaw
  double theta = get_best_yaw(C_R);

  // Compute rotation
  R = rot_z(theta);

  // Compute translation
  t.noalias() = p_gt_0 - R * p_es_0;
}

static double get_best_yaw(const Eigen::Matrix<double, 3, 3> &C) {
    double A = C(0, 1) - C(1, 0);
    double B = C(0, 0) + C(1, 1);
    // return M_PI_2 - atan2(B, A);
    return atan2(A, B);
}

static  Eigen::Matrix<double, 3, 1> get_mean(const std::vector<Eigen::Matrix<double, 3, 1>> &data) {
    Eigen::Matrix<double, 3, 1> mean = Eigen::Matrix<double, 3, 1>::Zero();
    for (size_t i = 0; i < data.size(); i++) {
      mean.noalias() += data[i];
    }
    mean /= data.size();
    return mean;
}

Eigen::Matrix<double, 3, 3> rot_z(double t) {
  Eigen::Matrix<double, 3, 3> r;
  double ct = cos(t);
  double st = sin(t);
  r << ct, -st, 0.0, st, ct, 0.0, 0.0, 0.0, 1.0;
  return r;
}

Eigen::Matrix<double, 4, 1> quat_multiply(const Eigen::Matrix<double, 4, 1> &q, const Eigen::Matrix<double, 4, 1> &p) {
  Eigen::Matrix<double, 4, 1> q_t;
  Eigen::Matrix<double, 4, 4> Qm;
  // create big L matrix
  Qm.block(0, 0, 3, 3) = q(3, 0) * Eigen::MatrixXd::Identity(3, 3) - skew_x(q.block(0, 0, 3, 1));
  Qm.block(0, 3, 3, 1) = q.block(0, 0, 3, 1);
  Qm.block(3, 0, 1, 3) = -q.block(0, 0, 3, 1).transpose();
  Qm(3, 3) = q(3, 0);
  q_t = Qm * p;
  // ensure unique by forcing q_4 to be >0
  if (q_t(3, 0) < 0) {
    q_t *= -1;
  }
  // normalize and return
  return q_t / q_t.norm();
}

double get_total_length(const std::vector<Eigen::Matrix<double, 7, 1>> &poses) {

  // Loop through every pose and append its segment
  double distance = 0.0;
  for (size_t i = 1; i < poses.size(); i++) {
    distance += (poses[i].block(0, 0, 3, 1) - poses[i - 1].block(0, 0, 3, 1)).norm();
  }

  // return the distance
  return distance;
}

Eigen::Matrix<double, 3, 1> log_so3(const Eigen::Matrix<double, 3, 3> &R) 
{
    // note switch to base 1
    double R11 = R(0, 0), R12 = R(0, 1), R13 = R(0, 2);
    double R21 = R(1, 0), R22 = R(1, 1), R23 = R(1, 2);
    double R31 = R(2, 0), R32 = R(2, 1), R33 = R(2, 2);
    // Get trace(R)
    const double tr = R.trace();
    Eigen::Vector3d omega;
    // when trace == -1, i.e., when theta = +-pi, +-3pi, +-5pi, etc.
    // we do something special
    if (tr + 1.0 < 1e-10) {
        if (std::abs(R33 + 1.0) > 1e-5)
        omega = (M_PI / sqrt(2.0 + 2.0 * R33)) * Eigen::Vector3d(R13, R23, 1.0 + R33);
        else if (std::abs(R22 + 1.0) > 1e-5)
        omega = (M_PI / sqrt(2.0 + 2.0 * R22)) * Eigen::Vector3d(R12, 1.0 + R22, R32);
        else
        // if(std::abs(R.r1_.x()+1.0) > 1e-5)  This is implicit
        omega = (M_PI / sqrt(2.0 + 2.0 * R11)) * Eigen::Vector3d(1.0 + R11, R21, R31);
    } else {
        double magnitude;
        const double tr_3 = tr - 3.0; // always negative
        if (tr_3 < -1e-7) {
        double theta = acos((tr - 1.0) / 2.0);
        magnitude = theta / (2.0 * sin(theta));
        } else {
        // when theta near 0, +-2pi, +-4pi, etc. (trace near 3.0)
        // use Taylor expansion: theta \approx 1/2-(t-3)/12 + O((t-3)^2)
        // see https://github.com/borglab/gtsam/issues/746 for details
        magnitude = 0.5 - tr_3 / 12.0;
        }
        omega = magnitude * Eigen::Vector3d(R32 - R23, R13 - R31, R21 - R12);
    }

    return omega;
}

Eigen::Matrix<double, 3, 3> skew_x(const Eigen::Matrix<double, 3, 1> &w) {
Eigen::Matrix<double, 3, 3> w_x;
w_x << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;
return w_x;
}

Eigen::Matrix<double, 4, 1> rot_2_quat(const Eigen::Matrix<double, 3, 3> &rot) {
  Eigen::Matrix<double, 4, 1> q;
  double T = rot.trace();
  if ((rot(0, 0) >= T) && (rot(0, 0) >= rot(1, 1)) && (rot(0, 0) >= rot(2, 2))) {
    q(0) = sqrt((1 + (2 * rot(0, 0)) - T) / 4);
    q(1) = (1 / (4 * q(0))) * (rot(0, 1) + rot(1, 0));
    q(2) = (1 / (4 * q(0))) * (rot(0, 2) + rot(2, 0));
    q(3) = (1 / (4 * q(0))) * (rot(1, 2) - rot(2, 1));

  } else if ((rot(1, 1) >= T) && (rot(1, 1) >= rot(0, 0)) && (rot(1, 1) >= rot(2, 2))) {
    q(1) = sqrt((1 + (2 * rot(1, 1)) - T) / 4);
    q(0) = (1 / (4 * q(1))) * (rot(0, 1) + rot(1, 0));
    q(2) = (1 / (4 * q(1))) * (rot(1, 2) + rot(2, 1));
    q(3) = (1 / (4 * q(1))) * (rot(2, 0) - rot(0, 2));
  } else if ((rot(2, 2) >= T) && (rot(2, 2) >= rot(0, 0)) && (rot(2, 2) >= rot(1, 1))) {
    q(2) = sqrt((1 + (2 * rot(2, 2)) - T) / 4);
    q(0) = (1 / (4 * q(2))) * (rot(0, 2) + rot(2, 0));
    q(1) = (1 / (4 * q(2))) * (rot(1, 2) + rot(2, 1));
    q(3) = (1 / (4 * q(2))) * (rot(0, 1) - rot(1, 0));
  } else {
    q(3) = sqrt((1 + T) / 4);
    q(0) = (1 / (4 * q(3))) * (rot(1, 2) - rot(2, 1));
    q(1) = (1 / (4 * q(3))) * (rot(2, 0) - rot(0, 2));
    q(2) = (1 / (4 * q(3))) * (rot(0, 1) - rot(1, 0));
  }
  if (q(3) < 0) {
    q = -q;
  }
  // normalize and return
  q = q / (q.norm());
  return q;
}

Eigen::Matrix<double, 3, 3> quat_2_Rot(const Eigen::Matrix<double, 4, 1> &q) 
{
    Eigen::Matrix<double, 3, 3> q_x = skew_x(q.block(0, 0, 3, 1));
    Eigen::MatrixXd Rot = (2 * std::pow(q(3, 0), 2) - 1) * Eigen::MatrixXd::Identity(3, 3) - 2 * q(3, 0) * q_x +
                            2 * q.block(0, 0, 3, 1) * (q.block(0, 0, 3, 1).transpose());
    return Rot;
}

Eigen::Matrix4d Inv_se3(const Eigen::Matrix4d &T) {
    Eigen::Matrix4d Tinv = Eigen::Matrix4d::Identity();
    Tinv.block(0, 0, 3, 3) = T.block(0, 0, 3, 3).transpose();
    Tinv.block(0, 3, 3, 1) = -Tinv.block(0, 0, 3, 3) * T.block(0, 3, 3, 1);
    return Tinv;
}

std::vector<int> compute_comparison_indices_length(std::vector<double> &distances, double distance, double max_dist_diff) {

    // Vector of end ids for our pose indexes
    std::vector<int> comparisons;

    // Loop through each pose in our trajectory (i.e. our distance vector generated from the trajectory).
    for (size_t idx = 0; idx < distances.size(); idx++) {

    // Loop through and find the pose that minimized the difference between
    // The desired trajectory distance and our current trajectory distance
    double distance_startpose = distances.at(idx);
    int best_idx = -1;
    double best_error = max_dist_diff;
    for (size_t i = idx; i < distances.size(); i++) {
        if (std::abs(distances.at(i) - (distance_startpose + distance)) < best_error) {
        best_idx = i;
        best_error = std::abs(distances.at(i) - (distance_startpose + distance));
        }
    }
    // If we have an end id that reached this trajectory distance then add it!
    // Else this isn't a valid segment, thus we shouldn't add it (we will try again at the next pose)
    // NOTE: just because we searched through all poses and didn't find a close one doesn't mean we have ended
    // NOTE: this could happen if there is a gap in the groundtruth poses and we just couldn't find a pose with low error
    comparisons.push_back(best_idx);
    }

    // Finally return the ids for each starting pose that have this distance
    return comparisons;
}

void calculate_ate(Statistics &error_ori, Statistics &error_pos) 
{
    // Clear any old data
    error_ori.clear();
    error_pos.clear();
    
    // Calculate the position and orientation error at every timestep
    for (size_t i = 0; i < est_poses_aignedtoGT.size(); i++) {

        // Calculate orientation error
        Eigen::Matrix3d e_R = quat_2_Rot(est_poses_aignedtoGT.at(i).block(3, 0, 4, 1)).transpose() *
                            quat_2_Rot(gt_poses.at(i).block(3, 0, 4, 1));
        double ori_err = 180.0 / M_PI * log_so3(e_R).norm();

        // Calculate position error
        double pos_err = (gt_poses.at(i).block(0, 0, 3, 1) - est_poses_aignedtoGT.at(i).block(0, 0, 3, 1)).norm();

        // Append this error!
        error_ori.timestamps.push_back(est_times.at(i));
        error_ori.values.push_back(ori_err);
        error_pos.timestamps.push_back(est_times.at(i));
        error_pos.values.push_back(pos_err);
    }
    // Update stat information
    error_ori.calculate();
    error_pos.calculate();
}

void calculate_rpe(const std::vector<double> &segment_lengths,
                                    std::map<double, std::pair<Statistics, Statistics>> &error_rpe) {

    // Distance at each point along the trajectory
    double average_pos_difference = 0;
    std::vector<double> accum_distances(gt_poses.size());
    accum_distances[0] = 0;
    for (size_t i = 1; i < gt_poses.size(); i++) {
        double pos_diff = (gt_poses[i].block(0, 0, 3, 1) - gt_poses[i - 1].block(0, 0, 3, 1)).norm();
        accum_distances[i] = accum_distances[i - 1] + pos_diff;
        average_pos_difference += pos_diff;
    }

    average_pos_difference /= gt_poses.size();

    // Warn if large pos difference
    double max_dist_diff = 0.5;
    if (average_pos_difference > max_dist_diff) {
        // PRINT_WARNING(YELLOW "[COMP]: average groundtruth position difference %.2f > %.2f is too large\n" RESET, average_pos_difference,
        //             max_dist_diff);
        // PRINT_WARNING(YELLOW "[COMP]: this will prevent the RPE from finding valid trajectory segments!!!\n" RESET);
        // PRINT_WARNING(YELLOW
        //             "[COMP]: the recommendation is to use a higher frequency groundtruth, or relax this trajectory segment logic...\n" RESET);
        cout << "average groundtruth position difference"<< endl;
    }

    // Loop through each segment length
    for (const double &distance : segment_lengths) {

        // Our stats for this length
        Statistics error_ori, error_pos;
        // Get end of subtrajectories for each possible starting point
        // NOTE: is there a better way to select which end pos is a valid segments that are of the correct lenght?
        // NOTE: right now this allows for longer segments to have bigger error in their start-end distance vs the desired segment length
        // std::vector<int> comparisons = compute_comparison_indices_length(accum_distances, distance, 0.1*distance);
        std::vector<int> comparisons = compute_comparison_indices_length(accum_distances, distance, max_dist_diff);
        assert(comparisons.size() == gt_poses.size());
        // Loop through each relative comparison
        for (size_t id_start = 0; id_start < comparisons.size(); id_start++) {
          // Get the end id (skip if we couldn't find an end)
          int id_end = comparisons[id_start];
          if (id_end == -1)
              continue;
          //===============================================================================
          // Get T I1 to world EST at beginning of subtrajectory (at state idx)
          Eigen::Matrix4d T_c1 = Eigen::Matrix4d::Identity();
          T_c1.block(0, 0, 3, 3) = quat_2_Rot(est_poses_aignedtoGT.at(id_start).block(3, 0, 4, 1)).transpose();
          T_c1.block(0, 3, 3, 1) = est_poses_aignedtoGT.at(id_start).block(0, 0, 3, 1);

          // Get T I2 to world EST at end of subtrajectory starting (at state comparisons[idx])
          Eigen::Matrix4d T_c2 = Eigen::Matrix4d::Identity();
          T_c2.block(0, 0, 3, 3) = quat_2_Rot(est_poses_aignedtoGT.at(id_end).block(3, 0, 4, 1)).transpose();
          T_c2.block(0, 3, 3, 1) = est_poses_aignedtoGT.at(id_end).block(0, 0, 3, 1);

          // Get T I2 to I1 EST
          Eigen::Matrix4d T_c1_c2 = Inv_se3(T_c1) * T_c2;

          //===============================================================================
          // Get T I1 to world GT at beginning of subtrajectory (at state idx)
          Eigen::Matrix4d T_m1 = Eigen::Matrix4d::Identity();
          T_m1.block(0, 0, 3, 3) = quat_2_Rot(gt_poses.at(id_start).block(3, 0, 4, 1)).transpose();
          T_m1.block(0, 3, 3, 1) = gt_poses.at(id_start).block(0, 0, 3, 1);

          // Get T I2 to world GT at end of subtrajectory starting (at state comparisons[idx])
          Eigen::Matrix4d T_m2 = Eigen::Matrix4d::Identity();
          T_m2.block(0, 0, 3, 3) = quat_2_Rot(gt_poses.at(id_end).block(3, 0, 4, 1)).transpose();
          T_m2.block(0, 3, 3, 1) = gt_poses.at(id_end).block(0, 0, 3, 1);

          // Get T I2 to I1 GT
          Eigen::Matrix4d T_m1_m2 = Inv_se3(T_m1) * T_m2;

          //===============================================================================
          // Compute error transform between EST and GT start-end transform
          Eigen::Matrix4d T_error_in_c2 = Inv_se3(T_m1_m2) * T_c1_c2;

          Eigen::Matrix4d T_c2_rot = Eigen::Matrix4d::Identity();
          T_c2_rot.block(0, 0, 3, 3) = T_c2.block(0, 0, 3, 3);
          Eigen::Matrix4d T_c2_rot_inv = Eigen::Matrix4d::Identity();
          T_c2_rot_inv.block(0, 0, 3, 3) = T_c2.block(0, 0, 3, 3).transpose();

          // Rotate rotation so that rotation error is in the global frame (allows us to look at yaw error)
          Eigen::Matrix4d T_error_in_w = T_c2_rot * T_error_in_c2 * T_c2_rot_inv;

          //===============================================================================
          // Compute error for position
          error_pos.timestamps.push_back(est_times.at(id_start));
          error_pos.values.push_back(T_error_in_w.block(0, 3, 3, 1).norm());

          // Calculate orientation error
          double ori_err = 180.0 / M_PI * log_so3(T_error_in_w.block(0, 0, 3, 3)).norm();
          error_ori.timestamps.push_back(est_times.at(id_start));
          error_ori.values.push_back(ori_err);
        }

        // Update stat information
        error_ori.calculate();
        error_pos.calculate();
        error_rpe.insert({distance, {error_ori, error_pos}});

    }

}


void load_data(std::string path_traj, std::vector<double> &times, std::vector<Eigen::Matrix<double, 7, 1>> &poses) {

  // Try to open our trajectory file
  std::ifstream file(path_traj);
  if (!file.is_open()) {
    printf("[LOAD]: Unable to open trajectory file...\n" );
    printf("[LOAD]: %s\n" , path_traj.c_str());
    std::exit(EXIT_FAILURE);
  }

  // Loop through each line of this file
  std::string current_line;
  while (std::getline(file, current_line)) {
    // Skip if we start with a comment
    if (!current_line.find("#"))
      continue;
    // Loop variables
    int i = 0;
    std::istringstream s(current_line);
    std::string field;
    Eigen::Matrix<double, 8, 1> data;

    // Loop through this line (timestamp(s) tx ty tz qx qy qz qw Pr11 Pr12 Pr13 Pr22 Pr23 Pr33 Pt11 Pt12 Pt13 Pt22 Pt23 Pt33)
    while (std::getline(s, field, ' ')) {
      // Skip if empty
      if (field.empty() || i >= data.rows())
        continue;
      // save the data to our vector
      data(i) = std::atof(field.c_str());
      i++;
    }

      // time and pose
      times.push_back(data(0));
      poses.push_back(data.block(1, 0, 7, 1));

  }
  // Finally close the file
  file.close();

  // Error if we don't have any data
  if (times.empty()) {
    printf("[LOAD]: Could not parse any data from the file!!\n" );
    printf("[LOAD]: %s\n" , path_traj.c_str());
    std::exit(EXIT_FAILURE);
  }

  // Assert that they are all equal
  if (times.size() != poses.size()) {
    printf( "[LOAD]: Parsing error, pose and timestamps do not match!!\n" );
    printf("[LOAD]: %s\n", path_traj.c_str());
    std::exit(EXIT_FAILURE);
  }
}


};

void load_data(std::string path_traj, std::vector<double> &times, std::vector<Eigen::Matrix<double, 7, 1>> &poses) {

  // Try to open our trajectory file
  std::ifstream file(path_traj);
  if (!file.is_open()) {
    printf("[LOAD]: Unable to open trajectory file...\n" );
    printf("[LOAD]: %s\n" , path_traj.c_str());
    std::exit(EXIT_FAILURE);
  }

  // Loop through each line of this file
  std::string current_line;
  while (std::getline(file, current_line)) {
    // Skip if we start with a comment
    if (!current_line.find("#"))
      continue;
    // Loop variables
    int i = 0;
    std::istringstream s(current_line);
    std::string field;
    Eigen::Matrix<double, 8, 1> data;

    // Loop through this line (timestamp(s) tx ty tz qx qy qz qw Pr11 Pr12 Pr13 Pr22 Pr23 Pr33 Pt11 Pt12 Pt13 Pt22 Pt23 Pt33)
    while (std::getline(s, field, ' ')) {
      // Skip if empty
      if (field.empty() || i >= data.rows())
        continue;
      // save the data to our vector
      data(i) = std::atof(field.c_str());
      i++;
    }

      // time and pose
      times.push_back(data(0));
      poses.push_back(data.block(1, 0, 7, 1));

  }
  // Finally close the file
  file.close();

  // Error if we don't have any data
  if (times.empty()) {
    printf("[LOAD]: Could not parse any data from the file!!\n" );
    printf("[LOAD]: %s\n" , path_traj.c_str());
    std::exit(EXIT_FAILURE);
  }

  // Assert that they are all equal
  if (times.size() != poses.size()) {
    printf( "[LOAD]: Parsing error, pose and timestamps do not match!!\n" );
    printf("[LOAD]: %s\n", path_traj.c_str());
    std::exit(EXIT_FAILURE);
  }
}
double get_total_length(const std::vector<Eigen::Matrix<double, 7, 1>> &poses) {

  // Loop through every pose and append its segment
  double distance = 0.0;
  for (size_t i = 1; i < poses.size(); i++) {
    distance += (poses[i].block(0, 0, 3, 1) - poses[i - 1].block(0, 0, 3, 1)).norm();
  }

  // return the distance
  return distance;
}


int main() {


  // List the groundtruth files in this folder
  std::string path_gts("/home/zhangyanyu/open_vins/traj_data/");
  std::vector<boost::filesystem::path> path_groundtruths;
  for (const auto &p : boost::filesystem::recursive_directory_iterator(path_gts)) {
    if (p.path().extension() == ".txt") {
      path_groundtruths.push_back(p.path());
    }
  }
  std::sort(path_groundtruths.begin(), path_groundtruths.end());

  // Try to load our paths
  for (size_t i = 0; i < path_groundtruths.size(); i++) {
    // Load it!
    std::vector<double> times;
    std::vector<Eigen::Matrix<double, 7, 1>> poses;


    load_data(path_groundtruths.at(i).string(), times, poses);
    // Print its length and stats
    double length = get_total_length(poses);
    printf("[COMP]: %d poses in %s => length of %.2f meters\n", (int)times.size(), path_groundtruths.at(i).filename().c_str(), length);
  }

  // Get the algorithms we will process
  // Also create empty statistic objects for each of our datasets
  std::string path_algos("/home/zhangyanyu/open_vins/alg/");
  std::vector<boost::filesystem::path> path_algorithms;
  for (const auto &entry : boost::filesystem::directory_iterator(path_algos)) {
    if (boost::filesystem::is_directory(entry)) {
      path_algorithms.push_back(entry.path());
    }
  }
  std::sort(path_algorithms.begin(), path_algorithms.end());

  //===============================================================================
  //===============================================================================
  //===============================================================================

  // ATE summery information
  std::map<std::string, std::vector<std::pair<Statistics, Statistics>>> algo_ate;
//   std::map<std::string, std::vector<std::pair<Statistics, Statistics>>> algo_ate_2d;
  for (const auto &p : path_algorithms) {
    std::vector<std::pair<Statistics, Statistics>> temp;
    for (size_t i = 0; i < path_groundtruths.size(); i++) {
      temp.push_back({Statistics(), Statistics()});
    }
    algo_ate.insert({p.filename().string(), temp});
    // algo_ate_2d.insert({p.filename().string(), temp});
  }

  // Relative pose error segment lengths
  std::vector<double> segments = {8.0, 16.0, 24.0, 32.0, 40.0, 48.0};
  // std::vector<double> segments = {10.0, 25.0, 50.0, 75.0, 120.0};
  // std::vector<double> segments = {40.0, 80.0, 120.0, 160.0, 200.0, 240.0};

  // The overall RPE error calculation for each algorithm type
  std::map<std::string, std::map<double, std::pair<Statistics, Statistics>>> algo_rpe;
  for (const auto &p : path_algorithms) {
    std::map<double, std::pair<Statistics, Statistics>> temp;
    for (const auto &len : segments) {
      temp.insert({len, {Statistics(), Statistics()}});
    }
    algo_rpe.insert({p.filename().string(), temp});
  }
  //===============================================================================
  //===============================================================================
  //===============================================================================
  // Loop through each algorithm type
  for (size_t i = 0; i < path_algorithms.size(); i++) {
    // Debug print
    printf("======================================\n");
    printf("[COMP]: processing %s algorithm\n", path_algorithms.at(i).filename().c_str());

    // Get the list of datasets this algorithm records
    std::map<std::string, boost::filesystem::path> path_algo_datasets;
    for (auto &entry : boost::filesystem::directory_iterator(path_algorithms.at(i))) {
      if (boost::filesystem::is_directory(entry)) {
        path_algo_datasets.insert({entry.path().filename().string(), entry.path()});
      }
    }
    // Loop through our list of groundtruth datasets, and see if we have it
    for (size_t j = 0; j < path_groundtruths.size(); j++) {

      // Check if we have runs for this dataset
      if (path_algo_datasets.find(path_groundtruths.at(j).stem().string()) == path_algo_datasets.end()) {
        printf( "[COMP]: %s dataset does not have any runs for %s!!!!!\n", path_algorithms.at(i).filename().c_str(),
                    path_groundtruths.at(j).stem().c_str());
        continue;
      }
      // Debug print
      printf("[COMP]: processing %s algorithm => %s dataset\n", path_algorithms.at(i).filename().c_str(),
                  path_groundtruths.at(j).stem().c_str());
      // Errors for this specific dataset (i.e. our averages over the total runs)
      Statistics ate_dataset_ori, ate_dataset_pos;
      Statistics ate_2d_dataset_ori, ate_2d_dataset_pos;
      std::map<double, std::pair<Statistics, Statistics>> rpe_dataset;
      for (const auto &len : segments) {
        rpe_dataset.insert({len, {Statistics(), Statistics()}});
      }
      // Loop though the different runs for this dataset
      std::vector<std::string> file_paths;
      for (auto &entry : boost::filesystem::directory_iterator(path_algo_datasets.at(path_groundtruths.at(j).stem().string()))) {
        if (entry.path().extension() != ".txt")
          continue;
        file_paths.push_back(entry.path().string());
      }
      std::sort(file_paths.begin(), file_paths.end());

      // Now loop through the sorted vector
      for (auto &path_esttxt : file_paths) {
        // Our paths
        std::string dataset = path_groundtruths.at(j).stem().string();
        std::string path_gttxt = path_groundtruths.at(j).string();

        // Create our trajectory object
        ResultTrajectory traj(path_esttxt, path_gttxt, "posyaw");

        // Calculate ATE error for this dataset
        Statistics error_ori, error_pos;
        // cout << "est_poses_aignedtoGT size before ATE" << est_poses_aignedtoGT.size() << endl;
        // cout << " gt_poses size before ATE" << gt_poses.size() << endl;
        traj.calculate_ate(error_ori, error_pos);
        // cout << "est_poses_aignedtoGT size after ATE" << est_poses_aignedtoGT.size() << endl;
        // cout << " gt_poses size after ATE" << gt_poses.size() << endl;
        ate_dataset_ori.values.push_back(error_ori.rmse);
        ate_dataset_pos.values.push_back(error_pos.rmse);

        // Calculate ATE 2D error for this dataset
        // Statistics error_ori_2d, error_pos_2d;
        // traj.calculate_ate_2d(error_ori_2d, error_pos_2d);
        // ate_2d_dataset_ori.values.push_back(error_ori_2d.rmse);
        // ate_2d_dataset_pos.values.push_back(error_pos_2d.rmse);

        // Calculate RPE error for this dataset
        std::map<double, std::pair<Statistics, Statistics>> error_rpe;
        // cout << "est_poses_aignedtoGT size before RPE" << est_poses_aignedtoGT.size() << endl;
        // cout << " gt_poses size before RPE" << gt_poses.size() << endl;
        traj.calculate_rpe(segments, error_rpe);
        // cout << "est_poses_aignedtoGT size after RPE" << est_poses_aignedtoGT.size() << endl;
        // cout << " gt_poses size after RPE" << gt_poses.size() << endl;
        for (const auto &elm : error_rpe) {
          rpe_dataset.at(elm.first).first.values.insert(rpe_dataset.at(elm.first).first.values.end(), elm.second.first.values.begin(),
                                                        elm.second.first.values.end());
          rpe_dataset.at(elm.first).first.timestamps.insert(rpe_dataset.at(elm.first).first.timestamps.end(),
                                                            elm.second.first.timestamps.begin(), elm.second.first.timestamps.end());
          rpe_dataset.at(elm.first).second.values.insert(rpe_dataset.at(elm.first).second.values.end(), elm.second.second.values.begin(),
                                                         elm.second.second.values.end());
          rpe_dataset.at(elm.first).second.timestamps.insert(rpe_dataset.at(elm.first).second.timestamps.end(),
                                                             elm.second.second.timestamps.begin(), elm.second.second.timestamps.end());
        }
      }
      // Compute our mean ATE score
      ate_dataset_ori.calculate();
      ate_dataset_pos.calculate();
    //   ate_2d_dataset_ori.calculate();
    //   ate_2d_dataset_pos.calculate();

      // Print stats for this specific dataset
      std::string prefix = (ate_dataset_ori.mean > 10 || ate_dataset_pos.mean > 10) ? "\033[31m" :  "";
      printf("%s\tATE: mean_ori = %.3f | mean_pos = %.3f (%d runs)\n" , prefix.c_str(), ate_dataset_ori.mean,
                  ate_dataset_pos.mean, (int)ate_dataset_pos.values.size());
      printf("\tATE: std_ori  = %.3f | std_pos  = %.3f\n", ate_dataset_ori.std, ate_dataset_pos.std);
    //   printf("\tATE 2D: mean_ori = %.3f | mean_pos = %.3f (%d runs)\n", ate_2d_dataset_ori.mean, ate_2d_dataset_pos.mean,
    //               (int)ate_2d_dataset_ori.values.size());
    //   printf("\tATE 2D: std_ori  = %.5f | std_pos  = %.5f\n", ate_2d_dataset_ori.std, ate_2d_dataset_pos.std);
      for (auto &seg : rpe_dataset) {
        seg.second.first.calculate();
        seg.second.second.calculate();
        // PRINT_DEBUG("\tRPE: seg %d - mean_ori = %.3f | mean_pos = %.3f (%d
        // samples)\n",(int)seg.first,seg.second.first.mean,seg.second.second.mean,(int)seg.second.second.values.size());
        printf("\tRPE: seg %d - median_ori = %.4f | median_pos = %.4f (%d samples)\n", (int)seg.first, seg.second.first.median,
                    seg.second.second.median, (int)seg.second.second.values.size());
        // PRINT_DEBUG("RPE: seg %d - std_ori  = %.3f | std_pos  = %.3f\n",(int)seg.first,seg.second.first.std,seg.second.second.std);
      }

      // Update the global ATE error stats
      std::string algo = path_algorithms.at(i).filename().string();
      algo_ate.at(algo).at(j).first = ate_dataset_ori;
      algo_ate.at(algo).at(j).second = ate_dataset_pos;
    //   algo_ate_2d.at(algo).at(j).first = ate_2d_dataset_ori;
    //   algo_ate_2d.at(algo).at(j).second = ate_2d_dataset_pos;

      // Update the global RPE error stats
      for (const auto &elm : rpe_dataset) {
        algo_rpe.at(algo).at(elm.first).first.values.insert(algo_rpe.at(algo).at(elm.first).first.values.end(),
                                                            elm.second.first.values.begin(), elm.second.first.values.end());
        algo_rpe.at(algo).at(elm.first).first.timestamps.insert(algo_rpe.at(algo).at(elm.first).first.timestamps.end(),
                                                                elm.second.first.timestamps.begin(), elm.second.first.timestamps.end());
        algo_rpe.at(algo).at(elm.first).second.values.insert(algo_rpe.at(algo).at(elm.first).second.values.end(),
                                                             elm.second.second.values.begin(), elm.second.second.values.end());
        algo_rpe.at(algo).at(elm.first).second.timestamps.insert(algo_rpe.at(algo).at(elm.first).second.timestamps.end(),
                                                                 elm.second.second.timestamps.begin(), elm.second.second.timestamps.end());
      }
    }
  }

  //===============================================================================
  //===============================================================================
  //===============================================================================

  // Finally print the ATE for all the runs
  printf("============================================\n");
  printf("ATE LATEX TABLE\n");
  printf("============================================\n");
  for (size_t i = 0; i < path_groundtruths.size(); i++) {
    std::string gtname = path_groundtruths.at(i).stem().string();
    boost::replace_all(gtname, "_", "\\_");
    printf(" & \\textbf{%s}", gtname.c_str());
  }
  printf(" & \\textbf{Average} \\\\\\hline\n");
  for (auto &algo : algo_ate) {
    std::string algoname = algo.first;
    boost::replace_all(algoname, "_", "\\_");
    printf(algoname.c_str());
    double sum_ori = 0.0;
    double sum_pos = 0.0;
    int sum_ct = 0;
    for (auto &seg : algo.second) {
      if (seg.first.values.empty() || seg.second.values.empty()) {
        printf(" & - / -");
      } else {
        seg.first.calculate();
        seg.second.calculate();
        printf(" & %.3f / %.3f", seg.first.mean, seg.second.mean);
        sum_ori += seg.first.mean;
        sum_pos += seg.second.mean;
        sum_ct++;
      }
    }
    printf(" & %.3f / %.3f \\\\\n", sum_ori / sum_ct, sum_pos / sum_ct);
  }
  printf("============================================\n");

  // Finally print the ATE_2D for all the runs
//   PRINT_INFO("============================================\n");
//   PRINT_INFO("ATE 2D LATEX TABLE\n");
//   PRINT_INFO("============================================\n");
//   for (size_t i = 0; i < path_groundtruths.size(); i++) {
//     std::string gtname = path_groundtruths.at(i).stem().string();
//     boost::replace_all(gtname, "_", "\\_");
//     PRINT_INFO(" & \\textbf{%s}", gtname.c_str());
//   }
//   PRINT_INFO(" & \\textbf{Average} \\\\\\hline\n");
//   for (auto &algo : algo_ate_2d) {
//     std::string algoname = algo.first;
//     boost::replace_all(algoname, "_", "\\_");
//     PRINT_INFO(algoname.c_str());
//     double sum_ori = 0.0;
//     double sum_pos = 0.0;
//     int sum_ct = 0;
//     for (auto &seg : algo.second) {
//       if (seg.first.values.empty() || seg.second.values.empty()) {
//         PRINT_INFO(" & - / -");
//       } else {
//         seg.first.calculate();
//         seg.second.calculate();
//         PRINT_INFO(" & %.3f / %.3f", seg.first.mean, seg.second.mean);
//         sum_ori += seg.first.mean;
//         sum_pos += seg.second.mean;
//         sum_ct++;
//       }
//     }
//     PRINT_INFO(" & %.3f / %.3f \\\\\n", sum_ori / sum_ct, sum_pos / sum_ct);
//   }
//   PRINT_INFO("============================================\n");

  // Finally print the RPE for all the runs
  printf("============================================\n");
  printf("RPE LATEX TABLE\n");
  printf("============================================\n");
  for (const auto &len : segments) {
    printf(" & \\textbf{%dm}", (int)len);
  }
  printf(" \\\\\\hline\n");
  for (auto &algo : algo_rpe) {
    std::string algoname = algo.first;
    boost::replace_all(algoname, "_", "\\_");
    printf(algoname.c_str());
    for (auto &seg : algo.second) {
      seg.second.first.calculate();
      seg.second.second.calculate();
      printf(" & %.3f / %.3f", seg.second.first.mean, seg.second.second.mean);
    }
    printf(" \\\\\n");
  }
  printf("============================================\n");

  // Plot line colors
  std::vector<std::string> colors = {"blue", "red", "black", "green", "cyan", "magenta"};
  std::vector<std::string> linestyle = {"-", "--", "-."};
  assert(algo_rpe.size() <= colors.size() * linestyle.size());

  // Parameters
  std::map<std::string, std::string> params_rpe;
  params_rpe.insert({"notch", "false"});
  params_rpe.insert({"sym", ""});

  // Plot this figure
  matplotlibcpp::figure_size(1000, 700);
  matplotlibcpp::subplot(2, 1, 1);
  matplotlibcpp::title("Relative Orientation Error");

  // Plot each RPE next to each other
  double width = 1.0 / (algo_rpe.size() + 1);
  double spacing = width / (algo_rpe.size() + 1);
  std::vector<double> xticks;
  std::vector<std::string> labels;
  int ct_algo = 0;
  double ct_pos = 0;
  for (auto &algo : algo_rpe) {
    // Start based on what algorithm we are doing
    xticks.clear();
    labels.clear();
    ct_pos = 1 + ct_algo * (width + spacing);
    // Loop through each length type
    for (auto &seg : algo.second) {
      xticks.push_back(ct_pos - (algo_rpe.size() * (width + spacing) - width) / 2);
      labels.push_back(std::to_string((int)seg.first));
      matplotlibcpp::boxplot(seg.second.first.values, ct_pos, width, colors.at(ct_algo % colors.size()),
                             linestyle.at(ct_algo / colors.size()), params_rpe);
      ct_pos += 1 + 3 * width;
    }
    // Move forward
    ct_algo++;
  }
  // Add "fake" plots for our legend
  ct_algo = 0;
  for (const auto &algo : algo_rpe) {
    std::map<std::string, std::string> params_empty;
    params_empty.insert({"label", algo.first});
    params_empty.insert({"linestyle", linestyle.at(ct_algo / colors.size())});
    params_empty.insert({"color", colors.at(ct_algo % colors.size())});
    std::vector<double> vec_empty;
    matplotlibcpp::plot(vec_empty, vec_empty, params_empty);
    ct_algo++;
  }
  // Plot each RPE next to each other
  matplotlibcpp::xlim(0.5, ct_pos - 0.5 - 3 * width);
  matplotlibcpp::xticks(xticks, labels);
  matplotlibcpp::ylabel("orientation error (deg)");
  matplotlibcpp::legend();
  matplotlibcpp::subplot(2, 1, 2);
  ct_algo = 0;
  ct_pos = 0;
  for (auto &algo : algo_rpe) {
    // Start based on what algorithm we are doing
    ct_pos = 1 + ct_algo * (width + spacing);
    // Loop through each length type
    for (auto &seg : algo.second) {
      matplotlibcpp::boxplot(seg.second.second.values, ct_pos, width, colors.at(ct_algo % colors.size()),
                             linestyle.at(ct_algo / colors.size()), params_rpe);
      ct_pos += 1 + 3 * width;
    }
    // Move forward
    ct_algo++;
  }
  // Add "fake" plots for our legend
  ct_algo = 0;
  for (const auto &algo : algo_rpe) {
    std::map<std::string, std::string> params_empty;
    params_empty.insert({"label", algo.first});
    params_empty.insert({"linestyle", linestyle.at(ct_algo / colors.size())});
    params_empty.insert({"color", colors.at(ct_algo % colors.size())});
    std::vector<double> vec_empty;
    matplotlibcpp::plot(vec_empty, vec_empty, params_empty);
    ct_algo++;
  }
  // Display to the user
  matplotlibcpp::xlim(0.5, ct_pos - 0.5 - 3 * width);
  matplotlibcpp::xticks(xticks, labels);
  matplotlibcpp::title("Relative Position Error");
  matplotlibcpp::ylabel("translational error (m)");
  matplotlibcpp::xlabel("sub-segment lengths (m)");
  matplotlibcpp::tight_layout();
  matplotlibcpp::show(true);
  // Done!
  return EXIT_SUCCESS;
}