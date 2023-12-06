#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

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

Eigen::Matrix<double, 3, 3> skew_x(const Eigen::Matrix<double, 3, 1> &w) {
  Eigen::Matrix<double, 3, 3> w_x;
  w_x << 0, -w(2), w(1), w(2), 0, -w(0), -w(1), w(0), 0;
  return w_x;
}

Eigen::Matrix<double, 3, 3> quat_2_Rot(const Eigen::Matrix<double, 4, 1> &q) {
  Eigen::Matrix<double, 3, 3> q_x = skew_x(q.block(0, 0, 3, 1));
  Eigen::MatrixXd Rot = (2 * std::pow(q(3, 0), 2) - 1) * Eigen::MatrixXd::Identity(3, 3) - 2 * q(3, 0) * q_x +
                        2 * q.block(0, 0, 3, 1) * (q.block(0, 0, 3, 1).transpose());
  return Rot;
}

Eigen::Matrix<double, 3, 3> rot_z(double t) {
  Eigen::Matrix<double, 3, 3> r;
  double ct = cos(t);
  double st = sin(t);
  r << ct, -st, 0.0, st, ct, 0.0, 0.0, 0.0, 1.0;
  return r;
}

void load_data_R(std::string path_traj, std::vector<double> &times, std::vector<Eigen::Matrix<double, 3, 1>> &poses,
                 std::vector<Eigen::Matrix<double, 3, 1>> &euler )
{ 

  // Try to open our trajectory file
  std::ifstream file(path_traj);
  if (!file.is_open()) {
    printf("[LOAD]: Unable to open trajectory file...\n" );
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
    Eigen::Matrix<double, 20, 1> data;

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
      poses.push_back(data.block(1, 0, 3, 1));
      euler.push_back(data.block(4, 0, 3, 1));

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


Eigen::Matrix<double, 4, 1> euler_to_q(Eigen::Matrix<double, 3, 1> data){
  Eigen::Matrix<double, 4, 1> q;
  double roll = data[0];
  double pitch = data[0];
  double yaw = data[0];

  q[0] = sin(roll/2) * cos(pitch/2) * cos(yaw/2) - cos(roll/2) * sin(pitch/2) * sin(yaw/2);
  q[1] = cos(roll/2) * sin(pitch/2) * cos(yaw/2) + sin(roll/2) * cos(pitch/2) * sin(yaw/2);
  q[2] = cos(roll/2) * cos(pitch/2) * sin(yaw/2) - sin(roll/2) * sin(pitch/2) * cos(yaw/2);
  q[3] = cos(roll/2) * cos(pitch/2) * cos(yaw/2) + sin(roll/2) * sin(pitch/2) * sin(yaw/2);

  return q;
}



int main(){
    // std::ifstream file_in("/home/zhangyanyu/gt.txt");
    // std::ofstream file_out("/home/zhangyanyu/gt_converted.txt");
    std::vector<double> gt_times;
    std::vector<Eigen::Matrix<double, 3, 1>>  gt_pose;
    // std::vector<Eigen::Matrix<double, 4, 1>>  gt_q;
    // std::vector<Eigen::Matrix<double, 3, 1>>  gt_R;

    std::vector<Eigen::Matrix<double, 3, 1>> euler;

    std::string path_in = "/home/zhangyanyu/open_vins/fast_lio.txt";
    std::string path_out = "/home/zhangyanyu/fast_lio_gt.txt";
    // load_data(path_in,gt_times, gt_pose, gt_q);
    load_data_R(path_in,gt_times, gt_pose, euler);


    // Eigen::Matrix<double, 3, 3> rotation_matrix = quat_2_Rot(gt_q[0]);
    // Eigen::Matrix<double, 3, 1> first_translation = gt_pose[0];

    std::vector<Eigen::Matrix<double, 3, 1>> converted_poses;
    std::vector<Eigen::Matrix<double, 4, 1>> converted_quaternions;

    // for (size_t i = 0; i < gt_times.size(); i++) {
    //     Eigen::Matrix<double, 3, 1> pose =  gt_pose[i] -  first_translation;
    //     Eigen::Matrix<double, 4, 1> quaternion = rot_2_quat( quat_2_Rot(gt_q[i]) * rotation_matrix.transpose());
    //     converted_poses.push_back(pose);
    //     converted_quaternions.push_back(quaternion);
    // }



    for (size_t i = 0; i < gt_times.size(); i++) {
        Eigen::Matrix<double, 3, 1> pose =  gt_pose[i];

        // Eigen::AngleAxisd rollAngle(roll[i](0,0), Eigen::Vector3d::UnitZ());
        // Eigen::AngleAxisd yawAngle(yaw[i](0,0), Eigen::Vector3d::UnitY());
        // Eigen::AngleAxisd pitchAngle(pitch[i](0,0), Eigen::Vector3d::UnitX());
        // Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
        // Eigen::Matrix3d rotationMatrix = q.matrix();

        // Eigen::Matrix<double, 3, 3> R = rot_z(yaw[i](0,0));
        // Eigen::Matrix<double, 4, 1> quaternion = rot_2_quat(R);
        Eigen::Matrix<double, 3, 1> euler_angle = euler[i];
        Eigen::Matrix<double, 4, 1> quaternion = euler_to_q(euler_angle);
        converted_poses.push_back(pose);
        converted_quaternions.push_back(quaternion);
    }
    
    std::ofstream file_out(path_out);

    // if (file_out.is_open()) {
    //     for (size_t i = 0; i < gt_times.size(); i++) {
    //         file_out << setprecision(12) << gt_times[i] << setprecision(9) << " " << converted_poses[i](0) << " " << converted_poses[i](1) << " " << converted_poses[i](2) << " "
    //                  << converted_quaternions[i](0) << " " << converted_quaternions[i](1) << " " << converted_quaternions[i](2) << " " << converted_quaternions[i](3) << "\n";
    //     }
    //     file_out.close();
    // } else {
    //     std::cerr << "Error opening the output file." << std::endl;
    //     return 1;
    // }


    //convert fast_lio R_2_q
        if (file_out.is_open()) {
        for (size_t i = 0; i < gt_times.size(); i++) {
            file_out << setprecision(12) << gt_times[i] << setprecision(12) << " " << converted_poses[i](0) << " " << converted_poses[i](1) << " " << converted_poses[i](2) << " "
                     << converted_quaternions[i](0) << " " << converted_quaternions[i](1) << " " << converted_quaternions[i](2) << " " << converted_quaternions[i](3) << "\n";
        }
        file_out.close();
    } else {
        std::cerr << "Error opening the output file." << std::endl;
        return 1;
    }

    return 0;


}