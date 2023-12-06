import numpy as np

def rot_2_quat(rot):
    q = np.zeros((4, 1), dtype=float)
    T = np.trace(rot)
    
    if (rot[0, 0] >= T) and (rot[0, 0] >= rot[1, 1]) and (rot[0, 0] >= rot[2, 2]):
        q[0] = np.sqrt((1 + 2 * rot[0, 0] - T) / 4)
        q[1] = (1 / (4 * q[0])) * (rot[0, 1] + rot[1, 0])
        q[2] = (1 / (4 * q[0])) * (rot[0, 2] + rot[2, 0])
        q[3] = (1 / (4 * q[0])) * (rot[1, 2] - rot[2, 1])
    elif (rot[1, 1] >= T) and (rot[1, 1] >= rot[0, 0]) and (rot[1, 1] >= rot[2, 2]):
        q[1] = np.sqrt((1 + 2 * rot[1, 1] - T) / 4)
        q[0] = (1 / (4 * q[1])) * (rot[0, 1] + rot[1, 0])
        q[2] = (1 / (4 * q[1])) * (rot[1, 2] + rot[2, 1])
        q[3] = (1 / (4 * q[1])) * (rot[2, 0] - rot[0, 2])
    elif (rot[2, 2] >= T) and (rot[2, 2] >= rot[0, 0]) and (rot[2, 2] >= rot[1, 1]):
        q[2] = np.sqrt((1 + 2 * rot[2, 2] - T) / 4)
        q[0] = (1 / (4 * q[2])) * (rot[0, 2] + rot[2, 0])
        q[1] = (1 / (4 * q[2])) * (rot[1, 2] + rot[2, 1])
        q[3] = (1 / (4 * q[2])) * (rot[0, 1] - rot[1, 0])
    else:
        q[3] = np.sqrt((1 + T) / 4)
        q[0] = (1 / (4 * q[3])) * (rot[1, 2] - rot[2, 1])
        q[1] = (1 / (4 * q[3])) * (rot[2, 0] - rot[0, 2])
        q[2] = (1 / (4 * q[3])) * (rot[0, 1] - rot[1, 0])
    
    if q[3] < 0:
        q = -q
    
    # Normalize and return
    q = q / np.linalg.norm(q)
    return q

def skew_x(w):
    w_x = np.array([[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]])
    return w_x

def quat_2_Rot(q):
    q_x = skew_x(q[:3])
    Rot = (2 * q[3]**2 - 1) * np.identity(3) - 2 * q[3] * q_x + 2 * q[:3].dot(q[:3].T)
    return Rot

def rot_z(t):
    ct = np.cos(t)
    st = np.sin(t)
    r = np.array([[ct, -st, 0.0], [st, ct, 0.0], [0.0, 0.0, 1.0]])
    return r

def load_data_R(path_traj, times, poses, euler):
    try:
        with open(path_traj, 'r') as file:
            for line in file:
                if line.strip().startswith('#'):
                    continue
                fields = line.strip().split()
                data = np.array([float(field) for field in fields])
                times.append(data[0])
                poses.append(data[1:4])
                euler.append(data[4:7])
    
    except IOError:
        print("[LOAD]: Unable to open trajectory file...")
        exit(1)
    
    if not times:
        print("[LOAD]: Could not parse any data from the file!!")
        print(f"[LOAD]: {path_traj}")
        exit(1)
    
    if len(times) != len(poses):
        print("[LOAD]: Parsing error, pose and timestamps do not match!!")
        print(f"[LOAD]: {path_traj}")
        exit(1)

def euler_to_q(data):
    q = np.zeros((4, 1), dtype=float)
    roll = data[0]
    pitch = data[1]
    yaw = data[2]

    q[0] = np.sin(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) - np.cos(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)
    q[1] = np.cos(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2)
    q[2] = np.cos(roll / 2) * np.cos(pitch / 2) * np.sin(yaw / 2) - np.sin(roll / 2) * np.sin(pitch / 2) * np.cos(yaw / 2)
    q[3] = np.cos(roll / 2) * np.cos(pitch / 2) * np.cos(yaw / 2) + np.sin(roll / 2) * np.sin(pitch / 2) * np.sin(yaw / 2)

    return q

def main():
    gt_times = []
    gt_pose = []
    euler = []
    path_in = "/home/zhangyanyu/open_vins/fast.txt"
    path_out = "/home/zhangyanyu/fast_lio_gt.txt"
    load_data_R(path_in, gt_times, gt_pose, euler)

    converted_poses = []
    converted_quaternions = []

    for i in range(len(gt_times)):
        pose = gt_pose[i]

        euler_angle = euler[i]
        quaternion = euler_to_q(euler_angle)
        converted_poses.append(pose)
        converted_quaternions.append(quaternion)

    with open(path_out, 'w') as file_out:
        if file_out:
            for i in range(len(gt_times)):
                file_out.write(f"{gt_times[i]:.12f} {converted_poses[i][0]:.9f} {converted_poses[i][1]:.9f} {converted_poses[i][2]:.9f} "
                            f"{converted_quaternions[i][0, 0]:.9f} {converted_quaternions[i][1, 0]:.9f} {converted_quaternions[i][2, 0]:.9f} {converted_quaternions[i][3, 0]:.9f}\n")
        else:
            print("Error opening the output file.")

if __name__ == "__main__":
    main()
