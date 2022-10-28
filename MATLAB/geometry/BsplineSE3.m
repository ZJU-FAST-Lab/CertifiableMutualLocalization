classdef BsplineSE3 < handle & matlab.mixin.Copyable
  %Pose Object with SE(3) data fields
  %   A point in the manifold SE(3) is created from the following data:
  %   - translation t and rotation R
  %
  %   Technical:
  %   This is a *handle* *copyable* object, that means that the object
  %   created by the constructor is a handle (pointer) to the data.
  %   If you want to get an *independent* copy of the object (clone)
  %   use the method *copy* in the object.
  %   Intuitively, the pose of an object is inherent to this
  %   and should be assigned to a different object, but another object
  %   could have the same value (copy).
  %
  %   See also handle, matlab.mixin.Copyable.
  
  properties
    dt
    timestamp_start
    timestamp_end
    control_points
    control_points_mat
    trajectory_points
    sim_freq_imu
    gravity_mag
    sigma_w
    sigma_a
    sigma_wb
    sigma_ab
    true_bias_gyro
    true_bias_accel
    omega_gt
    acceleration_gt
    omega
    acceleration
  end
  
  methods
      function traj = BsplineSE3( traj_points )
          % Find the average frequency to use as our uniform timesteps
          sumdt = 0;
          for i = 1:size(traj_points,2)-1
            sumdt = sumdt + traj_points(1,i+1) - traj_points(1,i);
          end
          dt = sumdt / (size(traj_points,2)-1);
          traj.dt = max(0.05,dt);
          
          traj.trajectory_points = containers.Map('KeyType','double','ValueType','any');
          traj.control_points = containers.Map('KeyType','double','ValueType','any');

          for i = 1:size(traj_points,2)
            traj.trajectory_points(traj_points(1,i)) = [bot_quat_to_matrix(traj_points(5:8,i)), traj_points(2:4,i); zeros(1,3), 1];
          end
          % then create spline control points
          times = keys(traj.trajectory_points);
          timestamp_curr = times{1,1};
          traj.timestamp_start = timestamp_curr + 2 * dt;
          traj.timestamp_end = times{1,end} - 2 * dt;
          while (true)
        
            %Get bounding posed for the current time
            [success,t0, pose0, t1, pose1] = traj.find_bounding_poses(timestamp_curr, traj.trajectory_points);
        
            %If we didn't find a bounding pose, then that means we are at the end of the dataset
            %Thus break out of this loop since we have created our max number of control points
            if (~success)
              break
            end
        
            %Linear interpolation and append to our control points
            lambda = (timestamp_curr - t0) / (t1 - t0)
            pose_interp = exp_se3(lambda * log_se3(pose1 * inv_se3(pose0))) * pose0;
            pose_interp - pose1
            traj.control_points(timestamp_curr) = pose_interp;
            traj.control_points_mat = [traj.control_points_mat pose_interp(1:3,4)];
            timestamp_curr = timestamp_curr + dt;
          end
          traj.sim_freq_imu = 100;
          traj.gravity_mag = 9.81;
          traj.sigma_w = 1.6968e-04;
          traj.sigma_a = 2.0000e-3;
          traj.sigma_wb = 1.9393e-05;
          traj.sigma_ab = 3.0000e-03;
          traj.true_bias_gyro = zeros(3,1);
          traj.true_bias_accel = zeros(3,1);

          traj.generate_imu_data();
      end

      function generate_imu_data(this)
          time_stamp = this.timestamp_start;
          while(true)
              
              [success, R_GtoI, p_IinG, w_IinI, v_IinG, alpha_IinI, a_IinG] = this.get_acceleration(time_stamp);
              if(~success)
                  break
              end

              % Transform omega and linear acceleration into imu frame
              omega_inI = w_IinI;
              gravity = [0;0;this.gravity_mag];
              accel_inI = R_GtoI * (a_IinG + gravity);
            
              % Now add noise to these measurements
              delta_t = 1.0 / this.sim_freq_imu;
              wm = omega_inI + this.true_bias_gyro + this.sigma_w / sqrt(delta_t) * rand(3,1);
              am = accel_inI + this.true_bias_accel + this.sigma_a / sqrt(delta_t) * rand(3,1);
            
              % Move the biases forward in time
              this.true_bias_gyro = this.true_bias_gyro + this.sigma_wb * sqrt(delta_t) * rand(3,1);
              this.true_bias_accel = this.true_bias_accel + this.sigma_ab * sqrt(delta_t) * rand(3,1);
    
              this.omega_gt = [this.omega_gt omega_inI];
              this.acceleration_gt = [this.acceleration_gt accel_inI];
              this.omega = [this.omega wm];
              this.acceleration = [this.acceleration am];

              time_stamp = time_stamp + delta_t;
          end
      end


      function [success, R_GtoI, p_IinG, w_IinI, v_IinG, alpha_IinI, a_IinG] = get_acceleration(this, time_stamp)
          R_GtoI=0; p_IinG=0; w_IinI=0; v_IinG=0; alpha_IinI=0; a_IinG=0;
          [success,t0, pose0, t1, pose1, t2, pose2, t3, pose3] = this.find_bounding_control_points(time_stamp, this.control_points);
          if(success)
              %Our De Boor-Cox matrix scalars
              DT = (t2 - t1);
              u = (time_stamp - t1) / DT;
              b0 = 1.0 / 6.0 * (5 + 3 * u - 3 * u * u + u * u * u);
              b1 = 1.0 / 6.0 * (1 + 3 * u + 3 * u * u - 2 * u * u * u);
              b2 = 1.0 / 6.0 * (u * u * u);
              b0dot = 1.0 / (6.0 * DT) * (3 - 6 * u + 3 * u * u);
              b1dot = 1.0 / (6.0 * DT) * (3 + 6 * u - 6 * u * u);
              b2dot = 1.0 / (6.0 * DT) * (3 * u * u);
              b0dotdot = 1.0 / (6.0 * DT * DT) * (-6 + 6 * u);
              b1dotdot = 1.0 / (6.0 * DT * DT) * (6 - 12 * u);
              b2dotdot = 1.0 / (6.0 * DT * DT) * (6 * u);
            
              %Cache some values we use alot
              omega_10 = log_se3(inv_se3(pose0) * pose1);
              omega_21 = log_se3(inv_se3(pose1) * pose2);
              omega_32 = log_se3(inv_se3(pose2) * pose3);
              omega_10_hat = hat_se3(omega_10);
              omega_21_hat = hat_se3(omega_21);
              omega_32_hat = hat_se3(omega_32);
            
              %Calculate interpolated poses
              A0 = exp_se3(b0 * omega_10);
              A1 = exp_se3(b1 * omega_21);
              A2 = exp_se3(b2 * omega_32);
              A0dot = b0dot * omega_10_hat * A0;
              A1dot = b1dot * omega_21_hat * A1;
              A2dot = b2dot * omega_32_hat * A2;
              A0dotdot = b0dot * omega_10_hat * A0dot + b0dotdot * omega_10_hat * A0;
              A1dotdot = b1dot * omega_21_hat * A1dot + b1dotdot * omega_21_hat * A1;
              A2dotdot = b2dot * omega_32_hat * A2dot + b2dotdot * omega_32_hat * A2;
              
              % Get the interpolated pose
              pose_interp = pose0 * A0 * A1 * A2;
              R_GtoI = pose_interp(1:3,1:3)';
              p_IinG = pose_interp(1:3,4);
            
              % Get the interpolated velocities
              % NOTE: Rdot = R*skew(omega) => R^T*Rdot = skew(omega)
              vel_interp = pose0 * (A0dot * A1 * A2 + A0 * A1dot * A2 + A0 * A1 * A2dot);
              w_IinI = vee(pose_interp(1:3,1:3)' * vel_interp(1:3,1:3));
              v_IinG = vel_interp(1:3,4);
            
              % Finally get the interpolated velocities
              % NOTE: Rdot = R*skew(omega)
              % NOTE: Rdotdot = Rdot*skew(omega) + R*skew(alpha) => R^T*(Rdotdot-Rdot*skew(omega))=skew(alpha)
              acc_interp = pose0 * (A0dotdot * A1 * A2 + A0 * A1dotdot * A2 + A0 * A1 * A2dotdot + ...
                  2 * A0dot * A1dot * A2 + 2 * A0 * A1dot * A2dot + 2 * A0dot * A1 * A2dot);
              omegaskew = pose_interp(1:3,1:3)' * vel_interp(1:3,1:3);
              alpha_IinI = vee(pose_interp(1:3,1:3)' * (acc_interp(1:3,1:3) - vel_interp(1:3,1:3) * omegaskew));
              a_IinG = acc_interp(1:3,4);
          end
      end

      function [success,t0, pose0, t1, pose1] = find_bounding_poses(this, timestamp_curr, trajectory_points)
          t0 = -1;
          i0 = -1;
          t1 = -1;
          i1 = -1;
          pose0 = zeros(4,4);
          pose1 = zeros(4,4);
          times = keys(trajectory_points);
          
          for i = 1:size(times,2)-1
            t_i = times{1,i};
            t_i_1 = times{1,i+1};
            if(t_i== timestamp_curr)
                t0 = t_i;
                i0 = i;
                pose0 = trajectory_points(t0);
                break
            end

            if(t_i < timestamp_curr && t_i_1 > timestamp_curr)
                t0 = t_i;
                i0 = i;
                pose0 = trajectory_points(t0);
                break
            end
          end

          for i = 1:size(times,2)
            t_i = times{1,i};
            if(t_i> timestamp_curr)
                t1 = t_i;
                i1 = i;
                pose1 = trajectory_points(t1);
                break
            end
          end
          success = (t0~=-1 && t1~=-1);
          if(success)
             assert(t0 < t1)
          end
      end

      function [success,t0, pose0, t1, pose1, t2, pose2, t3, pose3] = find_bounding_control_points(this, timestamp_curr, control_points)
          [success,t1, pose1, t2, pose2] = this.find_bounding_poses(timestamp_curr, control_points);
          times = cell2mat(keys(control_points));
          if(success)
             [~,iter_t1] = find(times==t1);
             [~,iter_t2] = find(times==t2);
             iter_t0 = iter_t1 - 1;
             iter_t3 = iter_t2 + 1;
             if(iter_t0 >= 1 && iter_t3 <= size(times,2))
                 t0 = times(iter_t0);
                 pose0 = control_points(t0);
                 t3 = times(iter_t3);
                 pose3 = control_points(t3);
             else
                 success = false;
             end
          end
          
          if(~success)
             t0 = -1;
             t3 = -1;
             pose0 = zeros(4,4);
             pose3 = zeros(4,4);
          end
      end
    
  end
   
end
