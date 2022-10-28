classdef Trajectory < handle & matlab.mixin.Copyable
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
    time_num;
    poses;
    bias;
  end
  
  methods
    function traj = Trajectory( varargin )
      if nargin == 1
          traj.time_num = varargin{1};
          traj.poses = repmat(Pose(),1,traj.time_num);
      end
    end
    
    function time_num = get.time_num(this)
      time_num = this.time_num;
    end
    function traj = get.poses(this)
      traj = this.poses;
    end
    
    function setPose(this, i, T)
      assert(i <= this.time_num && i >= 1);
      this.poses(i) = T;
    end
       
    % Additional interface properties
    function t = T(this, i)
      t = this.poses(i);
    end
    
    function line(this,data)
      for i = 1:this.time_num
          t = (i-1)*rand(1,1);
          T = [cos(t) -sin(t) 0 data.x + rand(1,1);
              sin(t) cos(t) 0 i-1 + rand(1,1);
              0 0 1 0;
              0 0 0 1];
          P = Pose(T);
          this.poses(i) =  P;
      end
      this.bias = zeros(3,1);
    end
    
    function random(this, param)
      for i = 1:this.time_num
          P = Pose.rand(param.only_yaw);
          this.poses(i) =  P;
      end
      if(param.if_bias)
          this.bias = randn(3,1)*5;
      else
          this.bias = zeros(3,1);
      end
    end

    function random_rotation(this)
      for i = 1:this.time_num
        this.poses(i).T(1:3,1:3) = rnd.rot();
      end
    end

    function random_translation(this, scale)
      for i = 1:this.time_num
        this.poses(i).T(1:3,4) = rand(3,1) * scale;
      end
    end
    
    function local_traj = scale_traj(traj, scale)
        local_traj = Trajectory(traj.time_num);
        local_traj.bias = traj.bias;
        pose_init = traj.T(1);
        pose_init_inv = Pose.inv(pose_init);
        for i = 1:traj.time_num
            local_pose = pose_init_inv.mtimes(traj.T(i));
            scaled_local_pose = Pose(scale * local_pose.t, local_pose.R);
            local_traj.setPose(i,scaled_local_pose);
        end
    end
  end
end
