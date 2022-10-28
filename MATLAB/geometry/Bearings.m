classdef Bearings < handle & matlab.mixin.Copyable
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
    bearings;
    depths;
  end
  
  methods
    function m = Bearings( varargin )
      
      if nargin == 1
          % time_num
          m.time_num = varargin{1};
          m.bearings = repmat(Point(),1,m.time_num);
          m.depths = zeros(1,m.time_num);
      elseif nargin == 3
          noise = varargin{3};
          %traj1 traj2
          traj1 = varargin{1};
          traj2 = varargin{2};
          m.time_num = traj1.time_num;

          m.bearings = repmat(Point(),1,m.time_num);
          m.depths = zeros(1,m.time_num);
          
          for i = 1:m.time_num
            pose1 = traj1.T(i);
            pose2 = traj2.T(i);
            rand_value = 1 + (rand(1,1) - 2) * noise;
            bias_world = pose2.R * traj2.bias + pose2.t;
            measurement = pose1.R'* (bias_world * rand_value - pose1.t);
            d = norm(measurement);
            m.depths(i) = d;
            m.bearings(i) = measurement / d;
          end
      end
    end
    
    function exchange(this, other_bearing)
        assert(this.time_num == other_bearing.time_num)
        num = floor(0.1*this.time_num);
        for i = 1:num
            index = randi(this.time_num);
            middile = other_bearing.bearings(index);
            other_bearing.bearings(index) = this.bearings(index);
            this.bearings(index) = middile;
        end
    end
       
    % Additional interface properties
    function B = bearing(this, i)
      assert((i <= this.time_num) && (i >= 1));
      B = this.bearings(i);
    end
    
    function D = depth(this, i)
      assert(i <= this.time_num && i >= 1);
      D = this.depths(i);
    end
  end
   
end
