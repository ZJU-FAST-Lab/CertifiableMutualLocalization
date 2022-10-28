
    base_X = R_opt_iteration' * R_opt_iteration;
    resolution = 100;
    Z_value = zeros(resolution, resolution);
    index_robot = randi(numRobots);
    rand_index = randi(3*numRobots, 2,2);
    
    for i = 1:resolution
        roll = abs(i - resolution / 2) * 2 * pi / 10;
        for j = 1:resolution
            pitch = abs(j - resolution / 2)  * 2 * pi / 10;

            rpy = [roll; pitch; 0];
            perturbance = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat(rpy));

            if(i > resolution / 2)
               index_i = rand_index(1,1); 
            else
               index_i = rand_index(1,2);
            end

            if(j > resolution / 2)
               index_j = rand_index(2,1); 
            else
               index_j = rand_index(2,2);
            end
            
            perturbance = zeros(3*numRobots, 3*numRobots);
            perturbance(index_i, :) = roll;
            perturbance(index_j, :) = pitch;
            for  k = 1:numRobots
                perturbance((k-1)*3+1 : (k-1)*3+3, (k-1)*3+1 : (k-1)*3+3) = zeros(3,3);
            end
            perturbance = 0.5 * (perturbance + perturbance');

            test = base_X + perturbance;
            Z_value(i,j) = evaluate_objective_WYJ_by_X(test, auxiliary_data) + trace(C*test);
        end
    end

   
    
    x = linspace(0,1,resolution); %设置x轴的范围
    y = x; %设置y轴范围
    [X,Y] = meshgrid(x,y); %将其x，y轴网格化
%     Fig = mesh(X,Y,Z_value);
    
    s = surf(X,Y,Z_value,'FaceAlpha',0.8);
    s.EdgeAlpha = 0.5;
    hold on;
    scatter3(X(resolution/2,resolution/2),Y(resolution/2,resolution/2),Z_value(resolution/2,resolution/2),'filled','MarkerFaceColor',[1 0 0], 'SizeData',500)
    hold off;
