function offspring = Polynomial_Mutation(x,ub,lb,D)

    r_pm = zeros(1,D);
    pm = 1/D;
    dis_idx = 20;
    length_bound = ub - lb;
    
    %polynomial mutation
        
%     for i = 1:D
%         r_rand = rand;
%         if r_rand < 0.5
%             r_pm(:,i) = (2 * r_rand)^(1/(dis_idx + 1))-1;
%         else
%             r_pm(:,i) = 1-(2-2 * r_rand)^(1/(dis_idx + 1));
%         end
%     end
% 
%     y_pm = r_pm .* length_bound;
%     pm_idx = rand(1,D) < pm;
%     x(:,pm_idx) = x(:,pm_idx) + y_pm(:,pm_idx);

    pm_idx = rand(1,D) < pm;
    r_rand = rand(1,D);

    for i = 1:D
        if r_rand(1,i) <= 0.5
            r_pm(:,i) = (2*r_rand(1,i) + (1-2*r_rand(1,i))*(1-(x(:,i)-lb(:,i))/length_bound(:,i))^(dis_idx+1))^(1/(dis_idx+1))-1;
        else
            r_pm(:,i) = 1-(2*(1-r_rand(1,i)) + 2*(r_rand(1,i)-0.5)*(1-(ub(:,i)-x(:,i))/length_bound(:,i))^(dis_idx+1))^(1/(dis_idx+1));
        end
    end

    y_pm = r_pm .* length_bound;
    x(:,pm_idx) = x(:,pm_idx) + y_pm(:,pm_idx);   

    %处理越界问题
    x=(x<lb).*(lb)+(x>=lb).*x;
    x=(x>ub).*(ub)+(x<=ub).*x;
        
    offspring = x;
    
end