function [num_FS, population] = MOEADMP_UDeV62(MOP,MaxEv,N)

    % MOP:multiobjective problem
    % SC:stopping criterion
    % N:the number of the subproblems considered in MOEA/D
    %增加外部种群
    %去除约束冲突更优的解更高概率被选中
    %向量改为多方向量
    %修改替换规则
    %添加定期删除规则
    %调整参数
    %修正CV_old和CV_new

%     s = RandStream.create('mrg32k3a','NumStreams',30,'StreamIndices',7);
%     RandStream.setGlobalStream(s);
    population = Init_Pop(MOP,N);
%     rate_T = 0.1;
%     T = ceil(N * rate_T);     
    T = 3;
    D =MOP.D;
    ub = MOP.max;
    lb = MOP.min;
%     MaxEP_p = N;
    H = 100;
    k1 = 1;
    k2 = 1; 
    k3 = 1;
    gen = 0;
    N_remove = 700;
    
    update_pop = 0.5;%按概率更新个体
    
%     sel_neighbor = 0.2;
    cal_pop = 0;%初始评估个体数cal_sc   

    % 生成均匀分布权重向量，存于W
    M = MOP.M;%M为目标函数维度，也就是目标函数包含的“小”目标个数
    DM = MOP.DM;
    W1 = WV_generator2(N^(1/2),M/DM);
    W = [];
    for (i = 1:1:N^(1/2))
        W = [W;[repmat(W1(i,:),N^(1/2),1),W1]];
    end
    N = size(W,1);
    
    population = population(1:N,:);
%     p_best = [0.571384707883747,-0.137992429952902,0.211521801089946,0.325427667152480,0.0491621364619185,0.103250940645195,0.175939433704619,-0.00159295101894481,-0.00934033659640696,0.0136105285512664,0.00220397432508599,-0.000669548569235403,-0.0589757404798028,-0.0685103826270797,-0.0684732648486362,-0.192160761036174,-0.178476597817853,-0.472488872945460,-0.352188232854153,-0.599297535077626,-0.279797555463019,-0.286670584282735,-0.367780510354391,-0.417517097155888,-0.385321352719058,-0.230749176626685,0.122840527279066,0.691971377574498,0.000391125000474750,0.000496723027696569,0.0619439842338422,0.995460844249653,0.996826712269382,0.913217737333581];

    FV = MOP.CalV(population);
    
    h_x = MOP.Cal_Con_h(population);
%     U_max = max(h_x,[],1);
%     eps_0 = max(sum(h_x.*(1./U_max),2)/sum(1./U_max));
%     eps_0 = max(sum(h_x,2));
%     alp = (-log(eps_0)-6)/log(1-gam);

    eps = 0;
        
    %初始化参考界面Z
    Z = min(FV,[],1);
    
    % 计算各权重向量相互间欧几里得距离,并找出每个向量距离最近的T个邻居向量（包括自身）
    D_of_W = pdist2(W,W);

    B = zeros(N,T);
    TempB = zeros(N,N);
    
    for i = 1:N
        [~,TempB(i,:)] = sort(D_of_W(i,:),2);
        B(i,:) = TempB(i,1:T);
    end
    
    eps_x = sum(h_x,2);
    CV_old = eps_x;
    Pi = ones(N,1);
    DELTA = zeros(N,1);
    EP_p = [];%外部种群保留近期被淘汰解
    
    % 初始化F, CR参数
    M_CR1 = 0.5*ones(1,H); M_CR3 = 0.5*ones(1,H);
    CR = ones(3,N);
    M_F1 = 0.5*ones(1,H); M_F2 = 0.5*ones(1,H); M_F3 = 0.5*ones(1,H);
    F = ones(3,N);
    
    % 开始循环
    %while gen < MaxG
    while cal_pop < MaxEv
        gen = gen+1;
        
        if rem(gen,1) == 0
            b_e = min(eps_x);
            %fprintf('评估代数：%d；最佳约束冲突值：%d\n',gen,b_e);
            fprintf('评估个体数：%d；最佳约束冲突值：%d\n',cal_pop,b_e);
        end
        
        %计算变化率
        if rem(gen,50) == 0
            CV_new = eps_x;
            for i = 1:N
                if CV_old(i,:) == 0
                    Pi(i) = 1;
                else
                    DELTA(i) = (CV_old(i,:) -  CV_new(i,:))/CV_old(i,:);
                    if DELTA(i) > 0.001
                        Pi(i) = 1;
                    else
                        Pi(i) = (0.95 + 0.05*DELTA(i)/0.001)*Pi(i);
                    end
                end
            end
            CV_old = CV_new;
        end
        
        
        %每隔一定代数，去掉CV变化率小的个体
        if gen >= 200 && rem(gen,50) == 0 && N > 200
            
%             CV_new = eps_x;
             
            [~,d_CV_rank] = sort(Pi);
            remove_ID = d_CV_rank(1:N_remove);
            p_remain = true(1,N);
            p_remain(remove_ID) = false;
            EP_p = [EP_p;population(remove_ID,:)];
            
            population = population(p_remain,:);
            FV = FV(p_remain,:);
            h_x = h_x(p_remain,:);
            eps_x = eps_x(p_remain,:);
            W = W(p_remain,:);
            Pi = Pi(p_remain,:);
            DELTA = DELTA(p_remain,:);
            
            [N,~] = size(population);
            
            D_of_W = pdist2(W,W);
            B = zeros(N,T);
            TempB = zeros(N,N);
            
            for i = 1:N
                [~,TempB(i,:)] = sort(D_of_W(i,:),2);
                B(i,:) = TempB(i,1:T);
            end
            
        end
        
        
%         debug_Pic = 1;
%         if debug_Pic
%             P_draw = population(:,[1 2]);
%             %             PF_each_p1 = PF_all{1,1};
%             %             PF_each_p2 = PF_all{1,2};          
%             plot(...
%                 P_draw(:,1),P_draw(:,2),'ob'...
%                 ,p_best(:,1),p_best(:,2),'or'...
%                 ,'markersize',10.0...
%                 )%,FV_best(:,1),FV_best(:,2),'py',FV_best(:,3),FV_best(:,4),'pm');
%             axis([-1 1 -1 1]);
%             title(cal_sc);
%             drawnow
%         end

%         gen = gen + 1;
%         if gen/max_gen <= gam
%             eps = eps_0*((1 - gen/max_gen)^alp);
%         else
%             
% %             population(200,:) = p_best;
%         end
        
%         prob = (N-1:-1:0)/N;
        SCR1 = zeros(N,1); SCR3 = zeros(N,1);
        SF1 = zeros(N,1); SF2 = zeros(N,1); SF3 = zeros(N,1);
        W1 = zeros(N,1); W2 = zeros(N,1); W3 = zeros(N,1);
        id1 = 0; id2 = 0; id3 = 0;
        
        
        
        MaxEP_p = N;
        
        EPp_size = size(EP_p,1);
        if (EPp_size > MaxEP_p)
            Select_N = randperm(EPp_size,MaxEP_p);
            EP_p = EP_p(Select_N,:);
            EPp_size = MaxEP_p;
        end

        [Num_of_P,~] = size(population);
        cal_pop = cal_pop + 3*Num_of_P;
%         [~,sortP] = sort(eps_x);
%         [~,sortsortP] = sort(sortP);
        for i = 1:N               

            % DE/rand/1 %%%%%%%%%%%%%%%%%%%
            r1 = randi(N);
%             while rand >= prob(sortsortP(r1)) || r1 == i
%                 r1 = randi(N);
%             end
            while r1 == i
                r1 = randi(N);
            end 

            r2 = randi(N);
%             while rand >= prob(sortsortP(r2)) || r2 == r1 || r2 == i
%                 r2 = randi(N);
%             end
            while r2 == r1 || r2 == i
                r2 = randi(N);
            end

            r3 = randi(N);
            while r3 == r1 || r3 == r2 || r3 == i
                r3 = randi(N);
            end

            temp_rand = randi(H);
            CR(1,i) = normrnd(M_CR1(temp_rand), 0.1);
            CR(1,i) = min(CR(1,i), 1); CR(1,i) = max(CR(1,i), 0);

            F(1,i) = cauchyrand(M_F1(temp_rand), 0.1);
            while F(1,i)<=0
                F(1,i) = cauchyrand(M_F1(temp_rand), 0.1);
            end
            F(1,i) = min(F(1,i), 1);

            % Mutation
            S1_Pop(i,:) = population(r1,:) + F(1,i)*(population(r2,:)-population(r3,:));


            % Binomial Crossover
            j_rand = randi(D);
            temp = rand(1, D);
            U1 = zeros(1,D);
            U1(temp <= CR(1,i)) = S1_Pop(i,temp <= CR(1,i));
            U1(j_rand) = S1_Pop(i,j_rand);
            U1(temp > CR(1,i)) = population(i, temp > CR(1,i));

            S1_Pop(i,:) = U1;
            S1_Pop(i,:)=min(max(S1_Pop(i,:),lb),ub);
            S1_Pop(i,:) = Polynomial_Mutation(S1_Pop(i,:),ub,lb,D);
            eps_S1Pop =sum(MOP.Cal_Con_h(S1_Pop(i,:)),2);               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % DE/current-to-rand/1 %%%%%%%%
            r1 = randi(N);
%             while rand >= prob(sortsortP(r1)) || r1 == i
%                 r1 = randi(N);
%             end
            while r1 == i
                r1 = randi(N);
            end

            r2 = randi(N);
%             while rand >= prob(sortsortP(r2)) || r2 == r1 || r2 == i
%                 r2 = randi(N);
%             end
            while r2 == r1 || r2 == i
                r2 = randi(N);
            end

            r3 = randi(N);
            while r3 == r1 || r3 == r2 || r3 == i
                r3 = randi(N);
            end

            temp_rand = randi(H);
            F(2,i) = cauchyrand(M_F2(temp_rand), 0.1);
            while F(2,i)<=0
                F(2,i) = cauchyrand(M_F2(temp_rand), 0.1);
            end
            F(2,i) = min(F(2,i), 1);

            % Mutation
            S2_Pop(i,:) = population(i,:) + F(2,i)*(population(r1,:)-population(i,:)) + F(2,i)*(population(r2,:)-population(r3,:));
            S2_Pop(i,:)=min(max(S2_Pop(i,:),lb),ub);
            S2_Pop(i,:) = Polynomial_Mutation(S2_Pop(i,:),ub,lb,D);
            eps_S2Pop =sum(MOP.Cal_Con_h(S2_Pop(i,:)),2); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % DE/current-to-pbest/1 %%%%%%%
            a = 5; b = 20;
            p = round((b-a)*rand + a);
            [~,choose_p] = sort(eps_x);%需要根据CV排序，选出前p%个体

            r1 = randi(p);
            while choose_p(r1) == i
                r1 = randi(p);
            end

            r2 = randi(N);
%             while rand >= prob(sortsortP(r2)) || r2 == r1 || r2 == i
%                 r2 = randi(N);
%             end
            while r2 == r1 || r2 == i
                r2 = randi(N);
            end                

            r3 = randi(N+EPp_size);
            while r3 == r1 || r3 == r2 || r3 == i
                r3 = randi(N+EPp_size);
            end


            temp_rand = randi(H);
            CR(3,i) = normrnd(M_CR3(temp_rand), 0.1);
            CR(3,i) = min(CR(3,i), 1); CR(3,i) = max(CR(3,i), 0);

            F(3,i) = cauchyrand(M_F3(temp_rand), 0.1);
            while F(3,i)<=0
                F(3,i) = cauchyrand(M_F3(temp_rand), 0.1);
            end
            F(3,i) = min(F(3,i), 1);

            % Mutation
            if (r3<=N)
                S3_Pop(i,:) = population(i,:) + F(3,i)*(population(choose_p(r1),:)-population(i,:)) + F(3,i)*(population(r2,:)-population(r3,:));
            else
                S3_Pop(i,:) = population(i,:) + F(3,i)*(population(choose_p(r1),:)-population(i,:)) + F(3,i)*(population(r2,:)-EP_p(r3-N,:));
            end

            % Binomial Crossover
            j_rand = randi(D);
            temp = rand(1, D);
            U1 = zeros(1,D);
            U1(temp <= CR(3,i)) = S3_Pop(i,temp <= CR(3,i));
            U1(j_rand) = S3_Pop(i,j_rand);
            U1(temp > CR(3,i)) = population(i, temp > CR(3,i));

            S3_Pop(i,:) = U1;
            S3_Pop(i,:)=min(max(S3_Pop(i,:),lb),ub);
            S3_Pop(i,:) = Polynomial_Mutation(S3_Pop(i,:),ub,lb,D);
            eps_S3Pop =sum(MOP.Cal_Con_h(S3_Pop(i,:)),2); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            F_S1 = MOP.CalV(S1_Pop(i,:));
            F_S2 = MOP.CalV(S2_Pop(i,:));
            F_S3 = MOP.CalV(S3_Pop(i,:));

            Z = min(Z,F_S1);
            Z = min(Z,F_S2);
            Z = min(Z,F_S3);

            g_S1 = max(abs(F_S1-Z).*W(i,:));
            g_S2 = max(abs(F_S2-Z).*W(i,:));
            g_S3 = max(abs(F_S3-Z).*W(i,:));

            g_x = max(abs(FV(i,:)-Z).*W(i,:));

            if(eps_S1Pop > 0 && eps_x(i) > 0)
                if(eps_S1Pop < eps_x(i))
                    id1 = id1 + 1;
                    SCR1(id1) = CR(1,i);
                    SF1(id1) = F(1,i);
                    W1(id1) = (eps_x(i) - (eps_S1Pop));
                end
            end

            if(eps_S1Pop == 0 && eps_x(i) > 0)
                id1 = id1 + 1;
                SCR1(id1) = CR(1,i);
                SF1(id1) = F(1,i);
                W1(id1) = (eps_x(i)) - 0;
            end

            if(eps_S1Pop == 0 && eps_x(i) == 0)
                if(g_S1 < g_x)
                    id1 = id1 + 1;
                    SCR1(id1) = CR(1,i);
                    SF1(id1) = F(1,i);
                    W1(id1) = g_x-g_S1;
                end
            end

            if(eps_S2Pop > 0 && eps_x(i) > 0)
                if(eps_S2Pop < eps_x(i))
                    id2 = id2 + 1;
                    SF2(id2) = F(2,i);
                    W2(id2) = (eps_x(i)) - (eps_S2Pop);
                end
            end

            if(eps_S2Pop == 0 && eps_x(i) > 0)
                id2 = id2 + 1;
                SF2(id2) = F(2,i);
                W2(id2) = (eps_x(i)) - 0;
            end

            if(eps_S2Pop == 0 && eps_x(i) == 0)
                if(g_S2 < g_x)
                    id2 = id2 + 1;
                    SF2(id2) = F(2,i);
                    W2(id2) = g_x-g_S2;
                end
            end

            if(eps_S3Pop > 0 && eps_x(i) > 0)
                if(eps_S3Pop < eps_x(i))
                    id3 = id3 + 1;
                    SCR3(id3) = CR(3,i);
                    SF3(id3) = F(3,i);
                    W3(id3) = (eps_x(i)) - (eps_S3Pop);
                end
            end

            if(eps_S3Pop == 0 && eps_x(i) > 0)
                id3 = id3 + 1;
                SCR3(id3) = CR(3,i);
                SF3(id3) = F(3,i);
                W3(id3) = (eps_x(i)) - 0;
            end

            if(eps_S3Pop == 0 && eps_x(i) == 0)
                if(g_S3 < g_x)
                    id3 = id3 + 1;
                    SCR3(id3) = CR(3,i);
                    SF3(id3) = F(3,i);
                    W3(id3) = g_x-g_S3;
                end
            end

            eps_temp = [eps_S1Pop,eps_S2Pop,eps_S3Pop];
            g_temp = [g_S1,g_S2,g_S3];
            Pop_temp = [S1_Pop(i,:);S2_Pop(i,:);S3_Pop(i,:)];
            sort_ID = sort1(g_temp,eps_temp);          

            %选择最优的用eC规则与原个体比较

            %按概率更新个体
            offspring_y = Pop_temp(sort_ID(1),:);
            h_y = MOP.Cal_Con_h(offspring_y);
            eps_y = sum(h_y,2);                
            F_y = MOP.CalV(offspring_y);
            Z = min(Z,F_y);
            for j = 1:T
                g_x = max(abs(FV(B(i,j),:)-Z).*W(B(i,j),:));
                g_y = max(abs(F_y-Z).*W(B(i,j),:));

                if (eps_y <= eps) && (eps_x(B(i,j),:) <= eps)
                    if (g_y <= g_x) 
                        if rand < update_pop
                            EP_p = [EP_p;population(B(i,j),:)];
                            population(B(i,j),:) = offspring_y;
                            FV(B(i,j),:) = F_y;
                            h_x(B(i,j),:) = h_y;
                            eps_x(B(i,j),:) = eps_y;
                        else
                            EP_p = [EP_p;offspring_y];
                        end
                    else
                        EP_p = [EP_p;offspring_y];
                    end
                elseif (eps_y == eps_x(B(i,j),:)) 
                    if (g_y <= g_x) 
                        if rand < update_pop
                            EP_p = [EP_p;population(B(i,j),:)];
                            population(B(i,j),:) = offspring_y;
                            FV(B(i,j),:) = F_y;
                            h_x(B(i,j),:) = h_y;
                            eps_x(B(i,j),:) = eps_y;
                        else
                            EP_p = [EP_p;offspring_y];
                        end
                    else
                        EP_p = [EP_p;offspring_y];
                    end
                elseif (eps_y < eps_x(B(i,j),:))
                    if rand < update_pop
                        EP_p = [EP_p;population(B(i,j),:)];
                        population(B(i,j),:) = offspring_y;
                        FV(B(i,j),:) = F_y;
                        h_x(B(i,j),:) = h_y;
                        eps_x(B(i,j),:) = eps_y;
                    else
                        EP_p = [EP_p;offspring_y];
                    end
                elseif (eps_y > eps_x(B(i,j),:))
                    EP_p = [EP_p;offspring_y];
                end                   
            end
        end
       
        
        EP_p = unique(EP_p,'rows');
        
        SCR1 = SCR1(1:id1); SF1 = SF1(1:id1); W1 = W1(1:id1);
        SF2 = SF2(1:id2); W2 = W2(1:id2);
        SCR3 = SCR3(1:id3); SF3 = SF3(1:id3); W3 = W3(1:id3);
        
        %自适应CR与F
        if id1 > 0
            M_CR1(k1) = update_MCR(SCR1, W1);
            M_F1(k1) = update_MF(SF1, W1);
            k1 = mod(k1,H) + 1;
        end
        
        if id2 > 0
            M_F2(k2) = update_MF(SF2, W2);
            k2 = mod(k2,H) + 1;
        end
        
        if id3 > 0
            M_CR3(k3) = update_MCR(SCR3, W3);
            M_F3(k3) = update_MF(SF3, W3);
            k3 = mod(k3,H) + 1;
        end
        
        
        
    end
    
    final_EP = MOP.CalV(population);
    final_h = MOP.Cal_Con_h(population);
%     final_eps = sum(final_h.*(1./U_max),2)/sum(1./U_max);
    final_eps = sum(final_h,2);
    num_FS = sum(final_eps == 0);
    best_eps = min(final_eps);
end
                
                    
                        
    
    
