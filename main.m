clc,clear,close all
warning off
% BA算法参数
maxiter = 200;  % 迭代次数
sizepop = 10;  % 种群数量
% 频率范围
popmin1 = -1;  popmax1 = 1; % x1  频率
popmin2 = -1;  popmax2 = 1; % x2  频率
Qmin = 0.1;    % 最小频率
Qmax = 0.5;    % 最大频率

A = 0.1;          % 音量 (不变或者减小)
impluse = 0.85;   % 脉冲率 (不变或增加)

Vmin = -1; % 最小速度
Vmax = 1;  % 最大速度
%% 初始化种群
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
    V(i,1)=0;
    V(i,2)=0;
end
% 记录一组最优值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % 全局最佳
gbest=pop;                % 个体最佳
fitnessgbest=fitness;     % 个体最佳适应度值
fitnesszbest=bestfitness; % 全局最佳适应度值
%% 迭代寻优
for i=1:maxiter
    for j=1:sizepop
        Q = Qmin + (Qmax-Qmin)*rand;
        V(j,:) = V(j,:) + (pop(j,:)-zbest)*Q;
        % V--x1
        if V(j,1)>Vmax
            V(j,1)=Vmax;
        end
        if V(j,1)<Vmin
            V(j,1)=Vmin;
        end
        % V--x2
        if V(j,2)>Vmax
            V(j,2)=Vmax;
        end
        if V(j,2)<Vmin
            V(j,2)=Vmin;
        end
        
        pop(j,:) = pop(j,:) + 0.5*V(j,:);
        
        % 脉冲率
        if rand>impluse
            pop(j,:) = zbest + A * randn(1,2);
%             x1 = popmin1 + (popmax1-popmin1)*rand;
%             x2 = popmin2 + (popmax2-popmin2)*rand;
%             pop(j,:) = [x1,x2];
        end
        
        % x1  越界限制
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2  越界限制
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % 适应度更新
        fitness(j) = fun(pop(j,:));
        
        % 比较  个体间比较
        if fitness(j)<fitnessgbest(j)
            fitnessgbest(j) = fitness(j);
            gbest(j,:) = pop(j,:);
        end
        if fitness(j)<bestfitness
            bestfitness = fitness(j);
            zbest =  pop(j,:);
        end
        
    end
    fitness_iter(i) = bestfitness;
end

disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight