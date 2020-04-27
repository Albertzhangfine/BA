clc,clear,close all
warning off
% BA�㷨����
maxiter = 200;  % ��������
sizepop = 10;  % ��Ⱥ����
% Ƶ�ʷ�Χ
popmin1 = -1;  popmax1 = 1; % x1  Ƶ��
popmin2 = -1;  popmax2 = 1; % x2  Ƶ��
Qmin = 0.1;    % ��СƵ��
Qmax = 0.5;    % ���Ƶ��

A = 0.1;          % ���� (������߼�С)
impluse = 0.85;   % ������ (���������)

Vmin = -1; % ��С�ٶ�
Vmax = 1;  % ����ٶ�
%% ��ʼ����Ⱥ
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
    V(i,1)=0;
    V(i,2)=0;
end
% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % ȫ�����
gbest=pop;                % �������
fitnessgbest=fitness;     % ���������Ӧ��ֵ
fitnesszbest=bestfitness; % ȫ�������Ӧ��ֵ
%% ����Ѱ��
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
        
        % ������
        if rand>impluse
            pop(j,:) = zbest + A * randn(1,2);
%             x1 = popmin1 + (popmax1-popmin1)*rand;
%             x2 = popmin2 + (popmax2-popmin2)*rand;
%             pop(j,:) = [x1,x2];
        end
        
        % x1  Խ������
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2  Խ������
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % ��Ӧ�ȸ���
        fitness(j) = fun(pop(j,:));
        
        % �Ƚ�  �����Ƚ�
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

disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight