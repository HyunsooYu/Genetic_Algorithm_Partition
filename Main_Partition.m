%Evolutionary Computation Assignment
clc
clear

fileLoc = 'F:\0_학교 일\2021-1\진화연산\Hw1\';
X = dir([fileLoc 'input*']);
degree={};
edge = {};
Best_scores=0;Avg_scores=0;score_stds=0;run_times=0; max_run=0;
Tn=table(Best_scores, Avg_scores, score_stds, run_times, max_run);
%% File Read
for i = 1:4
    name = X(i).name;
    if all(name(end-3:end)~= '.txt')
        movefile([fileLoc name],[fileLoc name '.txt'])
    end
    X = dir([fileLoc 'input*']);
    name = X(i).name;
    
    f_ID = fopen([fileLoc name]);
    while(1)
        text = fgetl(f_ID);
        if text == -1; break; end
        degree{i,str2num(text(1:5))} = str2num(text(27:29));
        edge{i, str2num(text(1:5))} = str2num(text(29:end));
    end
end
for j = [1, 2] % 1: min, 2 : max
    for K = [2 32]
        for i = 1:4 % input A, B, C, D
            %% Representation
            pop_size=150;
            Chromosome=cell(1,pop_size);
            for pop = 1:pop_size % Population
                if j==1 % min condition
                    rand_arr = [];
                    for joint = 1:floor(500/K)
                        rand_arr = [rand_arr randperm(K)];
                    end
                    rand_arr = [rand_arr randperm(mod(500,K))];
                    rand_arr = rand_arr(randperm(500));
                elseif j==2 % max condition
                    rand_arr = [];
                    for node = 1:500
                        rand_arr = [rand_arr randi(K,1,1)];
                    end
                end
                Chromosome{pop}=rand_arr;
            end
            %% Compute Fitness
            edge0={edge{i,:}};
            Fitness = Compute_fitness(Chromosome, edge0, K );
            max_run=100;
            Best_score=zeros(1,max_run);
            Avg_score=zeros(1,max_run);
            score_std=zeros(1,max_run);
            r_time=zeros(1,max_run);
            curv1 = zeros(1,max_run);
            curv2 = zeros(1,max_run);
            for run_time = 1:max_run
                tic
                %% Selection - Roulette wheel selection
                Prob = Fitness/sum(Fitness);
                num_Xover=pop_size/2;
                next_pop=cell(1,num_Xover);
                for num = 1:num_Xover %number of Crossover
                    parind = datasample(1:length(Prob),2,'Weights', Prob);
                    par1=Chromosome{parind(1)};
                    par2=Chromosome{parind(2)};
                    %% Crossover - 13 point crossover
                    div_point = sort(randi(length(par1)-1, [1,13]));
                    offspring=[par1(1:div_point(1)) par2(div_point(1)+1:div_point(2)) ...
                        par1(div_point(2)+1:div_point(3)) par2(div_point(3)+1:div_point(4) )...
                        par1(div_point(4)+1:div_point(5)) par2(div_point(5)+1:div_point(6) )...
                        par1(div_point(6)+1:div_point(7)) par2(div_point(7)+1:div_point(8) )...
                        par1(div_point(8)+1:div_point(9)) par2(div_point(9)+1:div_point(10) )...
                        par1(div_point(10)+1:div_point(11)) par2(div_point(11)+1:div_point(12) )...
                        par1(div_point(12)+1:div_point(13)) par2(div_point(13)+1:end )];

                    
                    %                     %% Crossover - One point crossover
                    %                     div_point = randi(length(par1)-1);
                    %                     offspring=[par1(1:div_point) par2(div_point+1:end) ];
                    
                    %                     %% Crossover - 3 point crossover
                    %                     div_point = sort(randi(length(par1)-1, [1,3]));
                    %                     offspring=[par1(1:div_point(1)) par2(div_point(1)+1:div_point(2)) ...
                    %                         par1(div_point(2)+1:div_point(3)) par2(div_point(3)+1:end) ];
                    
                    %                     %% Crossover - 21 point crossover
                    %                     div_point = sort(randi(length(par1)-1, [1,21]));
                    %                     offspring=[par1(1:div_point(1)) par2(div_point(1)+1:div_point(2)) ...
                    %                         par1(div_point(2)+1:div_point(3)) par2(div_point(3)+1:div_point(4) )...
                    %                         par1(div_point(4)+1:div_point(5)) par2(div_point(5)+1:div_point(6) )...
                    %                         par1(div_point(6)+1:div_point(7)) par2(div_point(7)+1:div_point(8) )...
                    %                         par1(div_point(8)+1:div_point(9)) par2(div_point(9)+1:div_point(10) )...
                    %                         par1(div_point(10)+1:div_point(11)) par2(div_point(11)+1:div_point(12) )...
                    %                         par1(div_point(12)+1:div_point(13)) par2(div_point(13)+1:div_point(14) )...
                    %                         par1(div_point(14)+1:div_point(15)) par2(div_point(15)+1:div_point(16) )...
                    %                         par1(div_point(16)+1:div_point(17)) par2(div_point(17)+1:div_point(18) )...
                    %                         par1(div_point(18)+1:div_point(19)) par2(div_point(19)+1:div_point(20) )...
                    %                         par1(div_point(20)+1:div_point(21)) par2(div_point(21)+1:end )];
                    
                    
                    
                    %% Mutation - Repair & Swap
                    % Repair
                    if j==1 % min condition
                        assign=zeros(1,32);
                        for m = 1:K % Repair on m-th subset
                            assign(m) = [length(find(offspring==m))];
                        end
                        m_flag=0;
                        for m=1:K % Fix Unbalanced condition
                            if ~( (assign(m)==ceil(length(offspring)/K) )||(assign(m)==floor(length(offspring)/K) ))
                                m_flag=1; % Unbalanced!
                            end
                        end
                        if m_flag
                            for m=1:K
                                if assign(m)>ceil(length(offspring)/K)
                                    y=randsample(find(offspring==m), assign(m)-ceil(length(offspring)/K) );
                                    for yy=1:length(y)
                                        offspring(yy)=randsample(find(assign==min(assign)),1);
                                    end
                                end
                            end
                            if length(find(assign==ceil(length(offspring)/K)))>20
                                y=randsample(find(assign==ceil(length(offspring)/K)), ...
                                    length(find(assign==ceil(length(offspring)/K)))-20);
                                for yy=1:length(y)
                                    offspring(yy)=randsample(find(assign==min(assign)),1);
                                end
                            end
                        end
                    end
                    % Swap mutation for permutation encoding
                    mu_rand = rand;
                    if mu_rand<0.001
                        y=randsample(length(offspring),2);
                        a1=offspring(y(1));
                        a2=offspring(y(2));
                        offspring(y(1))=a2;
                        offspring(y(2))=a1;
                    end
                    next_pop{num}=offspring;
                end
                %% Replacement - GENITOR style (Replace Worst)
                next_Fitness = Compute_fitness(next_pop, edge0, K);
                if j==2 % Replace worst on Maximization
                    while( min(Fitness) <max(next_Fitness))
                        [Fm1 Fm2] = min(Fitness);
                        [Fn1 Fn2] = max(next_Fitness);
                        Fitness(Fm2)=Fn1;
                        Chromosome{Fm2}=next_pop{Fn2};
                        next_Fitness(Fn2)=Fm1;
                    end
                elseif j==1 % Replace worst on minimization
                    while( max(Fitness) >min(next_Fitness))
                        [Fm1 Fm2] = max(Fitness);
                        [Fn1 Fn2] = min(next_Fitness);
                        Fitness(Fm2)=Fn1;
                        Chromosome{Fm2}=next_pop{Fn2};
                        next_Fitness(Fn2)=Fm1;
                    end
                end
                Final_Fitness = Compute_fitness(Chromosome, edge0, K);
                if j==2; Best_score(run_time)=max(Final_Fitness); %Maximization
                elseif j==1; Best_score(run_time)=min(Final_Fitness);end %minimization
                
                Avg_score(run_time)=mean(Final_Fitness);
                score_std(run_time)=std(Final_Fitness);
                curv1(run_time)= mean(Final_Fitness)+std(Final_Fitness);
                curv2(run_time)= mean(Final_Fitness)-std(Final_Fitness);
                r_time(run_time) = toc;
                rr_time = r_time(run_time)
                [run_time j K i]
                max(Final_Fitness);
            end
            plot(1:max_run, Best_score, 'b', 'LineWidth', 1.5);
            title('Best Score');
            xlabel('run time'); ylabel('score');
            saveas(gcf, [num2str(j) '_'  num2str(K) '_' num2str(i) '_Best score.png' ])
            close all
            
            plot(1:max_run, Avg_score);
            title('Average Score');
            xlabel('run time'); ylabel('score');
            saveas(gcf, [num2str(j) '_' num2str(K) '_' num2str(i) '_' 'Average score.png' ])
            close all
            
            x1=1:max_run;
            x2 = [x1, fliplr(x1)];
            inBetween = [curv1, fliplr(curv2)];
            patch(x2, inBetween, 'b', 'FaceAlpha', .3, 'LineStyle', 'none');hold on;
            plot(1:max_run, Avg_score, 'b', 'LineWidth', 1.5);
            xlabel('run time'); ylabel('score');
            saveas(gcf, [num2str(j) '_' num2str(K) '_' num2str(i) '_' 'Avg std score.png' ])
            close all
            
            Best_scores=Best_score(max_run);Avg_scores=Avg_score(max_run);score_stds=score_std(max_run);run_times=mean(r_time);
            max_run = max(r_time);
            Tnn = table(Best_scores, Avg_scores, score_stds, run_times, max_run);
            Tn=[Tn; Tnn];
        end
    end
end
save('Partition_GA_13_pop150', 'Tn');

