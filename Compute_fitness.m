
function Fitness = Compute_fitness(Chromosome, edge, K)
Fitness=cell(1,length(Chromosome));
for n=1:length(Chromosome)
    arr = Chromosome{1,n};  stack=[]; score=0;
    for sub = K:-1:1
        subset = find(arr==sub); %current subset component
        for ind = subset
            stack=[stack edge{ind}];%target subset component
            stack_n =[];
            if sub~=K
                for ind_n = prev_subset
                    stack_n = [stack_n edge{ind_n}];
                end
            end
            stack=setdiff(stack,stack_n); %Diffrence set(stack-stack_n) 
        end
        stack = setdiff(stack,subset);
        score = score + length(stack);
        prev_subset = subset;
    end
    Fitness{n} = score;
end
Fitness=cell2mat(Fitness);
end