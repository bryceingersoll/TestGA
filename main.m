clear; clc; close all;

xmin = 5;
xmax = 10;
ymin = 5;
ymax = 10;

% initialize population
pop_size = 50;

x = rand(pop_size,1)*2*xmax + xmin;
y = rand(pop_size,1)*2*ymax + ymin;

init_pop = [x,y];

% 100 or 1000 points
mating_pool = zeros(pop_size,2);

generations = 5000;

tic

for j = 1 : generations
    
    % generate mating pool
    for i = 1 : pop_size
        
        parent1 = ceil( rand(1)*pop_size);
        parent2 = ceil( rand(1)*pop_size);
        
        mating_pool(i,:) = [parent1, parent2];
        
    end
    
    % generate offspring
    offspring = zeros(pop_size*2,1);
    
    for i = 1 : pop_size
        
        offspring(2*i-1,1) = (init_pop(mating_pool(i,1),1) + init_pop(mating_pool(i,2),1))/2.0; % x value
        offspring(2*i-1,2) = (init_pop(mating_pool(i,1),2) + init_pop(mating_pool(i,2),2))/2.0; % y value
        
        parent1_value = Rosenbrock(init_pop(mating_pool(i,1),1),init_pop(mating_pool(i,1),2));
        parent2_value = Rosenbrock(init_pop(mating_pool(i,2),1),init_pop(mating_pool(i,2),2));
        
        if parent1_value < parent2_value
            offspring(2*i,1) = 0.75*init_pop(mating_pool(i,1),1) + 0.25*init_pop(mating_pool(i,2),1); % x value
            offspring(2*i,2) = 0.75*init_pop(mating_pool(i,1),2) + 0.25*init_pop(mating_pool(i,2),2); % y value
        else
            offspring(2*i,1) = 0.25*init_pop(mating_pool(i,1),1) + 0.75*init_pop(mating_pool(i,2),1); % x value
            offspring(2*i,2) = 0.25*init_pop(mating_pool(i,1),2) + 0.75*init_pop(mating_pool(i,2),2); % y value
        end
        
    end
    
    % allow for mutation
    for i = 1 : length(offspring)
        
        offspring(i,1) = offspring(i,1)*1.0 + rand(1)*0.05-0.025; %x value
        offspring(i,2) = offspring(i,1)*1.0 + rand(1)*0.05-0.025; %y value
        
    end
    
    %compute offspring's fitness
    for i = 1 : length(offspring)
        
        offspring_fitness(i) = Rosenbrock(offspring(i,1),offspring(i,2));
        
    end
    
    os_and_fitness = [offspring, offspring_fitness'];
    
    sorted_os_and_fitness = sortrows(os_and_fitness, 3);
    
    init_pop = sorted_os_and_fitness(1:pop_size,1:2);
    
    if j == generations
        
%         x1 = linspace(xmin,xmax,1000);
%         y1 = linspace(ymin,ymax,1000);
        x1 = linspace(-5,5,1000);
        y1 = linspace(-5,5,1000);
        F = zeros(length(x1),length(y1));
        for k = 1 : length(x1)
            for l = 1 : length(y1)
                F(k,l) = Rosenbrock(x1(k), y1(l));
            end
        end
        
        figure(j)
        hold on
        
        contour(x1,y1,F);
        %plot(sorted_os_and_fitness(:,1),sorted_os_and_fitness(:,2),sorted_os_and_fitness(:,3),'ok');
        plot(sorted_os_and_fitness(:,1),sorted_os_and_fitness(:,2),'ok');
        hold off
        
    end
    
end

toc

% plot ex


% x = -20 : 0.1 :10;
%
% for i = 1 : length(x)
%
%     F(i) = exp(x(i)) + x(i)^2;
% end
%
% plot(x,F)