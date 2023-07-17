suites=[1,3,4];%,4]%,2,3,4];
runs = [541];
fig = figure(1);
fig.Renderer = 'painters';
fig.Units = 'centimeters';
fig_size_x = 9.5;
fig_size_y = 9.5;
fig.Position = [10 10 fig_size_x fig_size_y];
for j=1:length(suites)
    
    for i= 1: length(runs)
        suite = suites(j);
        run_i=runs(i);
  
        subplot(2,2,1)
        plot(s(suite).eps_g{run_i}(s(suite).index_erupt_begin{run_i}),s(suite).M_dis{run_i}./s(suite).tot_Mass{run_i}(1),'v','MarkerSize',2)
        hold on
         ylabel('M_{dis}/M_o'); xlabel(' \epsilon_g at start of eruption')
        % title([num2str(s(suite).V{run_i}(1)/1e9) ', ' num2str(s(suite).mdot_in(run_i)) ', ' num2str(run_i)])
        %axis([0.065,.085,0.009,0.013])
        
        subplot(2,2,2)
        plot((s(suite).eps_g{run_i}(s(suite).index_erupt_end{run_i})-s(suite).eps_g{run_i}(s(suite).index_erupt_begin{run_i}))./(1),s(suite).M_dis{run_i}./s(suite).tot_Mass{run_i}(1),'v','MarkerSize',2)
        hold on
         ylabel('M_{dis}/M_o'); xlabel('Change in \epsilon_g during eruption')
        % title([num2str(s(suite).V{run_i}(1)/1e9) ', ' num2str(s(suite).mdot_in(run_i)) ', ' num2str(run_i)])
        axis([0.0045,0.0055,0.009,0.013])
        % legend('0.5% H_2O, 500 ppm CO_2','0.5% H_2O, 10,000 ppm CO_2','2% H_2O, 500 ppm CO_2')
    
    end
%         suite = suites(2);
%         run_i=runs(1);
%         
%         subplot(2,2,3)
%         plot(s(suite).eps_g{run_i}(s(suite).index_erupt_begin{run_i}),s(suite).M_dis{run_i}./s(suite).tot_Mass{run_i}(1),'.')
%        % hold on
%         %ylabel('M_{dis}/M_o'); xlabel(' \epsilon_g at start of eruption')
%        % title([num2str(s(suite).V{run_i}(1)/1e9) ', ' num2str(s(suite).mdot_in(run_i)) ', ' num2str(run_i)])
%         %axis([-inf,inf,-inf,inf])
%         
%         subplot(2,2,4)
%         plot(s(suite).eps_g{run_i}(s(suite).index_erupt_end{run_i})-s(suite).eps_g{run_i}(s(suite).index_erupt_begin{run_i}),s(suite).M_dis{run_i}./s(suite).tot_Mass{run_i}(1),'.')
%         %hold on
%         %ylabel('M_{dis}/M_o'); xlabel('Change in \epsilon_g during eruption')
%        % title([num2str(s(suite).V{run_i}(1)/1e9) ', ' num2str(s(suite).mdot_in(run_i)) ', ' num2str(run_i)])
%         %axis([-0.001,0.001,-inf,inf])
%         %legend('0.5% H_2O, 500 ppm CO_2','0.5% H_2O, 10,000 ppm CO_2','2% H_2O, 500 ppm CO_2')
    end