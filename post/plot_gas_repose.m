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
        start_repose_M=[s(suite).tot_Mass{run_i}(1);s(suite).tot_Mass{run_i}(s(suite).index_erupt_end{run_i}(1:end-1))];
        end_repose_M=[s(suite).tot_Mass{run_i}(s(suite).index_erupt_begin{run_i})];
        dM_repose=end_repose_M-start_repose_M;
        eps_g_repose_start =[s(suite).eps_g{run_i}(1);s(suite).eps_g{run_i}(s(suite).index_erupt_end{run_i}(1:end-1))];
        eps_g_repose_end =[s(suite).eps_g{run_i}(s(suite).index_erupt_begin{run_i}(1:end))];
        deps_g_repose =eps_g_repose_end- eps_g_repose_start;
        eps_g_repose_avg =(eps_g_repose_end+ eps_g_repose_start)./2;
        
        eps_x_repose_start =[s(suite).eps_x{run_i}(1);s(suite).eps_x{run_i}(s(suite).index_erupt_end{run_i}(1:end-1))];
        eps_x_repose_end =[s(suite).eps_x{run_i}(s(suite).index_erupt_begin{run_i}(1:end))];
        deps_x_repose =eps_x_repose_end- eps_x_repose_start;
        
        T_repose_start =[s(suite).T{run_i}(1);s(suite).T{run_i}(s(suite).index_erupt_end{run_i}(1:end-1))];
       T_repose_end =[s(suite).T{run_i}(s(suite).index_erupt_end{run_i}(1:end))];
        T_repose =T_repose_end- T_repose_start;
        
 X_co2_repose_start =[s(suite).X_co2{run_i}(1);s(suite).X_co2{run_i}(s(suite).index_erupt_end{run_i}(1:end-1))];
        X_co2_repose_end =[s(suite).X_co2{run_i}(s(suite).index_erupt_begin{run_i}(1:end))];
        dX_co2_repose =X_co2_repose_end- X_co2_repose_start;
        X_co2_repose_avg =(X_co2_repose_end+ X_co2_repose_start)./2;

        
        subplot(2,2,3)
        plot( eps_g_repose_start, dM_repose./s(suite).tot_Mass{run_i}(1),'v','MarkerSize',2)
        hold on
       % plot( eps_g_repose_avg, X_co2_repose_end,'.')
        
       % ylabel('M_{dis}/M_o'); xlabel(' \epsilon_g at start of eruption')
       % title([num2str(s(suite).V{run_i}(1)/1e9) ', ' num2str(s(suite).mdot_in(run_i)) ', ' num2str(run_i)])
       axis([0.07,.085,0.0105,0.0145])
       xlabel("g beginning repose")
        subplot(2,2,4)
        plot(deps_g_repose,dM_repose./s(suite).tot_Mass{run_i}(1),'v','MarkerSize',2)
        hold on
        %plot( deps_g_repose, X_co2_repose_end,'.')
        xlabel("\Delta g fraction repose")
         axis([-0.0055,-0.0045,0.0105,0.0145])
        legend('4% H_2O, 100 ppm CO_2','4% H_2O, 1000 ppm CO_2','6% H_2O, 100 ppm CO_2')
    end
end