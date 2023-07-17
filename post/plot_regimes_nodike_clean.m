%% Regime diagram plotting overview
% Always run this top section before running sections of the code
% List of plots (in order):
% General behavior, eruptive lifespan, number of eruptions, 
% eruption frequency, avg eruption volume, mass intrusive to extrusive ratio, 
% volume intrusive to extrusive ratio, mass change (Mfinal/M0), 
% average growth rate (mass), normalized average growth rate (mass),
% volume change (Vfinal/V0), average growth rate (volume), longevity,
% longevity comparison, mass change comparison

save_figures=1;
close all
composition= 'mafic'; %other option is silicic
if save_figures == 1
    folder_name = ['regimeDiagrams_' composition '/'];
    mkdir(folder_name)
end
conditions=[4];pointsize = 100;fig_size = 800;ftsize  = 14;
xlim_left = 2e-1;xlim_right = 2e2;
ylim_bottom = 2e-3;ylim_top = 1e2;
%% General behavior
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);   
    
        ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    ind_grow_erupt=[];
    for i=1:length(ind_erupt_events)
        j=ind_erupt_events(i);
        if s(h).mass_change(j)>1
            ind_grow_erupt = [ind_grow_erupt j];
        end
    end
    
        ind_grow_noerupt=[];
    for i=1:length(ind_no_erupt)
        j=ind_no_erupt(i);
        if s(h).mass_change(j)>1
            ind_grow_noerupt = [ind_grow_noerupt j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_no_erupt), s(h).theta2(ind_no_erupt),...
        pointsize, log10(s(h).num_eruption_events(ind_no_erupt)), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).num_eruption_events(ind_erupt_nolimit)), 'bo', 'filled', 'MarkerEdgeColor', 'k');    
    hand(3) = scatter(s(h).theta1(ind_grow_noerupt), s(h).theta2(ind_grow_noerupt),...
        pointsize, log10(s(h).num_eruption_events(ind_grow_noerupt)), 'mo', 'filled', 'MarkerEdgeColor', 'k');    
    hand(4) = scatter(s(h).theta1(ind_grow_erupt), s(h).theta2(ind_grow_erupt),...
        pointsize, log10(s(h).num_eruption_events(ind_grow_erupt)), 'ro', 'filled', 'MarkerEdgeColor', 'k');     
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    title(['No erupt, erupt, grow & erupt', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' MPa maxn 500'])    
    colorbar 
    cmin = 0;
    cmax = log10(500); 
    colormap(jet);
    caxis([cmin cmax])
    
    z=log10(s(h).num_eruption_events(ind_erupt_events));
    datatipRow = dataTipTextRow('C',z);
    hand(4).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_erupt_events));
    hand(4).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    
        z=log10(s(h).ave_freq(ind_erupt_events)*(365*24*3600*1e3));
    datatipRow = dataTipTextRow('C',z);
    hand(4).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_erupt_events));
    hand(4).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    
    
        z=log10(s(h).num_eruption_events(ind_grow_noerupt));
    datatipRow = dataTipTextRow('C',z);
    hand(3).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_grow_noerupt));
    hand(3).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    
    hand(3).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_grow_noerupt));
    hand(3).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    box on   
    hold on
    if save_figures ==1
        saveas(gca,[folder_name, 'er_and_grow_', composition, '_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6),'_maxn_500'],'epsc')
        saveas(gca,[folder_name, 'er_and_grow_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6),'_maxn_500'],'jpeg')  
    end
end
%% Eruptive lifespan
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    
    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);
    
    ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
     
    hand(1) = scatter(s(h).theta1(ind_no_erupt),...
        s(h).theta2(ind_no_erupt), pointsize,...
        log10(s(h).avg_vol_erupt(ind_no_erupt)./1e9), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit),...
        s(h).theta2(ind_erupt_nolimit), pointsize,...
        log10(s(h).eruptive_lifespan(ind_erupt_nolimit)./(3600*24*365)), 'o', 'filled', 'MarkerEdgeColor', 'k');
    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
    set(gca, 'XLim',[xlim_left,xlim_right])
    set(gca,'YLim',[ylim_bottom,ylim_top])
    
    title(['log Eruptive lifespan (years):', ' H_2O ',...
        num2str(s(h).H2O_content(1)*100),' wt%, ',...
        ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ',...
        num2str(s(h).P_lit(1)/1e6),' MPa'])
    
     colorbar    
     cmin = 2;
     cmax = 4;     
     colormap(jet);
     caxis([cmin cmax])
    
    if save_figures ==1
        saveas(gca,[folder_name, 'eruptive_lifespan_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'eruptive_lifespan_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
end
%% number of eruptions
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);    
    
        ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_no_erupt), s(h).theta2(ind_no_erupt),...
        pointsize, log10(s(h).num_eruption_events(ind_no_erupt)), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).num_eruption_events(ind_erupt_nolimit)), 'o', 'filled', 'MarkerEdgeColor', 'k');    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
    set(gca, 'XLim',[xlim_left,xlim_right])
    set(gca,'YLim',[ylim_bottom,ylim_top])
    title(['log Number of Eruptions:', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' MPa maxn 500'])    
    colorbar 
   % cmin = 0;
   % cmax = log10(1000); 
    colormap(jet);
   % caxis([cmin cmax])
    
    z=log10(s(h).num_eruption_events(ind_erupt_events));
    datatipRow = dataTipTextRow('C',z);
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_erupt_events));
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    
    box on   
    hold on
    if save_figures ==1
        saveas(gca,[folder_name, 'numer_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6),'_maxn_500'],'epsc')
        saveas(gca,[folder_name, 'numer_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6),'_maxn_500'],'jpeg')  
    end
end
%% Frequency of eruptions
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])    
    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);   
    
        
        ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_no_erupt),...
        s(h).theta2(ind_no_erupt), pointsize,...
        log10(s(h).ave_freq(ind_no_erupt)*(365*24*3600)), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit),...
        s(h).theta2(ind_erupt_nolimit), pointsize,...
        log10(s(h).ave_freq(ind_erupt_nolimit)*(365*24*3600*1e3)), 'o', 'filled', 'MarkerEdgeColor', 'k');    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    title(['log Avg Eruption Frequency (1/kyr):', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' MPa'])
    colorbar
    
    cmin = 0;
    cmax = 3;
    
    colormap(jet);
    caxis([cmin cmax])
    
    z=log10(s(h).ave_freq(ind_erupt_events)*(365*24*3600*1e3));
    datatipRow = dataTipTextRow('C',z);
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_erupt_events));
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    
    box on
    
    if save_figures ==1
        saveas(gca,[folder_name, 'freq_ave_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'freq_ave_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
end
%% Size of eruptions (volume)
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])    
 
    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);

        ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_no_erupt),...
        s(h).theta2(ind_no_erupt), pointsize,...
        log10(s(h).avg_vol_erupt(ind_no_erupt)./1e9), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit),...
        s(h).theta2(ind_erupt_nolimit), pointsize,...
        log10(s(h).avg_vol_erupt(ind_erupt_nolimit)./1e9), 'o', 'filled', 'MarkerEdgeColor', 'k');
    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    
    title(['log Avg Eruption Volume (km^3):', ' H_2O ',...
        num2str(s(h).H2O_content(1)*100),' wt%, ',...
        ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ',...
        num2str(s(h).P_lit(1)/1e6),' MPa'])
    
    colorbar
    
    cmin = -2.5;
    cmax = -0.5;
    
    colormap(jet);
    caxis([cmin cmax])
    
    z=log10(s(h).avg_vol_erupt(ind_erupt_events)./1e9);
    datatipRow = dataTipTextRow('C',z);
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_erupt_events));
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    box on
    
    if save_figures ==1
        saveas(gca,[folder_name, 'size_er_ave_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'size_er_ave_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
end
%% Mass Erupted:Intruded (intruded means chamber plus sum of all dikes,
%   erupted means sum of all erupted mass)
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])  

    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);
    
            ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    
    hand(1) = scatter(s(h).theta1(ind_no_erupt),...
        s(h).theta2(ind_no_erupt), pointsize,...
        s(h).mass_erupted(ind_no_erupt)./s(h).mass_intruded(ind_no_erupt),...
        'ko', 'filled', 'MarkerEdgeColor', 'k');
   hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit),...
        s(h).theta2(ind_erupt_nolimit), pointsize,...
        log10(s(h).mass_erupted(ind_erupt_nolimit)./s(h).mass_intruded(ind_erupt_nolimit)),...
        'o', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    
    title(['log Mass erupted/intruded:', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' MPa'])
    
    colorbar
    
    
    cmin = -4;
    cmax = 1;
    
    colormap(jet);
    caxis([cmin cmax])
    
    if save_figures ==1
        saveas(gca,[folder_name, 'M_er2in_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'M_er2in_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
end
%% Volume Erupted:Intruded (intruded means chamber plus sum of all dikes,
%   erupted means sum of all erupted mass)
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    
    ind_erupt_events = find(s(h).num_eruption_events>0);
    ind_no_erupt = find(s(h).num_eruption_events==0);
    
    ind_erupt_nolimit=[];
    for i = 1:length(ind_erupt_events)
        j = ind_erupt_events(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_no_erupt),...
        s(h).theta2(ind_no_erupt), pointsize,...
        s(h).vol_erupted(ind_no_erupt)./s(h).vol_intruded(ind_no_erupt),...
        'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit),...
        s(h).theta2(ind_erupt_nolimit), pointsize,...
        log10(s(h).vol_erupted(ind_erupt_nolimit)./s(h).vol_intruded(ind_erupt_nolimit)),...
        'o', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    
    title(['log Volume erupted/intruded:', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' MPa'])
    
    colorbar
    
    
    cmin = -3;
    cmax = 0;
    
    colormap(jet);
    caxis([cmin cmax])
    
    
    z=log10(s(h).vol_erupted(ind_erupt_events)./s(h).vol_intruded(ind_erupt_events));
    datatipRow = dataTipTextRow('C',z);
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_erupt_events));
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    
    box on
    
    if save_figures ==1
        saveas(gca,[folder_name, 'vol_er2in_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'vol_er2in_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
end
%% Mass Change
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    %subplot(3,3,h)
    ind_shrink = find(s(h).mass_change<1);
    ind_same = find(s(h).mass_change==1);
    ind_grow = find(s(h).mass_change>1);
    
    ind_erupt_nolimit=[];
    for i = 1:length(ind_grow)
        j = ind_grow(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_shrink), s(h).theta2(ind_shrink),...
        pointsize, s(h).mass_change(ind_shrink), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).mass_change(ind_erupt_nolimit)), 'o', 'filled', 'MarkerEdgeColor', 'k');    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',12);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])   
    title(['log Mass/Mass_{initial}:', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' kbar'])
    
   colormap('jet')
    colorbar
    
        z=log10(s(h).mass_change(ind_grow));
    datatipRow = dataTipTextRow('C',z);
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_grow));
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
    box on
    
  cmin = 0;
  cmax = 1;
  caxis([cmin cmax])
    

   
    
    box on
    if save_figures ==1
        saveas(gca,[folder_name, 'mass_growth_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'mass_growth_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% Mass Growth Rate
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    for ggg=1:900
        m_0(ggg)=s(h).tot_Mass{ggg}(1);
    end
    %subplot(3,3,h)
    ind_shrink = find(s(h).mass_change<1);
    ind_same = find(s(h).mass_change==1);
    ind_grow = find(s(h).mass_change>1);
    
        ind_erupt_nolimit=[];
    for i = 1:length(ind_grow)
        j = ind_grow(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_shrink), s(h).theta2(ind_shrink),...
        pointsize, log10(s(h).mass_growth_rate(ind_shrink).*(3600*24*365)./m_0(ind_shrink)'), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).mass_growth_rate(ind_erupt_nolimit).*(3600*24*365)./m_0(ind_erupt_nolimit)'), 'o', 'filled', 'MarkerEdgeColor', 'k');  
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',12);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])  
    title(['log Mass growth rate (kg/yr):', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' kbar'])
    
    colorbar
    
     cmin = -6;
     cmax = -3;
    caxis([cmin cmax])

    
    box on
    if save_figures ==1
        saveas(gca,[folder_name, 'mass_growth_rate_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'mass_growth_rate_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% Mass Growth Rate norm
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    for j=1:runs_per_folder
        longevity(j)=s(h).time{j}(end)./(60*60*24*365);   
    end
    %subplot(3,3,h)
      ind_shrink = find(s(h).mass_change<1);
    ind_same = find(s(h).mass_change==1);
    ind_grow = find(s(h).mass_change>1);
    
        ind_erupt_nolimit=[];
    for i = 1:length(ind_grow)
        j = ind_grow(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_shrink), s(h).theta2(ind_shrink),...
        pointsize, log10(s(h).mass_change(ind_shrink).*(3600*24*365)./longevity(ind_shrink)'), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).mass_change(ind_erupt_nolimit)./longevity(ind_erupt_nolimit)'), 'o', 'filled', 'MarkerEdgeColor', 'k');  
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',12);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])  
    title(['log Mass growth rate (kg/yr):', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' kbar'])
    
    colorbar
    
     cmin = -4;
     cmax = -3;
    caxis([cmin cmax])

    
    box on
    if save_figures ==1
        saveas(gca,[folder_name, 'mass_growth_rate_norm_', composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'mass_growth_rate_norm_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% Volume Change
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    
    
    ind_erupt_nolimit=[];
    for i = 1:length(s(h).num_eruption_events)
        if s(h).num_eruption_events(i)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit i];
        end
    end
    
    
    hand(1) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, (s(h).vol_change(ind_erupt_nolimit)), 'o', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',12);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    title(['V/V_{initial}:', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' kbar'])
    
    colorbar
    
    cmin = .9;
    cmax = 3;
    cwhite = 1;
    clength = (cmax-cmin);
    cnum = 5000;
    cbh = colorbar;
    cbh.Ticks = [cmin  1 cmax];
    blue = cnum*(cwhite-cmin)/clength;
    red = cnum*(cmax-cwhite)/clength;
    
    
    bmapx1 = linspace(0,0,blue/2)';
    bmapx2 = linspace(0,1,blue/2)';
    
    bmapy1 = linspace(0,0,blue/2)';
    bmapy2 = linspace(0,1,blue/2)';
    
    bmapz1 = linspace(0.5,1,blue/2)';
    bmapz2 = linspace(1,1,blue/2)';
    
    bmapx = [bmapx1;bmapx2];
    bmapy = [bmapy1;bmapy2];
    bmapz = [bmapz1;bmapz2];
    
    rmapx1 = linspace(1,1,red/2)';
    rmapx2 = linspace(1,0.5, red/2)';
    
    rmapy1 = linspace(1,0,red/2)';
    rmapy2 = linspace(0,0,red/2)';
    
    rmapz1 = linspace(1,0,red/2)';
    rmapz2 = linspace(0,0,red/2)';
    
    rmapy = [rmapy1;rmapy2];
    rmapx = [rmapx1;rmapx2];
    rmapz = [rmapz1; rmapz2];
    
    rx = linspace(1,1,red)';
    ry1 = linspace(1,1,red/12)';
    ry2 = linspace(1,0,11*red/12)';
    rz1 = linspace(1,0,red/12)';
    rz2 = linspace(0,0,11*red/12)';
    
    mapx = [bmapx; 1; rmapx];
    mapy = [bmapy;1;rmapy];
    mapz = [bmapz; 1;rmapz];
    
    newMap = [mapx, mapy, mapz];
    
    
    colormap(newMap);
    caxis([cmin cmax])
    
    box on
    if save_figures ==1
        saveas(gca,[folder_name, 'v_change_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'v_change_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% Volume Growth Rate
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
       ind_shrink = find(s(h).vol_change<1);
    ind_same = find(s(h).vol_change==1);
    ind_grow = find(s(h).vol_change>1);
    
            ind_erupt_nolimit=[];
    for i = 1:length(ind_grow)
        j = ind_grow(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_shrink), s(h).theta2(ind_shrink),...
        pointsize, s(h).vol_growth_rate(ind_shrink).*(3600*24*365), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).vol_growth_rate(ind_erupt_nolimit).*(3600*24*365)./1e9), 'o', 'filled', 'MarkerEdgeColor', 'k');  
     xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',12);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    title(['log V growth rate (km^3/yr):', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' kbar'])
    
    colorbar
    colormap('jet')
    
     cmin = -4;
     cmax = -2;
          caxis([cmin cmax]) 
    box on
    if save_figures ==1
        saveas(gca,[folder_name, 'v_growth_rate_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'v_growth_rate_',composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% longevity
for h=1:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    
    for j=1:runs_per_folder
        longevity(j)=s(h).time{j}(end)./(60*60*24*365);   
    end
   
    ind_erupt_nolimit=[];
    for i = 1:length(s(h).num_eruption_events)
        if s(h).num_eruption_events(i)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit i];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(longevity(ind_erupt_nolimit)), 'o', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    
    title(['Longevity (log(yrs)):', ' H_2O ', num2str(s(h).H2O_content(1)*100),...
        ' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ',...
        num2str(s(h).P_lit(1)/1e6),' MPa'])
    
    colorbar
    
    cmin =3;
    cmax = 5;
    
    colormap(jet);
    caxis([cmin cmax])
    
    z=log10(longevity);
    datatipRow = dataTipTextRow('C',z);
    hand(1).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run);
    hand(1).DataTipTemplate.DataTipRows(end+1) = datatipRow2;

        box on
    if save_figures ==1
        saveas(gca,[folder_name, 'longevity_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'longevity_', composition,'_H2O_',num2str(s(h).H2O_content(1)*100),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% longevity Comparison
for h=3:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    
    for j=1:runs_per_folder
        longevity_1(j)=s(1).time{j}(end)./(60*60*24*365);   
        longevity(j)=s(h).time{j}(end)./(60*60*24*365);   
    end
   

    hand(1) = scatter(s(h).theta1(), s(h).theta2(),...
        pointsize, (longevity./longevity_1), 'o', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',ftsize);
set(gca, 'XLim',[xlim_left,xlim_right])
set(gca,'YLim',[ylim_bottom,ylim_top])
    
    title(['Longevity normalized by longevity in low water low co2 case :', ' H_2O ', num2str(s(h).H2O_content(1)*100),...
        ' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ',...
        num2str(s(h).P_lit(1)/1e6),' MPa'])
    
    colorbar
    
    if h==3
        cmin = 0.9;
        cmax = 1.1;
    end
    
    if h==4
        cmin = 2;
        cmax = 15;
    end
    
    colormap(jet);
    caxis([cmin cmax])
    
    z=longevity./longevity_1;
    datatipRow = dataTipTextRow('C',z);
    hand(1).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run);
    hand(1).DataTipTemplate.DataTipRows(end+1) = datatipRow2;

        box on
    if save_figures ==1
        saveas(gca,[folder_name, 'longevitycomparison_', composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'longevitycomparison_', composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end
%% Mass Change Comparison
for h=3:conditions
    figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
    %subplot(3,3,h)
    ind_shrink = find(s(h).mass_change<1);
    ind_same = find(s(h).mass_change==1);
    ind_grow = find(s(h).mass_change>1);
    
    ind_erupt_nolimit=[];
    for i = 1:length(ind_grow)
        j = ind_grow(i);
        if s(h).num_eruption_events(j)~=999
            ind_erupt_nolimit = [ind_erupt_nolimit j];
        end
    end
    
    hand(1) = scatter(s(h).theta1(ind_shrink), s(h).theta2(ind_shrink),...
        pointsize, s(h).mass_change(ind_shrink), 'ko', 'filled', 'MarkerEdgeColor', 'k');
    hold on
    hand(2) = scatter(s(h).theta1(ind_erupt_nolimit), s(h).theta2(ind_erupt_nolimit),...
        pointsize, log10(s(h).mass_change(ind_grow)./s(1).mass_change(ind_grow)), 'o', 'filled', 'MarkerEdgeColor', 'k');    
    xlabel('\Theta_{1} = \tau_{cool} / \tau_{in}'), ylabel('\Theta_{2} = \tau_{relax} / \tau_{in}')
    axis('square')
    set(gca,'xscale','log', 'yscale', 'log')
    set(gca,'FontSize',12);
    set(gca, 'XLim',[xlim_left,xlim_right])
    set(gca,'YLim',[ylim_bottom,ylim_top])   
    title(['Mass Change compare log:', ' H_2O ', num2str(s(h).H2O_content(1)*100),' wt%, ', ' CO_2 ', num2str(s(h).CO2_content(1)*1e6) ,' ppm ', num2str(s(h).P_lit(1)/1e6),' kbar'])
    
   colormap('jet')
    colorbar
    if h==3  
  cmin = -0.1;
  cmax = 0.1;
  end
  
  if h==4  
  cmin = 0.1;
  cmax = 2;
  end
  caxis([cmin cmax])
        z=log10(s(h).mass_change(ind_grow)./s(1).mass_change(ind_grow));
    datatipRow = dataTipTextRow('C',z);
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow;
        datatipRow2 = dataTipTextRow('run',s(h).run(ind_grow));
    hand(2).DataTipTemplate.DataTipRows(end+1) = datatipRow2;
  
  
    

   
    
    box on
    if save_figures ==1
        saveas(gca,[folder_name, 'mass_growth_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'epsc')
        saveas(gca,[folder_name, 'mass_growth_',composition,'_H2O_',num2str(s(h).H2O_content(1)*1000),'_CO2_', num2str(s(h).CO2_content(1)*1e6)],'jpeg')
    end
    
end

