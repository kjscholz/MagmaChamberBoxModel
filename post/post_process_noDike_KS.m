clear all
tic
sf1 = 'output/';
% sf2 = ["silicic_H2O_3_CO2_50/";...
%        "silicic_H2O_3_CO2_500/";...
%        "silicic_H2O_3_CO2_1000/";...
%        "silicic_H2O_45_CO2_50/";...
%        "silicic_H2O_45_CO2_500/";...
%        "silicic_H2O_45_CO2_1000/";...
%        "silicic_H2O_6_CO2_50/";...
%        "silicic_H2O_6_CO2_500/";...
%        "silicic_H2O_6_CO2_1000/"];
%
sf2 = ["mafic_H2O_0.5_CO2_500/";
    "mafic_H2O_0.5_CO2_5000/";
    "mafic_H2O_0.5_CO2_10000/";
    "mafic_H2O_2_CO2_500/"];
% 
% sf2 = ["silicic_H2O_4_CO2_100/";
 %    "silicic_H2O_4_CO2_1000/"];
 %        "silicic_H2O_4_CO2_1000/"];
  %      "silicic_H2O_6_CO2_100/";
  %      "silicic_H2O_6_CO2_1000/"];

runs_per_folder =900;
N = runs_per_folder; %runs per folder
M = length(sf2);
for q=1:M
   %for u = 1343
     for u=1:runs_per_folder
        if q==5 %&& u>440
            run_num=u+1000;
        else
            run_num=u+1000;
        end
        FN((q-1)*runs_per_folder+u,:)=string([sf1 char(sf2(q)) 'run_' num2str(run_num) '.mat']);
    end
end


h=0;
for f=1:M
    h=h+1;
    j = 0;
    
 %  for iii = 1343 
    for iii = 1:runs_per_folder
        choices=[271,333,686,511,799];
        %iii=choices(hhh);
        load(FN((f-1)*runs_per_folder+iii))
        j = j+1;
        run_n      = iii;
        time_dis = []; dur_dis  = []; mass_dis = []; vol_dis  = [];
        
        % get input
        s(h).run(j,1)=iii;
        s(h).H2O_content(j,1)=range_water;
        s(h).CO2_content(j,1)=range_co2;
        s(h).P_lit(j,1) = P_lit;
        s(h).T_0(j,1)  = T_0;
        s(h).eps_g0(j,1) = eps_g(1);
        s(h).eps_gf(j,1) = eps_g(end);
        s(h).eps_xf(j,1) = eps_x(end);
        s(h).R(j,1) = (V_0*(3/4)*(1/pi))^(1/3);
        s(h).rho_m0(j,1) = rho_m(1);
        s(h).rho_x0(j,1) = rho_x(1);
        s(h).mdot_in(j,1) = mdot_in;
        s(h).mdot_out(j,1) = mdot_out;
        s(h).rho_0(j,1)  = rho(1);
        s(h).beta_r(j,1) = param.beta_r;
        s(h).beta_m(j,1) = param.beta_m;
        s(h).beta_x(j,1) = param.beta_x;
        s(h).tau_inj(j,1) = tau_inj;
        s(h).tau_cooling(j,1) = tau_cooling;
        s(h).tau_visco_elastic(j,1) = tau_visco_elastic;
        s(h).Tb(j,1) = Tb;
        s(h).P_crit(j,1) = P_crit;
        s(h).rho_r(j,1) = param.rho_r;
        s(h).m_eq_in(j,1) = m_eq_in;
      %  s(h).eta_r0(j,1) = param.eta_r0;
        s(h).X_co20(j,1) = X_co2(1);
        [~,~,~,~,s(h).C_co20(j,1),~,~,~] = exsolve_mafic(s(h).P_lit(j,1),s(h).T_0(j,1),s(h).X_co20(j,1));
     %   s(h).eta(j,1)    =   param.eta;       % Magma viscosity (Pa s)
    %    s(h).mu(j,1)     = param.mu;       % Shear modulus crust (Pa)
   %     s(h).pr(j,1) = param.pr;      % Poisson's ratio
       % s(h).r(j,1) = param.r;         % radius of conduit connecting dike and chamber (m)
       % s(h).L(j,1) =    param.L;         % length of conduit connecting dike and chamber (m)
       % s(h).w(j,1) =   param.w;       % fissure length as fraction of dike half-length
       % s(h).con(j,1) = param.con; % hydraulic connectivity between chamber and dike
        
        s(h).time_of_run(j,1)= time_of_run;
        
        s(h).eps_x{j}=eps_x;
        
        s(h).time{j}        = time;
        s(h).Pc{j}          = P;
        s(h).T{j}           = T;
        s(h).eps_g{j}       = eps_g;
        s(h).V{j}           = V;
        s(h).rho_m{j}       = rho_m;
        s(h).rho_x{j}       = rho_x;
      %  s(h).Pd{j}          = Pd;
      %  s(h).a{j}           = a;
      %  s(h).b{j}           = b;
        
        
        
        s(h).m_eq{j}        = m_eq;
        s(h).rho_g{j}       = rho_g;
        s(h).drho_g_dP{j}   = drho_g_dP;
        s(h).rho{j}         = rho;
        s(h).beta{j}        = ((1./rho).*((1-eps_x-eps_g).*rho_m./param.beta_m...
            +eps_x.*rho_x./param.beta_x+eps_g.*drho_g_dP)).^(-1);
        s(h).beta_0(j,1)    = s(h).beta{j}(1);
        s(h).c{j}           = c;
        s(h).X_co2{j}       = X_co2;
        s(h).tot_Mass{j}    = tot_Mass;
        s(h).tot_Mass_H2O{j}= tot_Mass_H2O;
        s(h).tot_Mass_CO2{j}= tot_Mass_CO2;
        
        s(h).time_erupt{j} = time_erupt;
        
        
      %  s(h).vol_erupt{j} = vol_erupt;
       % s(h).mass_erupt{j} = mass_erupt;
      %  s(h).vol_dike{j} = vol_dike;
      %  s(h).mass_dike{j} = mass_dike;
     %   s(h).dur_erupt{j} = dur_erupt;
      %  s(h).num_dike_events(j,1) = length(time_erupt);
        s(h).num_eruption_events(j,1) = length(time_erupt);
    %    s(h).avg_mass_dike(j,1)=mean(mass_dike);
      s(h).mass_growth_rate(j,1) = (tot_Mass(end)-tot_Mass(1))/time(end); % kg/s
        s(h).vol_growth_rate(j,1) = (V(end)-V(1))/time(end); % m^3 / s
   
    if ~isempty(time_erupt)
        s(h).eruptive_lifespan(j,1) = time_erupt(end)-time_erupt(1); % seconds

        s(h).avg_mass_erupt(j,1)=mean(mass_erupt);
        s(h).avg_vol_erupt(j,1)=mean(vol_erupt);
     %   s(h).avg_vol_dike(j,1)=mean(vol_dike);  
   %     s(h).avg_mass_dis(j,1) = mean(mass_dike)+mean(mass_erupt);
     %   s(h).avg_vol_dis(j,1) = mean(vol_dike)+mean(vol_erupt);  
       s(h).mass_erupted(j,1) = sum(mass_erupt);
       s(h).vol_erupted(j,1) = sum(vol_erupt);
       
    else
        
        s(h).eruptive_lifespan(j,1) = 0; % seconds
        s(h).avg_mass_erupt(j,1)=0;
        s(h).avg_vol_erupt(j,1)=0;
       s(h).mass_erupted(j,1) = 0;
       s(h).vol_erupted(j,1) = 0;
    end
       
        s(h).mass_intruded(j,1) = tot_Mass(end);
     
        s(h).vol_intruded(j,1) = V(end);
        
        s(h).repose_times{j} = diff(time_erupt);
        s(h).ave_freq(j,1) = 1/mean(diff(time_erupt)); % 1/second
        
        
     
                er_P=[];er_ind=[]; end_er_tf=[]; er_end=[]; er_end_ind=[];
                ind_array=1:1:length(P);
                [er_P_prelim, er_ind_prelim]=findpeaks(s(h).Pc{j});
                end_er_tf=islocalmin(s(h).Pc{j});
                [er_end_prelim,er_end_ind_prelim]=findpeaks(double(end_er_tf));
        
                gg=1;
                for kk=1:length(er_P_prelim)
                    if er_P_prelim(kk)>(P_lit+P_crit-1000)
                        er_P(gg,1)=er_P_prelim(kk);
                        er_ind(gg,1)=er_ind_prelim(kk);
                        gg=gg+1;
                    end
                end
                ll=1;
                for ss=1:length(er_end_prelim)
                    P(er_end_ind_prelim(ss));
                    if abs(P(er_end_ind_prelim(ss))-(P_lit))<1000
                        er_end_ind(ll,1)=er_end_ind_prelim(ss);
                        ll=ll+1;
                    end
        
                end
        
        
                if ~isempty(er_end_ind)&& ~isempty(er_ind)
        
                    if er_end_ind(1)<er_ind(1)
                        er_end_ind=er_end_ind(2:end);
                    end
                    if length(er_end_ind) < length(er_ind)
                        er_end_ind=[er_end_ind; ind_array(end)];
                    end
                end
        
                if isempty(er_end_ind) && ~isempty(er_ind)
                        er_end_ind=[ind_array(end)];
                end
        
        
        
                time_dis=[];dur_dis=[];mass_dis=[]; vol_dis=[];
                time_dis=s(h).time{j}(er_ind);
                time_dis_end=s(h).time{j}(er_end_ind);
                dur_dis=s(h).time{j}(er_end_ind)-s(h).time{j}(er_ind);
                vol_dis=s(h).V{j}(er_ind)-s(h).V{j}(er_end_ind);
                mass_dis=s(h).tot_Mass{j}(er_ind)-s(h).tot_Mass{j}(er_end_ind);
                rho_dis = s(h).rho{j}(er_ind);
                vol_dis2=mass_dis./rho_dis;
                 s(h).index_erupt_begin{j}=er_ind;
                 s(h).index_erupt_end{j}=er_end_ind;
                 s(h).M_dis{j}                 = mass_dis;
                 s(h).V_dis{j}                 = vol_dis2;
                 
                 s(h).mass_erupted2(j,1) = sum(mass_dis);
                 s(h).vol_erupted2(j,1)  = sum(vol_dis2);
                 
                 s(h).avg_mass_erupt2(j,1) = mean(mass_dis);
                 s(h).avg_vol_erupt2(j,1) = mean(vol_dis2);
                 
                 s(h).dur_dis{j}               = dur_dis;
                 s(h).tot_Mass_begin_dis{j}    = tot_Mass(er_ind);
                 s(h).V_begin_dis{j}           = V(er_ind);
                 s(h).beta_begin_dis{j}        = s(h).beta{j}(er_ind);
                 s(h).time_begin_dis{j}        = s(h).time{j}(er_ind);
                 s(h).time_end_dis{j}          = s(h).time{j}(er_end_ind);
                 s(h).rho_begin_dis{j}         = s(h).rho{j}(er_ind);
                 s(h).eps_g_begin_dis{j}       = s(h).eps_g{j}(er_ind);
                 s(h).V_ch_scaled_begin_dis{j} = s(h).V_begin_dis{j}.*(s(h).beta_begin_dis{j}).^(-1).*(s(h).P_crit(j,1));
      
                if ~isempty(er_ind)
                    [Vd,Ver,W_D]=dikeVolume(s(h).V_ch_scaled_begin_dis{j}./1e9,-s(h).V_dis{j}./1e9,s(h).eps_g_begin_dis{j});
                    s(h).V_d_dis_end{j} = Vd;
                    s(h).V_er_dis_end{j}= Ver;
                else
                    s(h).V_d_dis_end{j} = [];
                    s(h).V_er_dis_end{j}= [];
                end
        
        
        s(h).mass_change(j,1) = tot_Mass(end)./tot_Mass(1);
        s(h).eps_g_change(j,1)=eps_g(end)./eps_g(1);
        s(h).eps_g_initial(j,1)=eps_g(1);
        s(h).eps_g_final(j,1)=eps_g(end);
        s(h).vol_change(j,1) = V(end)/V(1);
        %         % Number of eruptions
        %         s(h).num_dis(j,1) = length(time_dis);
        %         s(h).num_er(j,1)  = length(nonzeros(s(h).V_er_dis_end{j}));
        %
        
        %     % Number of phases (2, 3, or 23 = transition)
        if isempty(find(s(h).eps_g{j} < 1e-10))
            s(h).phase(j,1) = 3;
        elseif isempty(find(s(h).eps_g{j} > 1e-5))
            s(h).phase(j,1) = 2;
        else
            s(h).phase(j,1) = 23;
        end
        
        %     % Viscosity of host rock
        s(h).eta_r_all{j}    = eta_r_vec; % crustal viscosity
        s(h).eta_time_all{j} = eta_t_vec; % times at which viscosity is evaluated
        
        %     % Timescale for growth of chamber
        s(h).growrate_V(j,1) = (V(end)-V(1))/time(end);
        s(h).growrate_M(j,1) = (tot_Mass(end)-tot_Mass(1))/time(end);
        
        %     % erupted:added
        
        %
        %         s(h).avg_dis_M(j,1) = mean(mass_dis);
        %         s(h).avg_dis_V(j,1) = mean(vol_dis);
        %         s(h).dis2add_M(j,1) = -1*sum(mass_dis)/(mdot_in*time(end));
        %         s(h).dis2add_V(j,1) = -1*sum(vol_dis)/(mdot_in*time(end)./rho_0);
        %
        %         s(h).ext2int_V(j,1) = sum(s(h).V_er_dis_end{j}.*1e9)/(sum(s(h).V_d_dis_end{j}.*1e9)+s(h).V{j}(end));
        %         %s(h).ext2int_M(j,1) = sum(s(h).V_er_dis_end{j}.*1e9.*s(h).rho_begin_dis{j})/(sum((s(h).V_d_dis_end{j}.*1e9.*s(h).rho_begin_dis{j}))+(s(h)V{j}(end).*rho(end)));
        %         s(h).ext2dike_M(j,1) = sum(s(h).V_er_dis_end{j}.*1e9.*s(h).rho_begin_dis{j})/(sum(s(h).V_d_dis_end{j}.*1e9.*s(h).rho_begin_dis{j}));
        %         s(h).ext2dike_V(j,1) = sum(s(h).V_er_dis_end{j}.*1e9)/(sum(s(h).V_d_dis_end{j}.*1e9)+s(h).V{j}(end));
        %
        %         % repose times and eruption frequency
        %         if length(time_dis)<2
        %             s(h).repose_times{j} = NaN;
        %         else
        %             for q = 1:length(time_dis)-1
        %                 s(h).repose_times{j}(q) = time_dis(q+1)-time_dis(q);
        %             end
        %         end
        %
        %         s(h).repose_ave(j,1) = mean(s(h).repose_times{j}); % mean repose time
        %         s(h).freq_ave(j,1) = 1/s(h).repose_ave(j); % frequency in seconds^-1
        %
        s(h).tau_inj_vec{j}=(rho.*V)./mdot_in;
        s(h).tau_cooling_vec{j}=(((3/(4*pi)).*V).^(1/3)).^2./param.kappa;
        s(h).tau_visco_elastic_vec{j}=eta_r_vec/P_crit;
        % Thetas
        s(h).theta1(j,1) = s(h).tau_cooling(j)/s(h).tau_inj(j);
        s(h).theta2(j,1) = s(h).tau_visco_elastic(j)/s(h).tau_inj(j);
        
        s(h).theta1_final(j,1) = s(h).tau_cooling_vec{j}(end)/s(h).tau_inj_vec{j}(end);
        s(h).theta2_final(j,1) = s(h).tau_visco_elastic_vec{j}(end)/s(h).tau_inj_vec{j}(end);
        
        
    end
    
    
end
%clearvars -except s
toc