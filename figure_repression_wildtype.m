clear
close all

addpath(genpath('./lib'))


%% create color map
%Green color map
cmap_green = [[247,252,253];[229,245,249];[204,236,230];[153,216,201];[102,194,164];[65,174,118];[35,139,69];[0,109,44];[0,68,27]];

old_stepNum = size(cmap_green,1);
new_stepNum = 256;

cmap_green_1 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,1)),linspace(0,1,new_stepNum));
cmap_green_2 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,2)),linspace(0,1,new_stepNum));
cmap_green_3 = interp1(linspace(0,1,old_stepNum),double(cmap_green(:,3)),linspace(0,1,new_stepNum));
cmap_green_new = [cmap_green_1' cmap_green_2' cmap_green_3']/256;


%Red color map
cmap_red = [[255,255,229];[255,247,188];[254,227,145];[254,196,79];[254,153,41];[236,112,20];[204,76,2];[153,52,4];[102,37,6]];

old_stepNum = size(cmap_red,1);
new_stepNum = 256;

cmap_red_1 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,1)),linspace(0,1,new_stepNum));
cmap_red_2 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,2)),linspace(0,1,new_stepNum));
cmap_red_3 = interp1(linspace(0,1,old_stepNum),double(cmap_red(:,3)),linspace(0,1,new_stepNum));

cmap_red_new = [cmap_red_1' cmap_red_2' cmap_red_3']/256;

%Magenta color map
cmap_magenta = [217 82 147];
cmap_white = [255 255 255];

cmap_magenta_1 = interp1([0 1],[double(cmap_white(1)) double(cmap_magenta(1))],linspace(0,1,new_stepNum));
cmap_magenta_2 = interp1([0 1],[double(cmap_white(2)) double(cmap_magenta(2))],linspace(0,1,new_stepNum));
cmap_magenta_3 = interp1([0 1],[double(cmap_white(3)) double(cmap_magenta(3))],linspace(0,1,new_stepNum));

cmap_magenta_new = [cmap_magenta_1' cmap_magenta_2' cmap_magenta_3']/256;

beta = -0.5;
cmap_magenta_updated = brighten(cmap_magenta_new,beta);

%Color for lines
color_green = [38 143 75]/256;
green = [122 168 116]/256;
mRNA_red = brighten([212 100 39]/256,.2);
magenta = [217 82 147]/256;


%% parameters
knirps_offset = 375698.13;
ap_lim = 0.02;


%% initialization

projectName = 'optokni_eve4+6_WT'; 

% load data
load(['data' filesep 'main_analysis' filesep projectName filesep 'spot_struct.mat'])

ever_on_vec = zeros(size(spot_struct));
% initialize longform vectors for regression
ap_vec_long = [];
time_vec_long = [];
knirps_vec_long_raw = [];
fluo_raw_long = [];
fluo_zeros_long = [];
mRNA_vec_long = [];

post_turn_on_flags = [];
post_turn_off_flags = [];
ever_on_flags = [];

for i = 1:length(spot_struct)
  
    % extract core vectors 
    fluo_vec = spot_struct(i).fluo;
    time_vec = spot_struct(i).time;
    knirps_vec = spot_struct(i).rawNCProtein;
    ap_vec = spot_struct(i).APPosNucleus;
    
    if time_vec(end) - time_vec(1) >=30*60
        
        % make average vectors        
      
        % get off and on indices
        ever_on_vec(i) = any(~isnan(fluo_vec));
        
        post_on_vec = zeros(size(ap_vec));
        post_off_vec = zeros(size(ap_vec));
        if ever_on_vec(i)
            start_i = find(~isnan(fluo_vec),1);
            stop_i = find(~isnan(fluo_vec),1,'last');
            if true%stop_i < length(fluo_vec)-10
                post_off_vec(stop_i+1:end) = 1;
            end
            
            if start_i > 1
                post_on_vec(start_i+1:end) = 1;
            end
        end
        
        % make regression vectors
        post_turn_on_flags = [post_turn_on_flags post_on_vec];
        post_turn_off_flags = [post_turn_off_flags post_off_vec];
        
        fluo_zeros = fluo_vec;
        all_zeros = fluo_vec;
        fluo_zeros(post_on_vec&~post_off_vec&isnan(fluo_vec)) = 0;
        
        all_zeros(isnan(all_zeros)) = 0;
        mRNA_vec_long = [mRNA_vec_long all_zeros];
        
        mean_fluo(i) = nanmean(fluo_zeros);
        
        knirps_vec_long_raw = [knirps_vec_long_raw knirps_vec];
        ap_vec_long = [ap_vec_long ap_vec];
        time_vec_long = [time_vec_long time_vec];
        
    end
end


%% calculate long vectors

%set parameters
timeBins = 61;

ap_bins = linspace(-0.2,0.2,31);
ap_bins_plot = (ap_bins(1:end-1) + ap_bins(2:end))/2;
ap_axis = 100*(ap_bins(1:end-1) + diff(ap_bins)/2);
time_bins = linspace(0,35*60,timeBins);
time_bins_plot = (time_bins(1:end-1)+time_bins(2:end))/120;

% calculate mean vectors
knirps_vec_long = knirps_vec_long_raw - knirps_offset;

ap_groups = discretize(ap_vec_long,ap_bins); 
time_groups = discretize(time_vec_long,time_bins); 

frac_inst_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_mean = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_inst_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_std = NaN(length(time_bins)-1,length(ap_bins)-1);

frac_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_inst_on_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_num = NaN(length(time_bins)-1,length(ap_bins)-1);


frac_on_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
frac_on_knirps_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
eve_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);
knirps_time_array_ste = NaN(length(time_bins)-1,length(ap_bins)-1);

for t = 1:length(time_bins)
    for a = 1:length(ap_bins)
        time_window_filter = time_groups==t & ap_groups==a & post_turn_on_flags;
        if sum(time_window_filter) > 10
            still_on_flags_time = ~post_turn_off_flags(time_window_filter);
            inst_on_flags_time = mRNA_vec_long(time_window_filter)>0;

            frac_on_time_array_mean(t,a) = nanmean(still_on_flags_time);
            frac_on_time_array_num(t,a) = sum(~isnan(still_on_flags_time));
            frac_inst_on_time_array_mean(t,a) = mean(inst_on_flags_time);
            frac_inst_on_time_array_num(t,a) = length(inst_on_flags_time);
            
            eve_sample = mRNA_vec_long(time_window_filter);
            eve_sample(isnan(eve_sample)) = 0;
            eve_time_array_mean(t,a) = mean(eve_sample);
            eve_time_array_num(t,a) = sum(~isnan(eve_sample));
            eve_time_array_std(t,a) = std(eve_sample);
            eve_time_array_ste(t,a) = eve_time_array_std(t,a)/sqrt(eve_time_array_num(t,a));
            
            knirps_sample = knirps_vec_long(time_window_filter);
            knirps_sample(isnan(knirps_sample)) = 0;
            knirps_time_array_mean(t,a) = mean(knirps_sample);
            knirps_time_array_num(t,a) = sum(~isnan(knirps_sample));
            knirps_time_array_std(t,a) = std(knirps_sample);
            knirps_time_array_ste(t,a) = knirps_time_array_std(t,a)/sqrt(knirps_time_array_num(t,a));
            
        end

    end
end
knirps_time_array_mean(isnan(knirps_time_array_mean)) = 0;
eve_time_array_mean(isnan(eve_time_array_mean)) = 0;
frac_inst_on_time_array_mean(isnan(frac_inst_on_time_array_mean)) = 0;

knirps_time_array_mean = movmean(knirps_time_array_mean,5,1);
eve_time_array_mean = movmean(eve_time_array_mean,5,1);
frac_inst_on_time_array_mean = movmean(frac_inst_on_time_array_mean,5,1);

%% plot calculated mean vector

WT.knirps_mean = knirps_time_array_mean;
WT.eve_mean = eve_time_array_mean;
WT.inst_on = frac_inst_on_time_array_mean;


% plot for wildtype
WT_mean_kni_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.knirps_mean)
colorbar
colormap(cmap_green_new)
caxis([0 12E5])
xlim([-0.18 0.18])
ylim([2.5 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')
pbaspect([3 2 1])


WT_mean_eve_fig = figure;
imagesc(ap_bins_plot,time_bins_plot,WT.eve_mean)
colorbar
colormap(cmap_magenta_updated)
%caxis([0 2E5])
xlim([-0.18 0.18])
ylim([2.5 35])
xlabel('AP position (% embryo length)')
ylabel('time (min)')
pbaspect([3 2 1])


%% Find Knirps center

for i = 1:length(time_bins_plot)

    a = ap_bins_plot';
    b = WT.knirps_mean(i,:)';
    b(b<0) = 0;
    
    
    % define initial conditions and bounds
    init_vec = [7e5 0 0.05];
    ub_vec = [1e8 ap_bins_plot(end) 1];
    lb_vec = [1e3 ap_bins_plot(1) 0];
    
    options = optimoptions(@lsqnonlin,'Display','off');
    
    % define function
    gauss_fun = @(x) x(1)*exp(-0.5*((x(2)-a)./x(3)).^2);
    ob_fun = @(x) gauss_fun(x)-b;
        
    % perform fit
    fitted_param = lsqnonlin(ob_fun,init_vec,lb_vec,ub_vec,options);
    
    % plot fitted center
    ap_plot = linspace(ap_bins_plot(1),ap_bins_plot(end),101);

    center_pos(i) = fitted_param(2);
    

    %fig = figure;
    %hold on
    %plot(a,b)
    %plot(ap_plot,fitted_param(1)*exp(-0.5*((fitted_param(2)-ap_plot)./fitted_param(3)).^2))

end


fig = figure;
%plot(time_bins_plot,center_pos,'- .','MarkerSize',10,'LineWidth',2) 
plot(time_bins_plot,center_pos)
ylim([-0.05 0.05])
xlabel('time (min)')
ylabel('Knirps center')


%{
center_pos = zeros(length(time_bins_plot),1);

for i = 1:length(time_bins_plot)

    x = WT.knirps_mean(i,5:30);

    [pks,locs] = findpeaks(x);
    locs = locs((locs>=10) & (locs<=15));
    center_pos(i) = locs;

end

fig = figure;
plot(center_pos) 
ylim([0 20])
%}
    

%% Figure: plot mean fluorescence vs time (not aligned)
 
nBoots = 100;
 
ever_on_vec = [];
mean_ap = [];
time_orig_long = [];
knirps_orig_long = [];
fluo_orig_long = [];
off_time_long = [];

 
fig  = figure;
hold on
 
for i = 1:length(spot_struct)
    
    if (spot_struct(i).TraceQCFlag == 1)
        
        % extract core vectors 
        fluo_vec_orig = spot_struct(i).fluo;
        time_vec_orig = spot_struct(i).time;
        knirps_vec_orig = spot_struct(i).rawNCProtein-knirps_offset;
        ap_vec_orig = spot_struct(i).APPosNucleus;
        
        % get off and on indices
        ever_on_orig = any(~isnan(fluo_vec_orig));
        mean_ap_orig = nanmean(ap_vec_orig);
        
        mean_ap = [mean_ap nanmean(ap_vec)];
        
        if ever_on_orig
            start_i = find(~isnan(fluo_vec_orig),1);
            stop_i = find(~isnan(fluo_vec_orig),1,'last');
            off_time_orig = time_vec_orig(stop_i);
            off_knirps_vec_orig = knirps_vec_orig(stop_i);
            off_spot_fluo_orig = fluo_vec_orig(stop_i);
            off_ap_orig = ap_vec_orig(stop_i);
        end
        
        if (mean_ap_orig > -ap_lim) && (mean_ap_orig < ap_lim)
           time_orig_long = [time_orig_long time_vec_orig];
           knirps_orig_long = [knirps_orig_long knirps_vec_orig];
           fluo_orig_long = [fluo_orig_long fluo_vec_orig];
           off_time_long = [off_time_long off_time_orig];
           
           %plot((time_vec_orig)/60,fluo_vec_orig,'Color', [175 175 175]/255);
        end
        
    end
    
end

fluo_orig_long(isnan(fluo_orig_long)) = 0;

time_bin = linspace(0,50,51);
time_plot = (time_bin(1:end-1) + time_bin(2:end))/2;
time_groups = discretize(time_orig_long/60,time_bin);
 
knirps_vec_mean = zeros(length(time_bin)-1,1);
knirps_vec_ste = zeros(length(time_bin)-1,1);
fluo_vec_mean = zeros(length(time_bin)-1,1);
fluo_vec_ste = zeros(length(time_bin)-1,1);
 
for i = 1:length(time_plot)-1
 
    time_filter_long = time_groups==i;
    
    if sum(time_filter_long) > 10   
        boot_samples_fluo = bootstrp(nBoots,@nanmean,fluo_orig_long(time_filter_long));
        fluo_vec_mean(i) = nanmean(boot_samples_fluo);
        fluo_vec_ste(i) = std(boot_samples_fluo);
 
        boot_samples_knirps = bootstrp(nBoots,@nanmean,knirps_orig_long(time_filter_long));
        knirps_vec_mean(i) = nanmean(boot_samples_knirps);
        knirps_vec_ste(i) = std(boot_samples_knirps);
    end
    
end  

fluo_vec_mean = movmean(fluo_vec_mean,3,1);
%knirps_vec_mean = movmean(knirps_vec_mean,3,1);


hold on
yyaxis left

%plot(time_plot,knirps_vec_mean,'LineWidth',5)
errorbar(time_plot,knirps_vec_mean,knirps_vec_ste,'Color','k','CapSize',0);
plot(time_plot,knirps_vec_mean,'-k','LineWidth',1)
scatter(time_plot,knirps_vec_mean,50,'MarkerFaceColor',green,'MarkerEdgeColor','k')
ylim([1.25E5 11.25E5])
yyaxis right
errorbar(time_plot,fluo_vec_mean,fluo_vec_ste,'Color','k','CapSize',0);
plot(time_plot,fluo_vec_mean,'-k','LineWidth',1)
scatter(time_plot,fluo_vec_mean,50,'MarkerFaceColor',magenta,'MarkerEdgeColor','k')

xlim([5 30])
ylim([-5E3 2E5])
xlabel(['time (min) into nc14'])

ax = gca;
ax.YAxis(1).Color = green;
ax.YAxis(2).Color = magenta;

pbaspect([2 1 1])


%% Plot supplemental sample traces

spot_num = [2053,2020,2008,2002,1158,1142,1127,657,652,477,463,460,456,452,448,445,444,442,437,436,431,428,424,397];


time_vec_final = 0:20:1800;
fluo_interp_vec_final = zeros(size(time_vec_final));

for i = 1:length(spot_num)
    fluo_interp_vec_final = zeros(size(time_vec_final));
    % extract core vectors 
    fluo_vec = spot_struct(spot_num(i)).fluo;
    time_vec = spot_struct(spot_num(i)).time;
    knirps_vec = spot_struct(spot_num(i)).rawNCProtein-knirps_offset;

    ap_interp_vec = spot_struct(spot_num(i)).APPosNucleus;  
    mean_ap = nanmean(ap_interp_vec);
    
    fluo_interp_vec = spot_struct(spot_num(i)).fluoInterp;
    time_interp_vec = spot_struct(spot_num(i)).timeInterp;
    knirps_interp_vec = spot_struct(spot_num(i)).rawNCProteinInterp-knirps_offset;
    
    start_index = time_interp_vec(1)/20+1;
    final_index = time_interp_vec(end)/20+1;
    
    fluo_interp_vec_final(start_index:final_index) = fluo_interp_vec;
    
    knirps_vec = movmean(knirps_vec,3);
    fluo_vec_final = movmean(fluo_vec,3);
    fluo_interp_vec_final = movmean(fluo_interp_vec_final,3);
   
    if abs(mean_ap)<0.02
        fig1 = figure;
        hold on
        yyaxis left
        plot(time_vec/60,knirps_vec,'-','LineWidth',2,'MarkerSize',15,'Color',green)
        scatter(time_vec/60,knirps_vec,50,'MarkerFaceColor',green,'MarkerEdgeColor','k','Color',magenta)
        ylim([1.25E5 13.25E5])
        ylabel('[Knirps] (AU)')
        yyaxis right
        plot(time_vec_final/60,fluo_interp_vec_final/max(fluo_interp_vec_final),'-','LineWidth',2,'MarkerSize',15,'Color',magenta)
        scatter(time_vec_final/60,fluo_interp_vec_final/max(fluo_interp_vec_final),50,'MarkerFaceColor',magenta,'MarkerEdgeColor','k')
        
        xlim([5 27.5])
        ylim([-0.05 1.1])
        ylabel('transcription rate (normalized)')
        
        plot(time_plot,fluo_vec_mean/fluo_vec_mean(7),'-k','LineWidth',1)
        pbaspect([3 1 1])
    
        ax = gca;
        ax.YAxis(1).Color = green;
        ax.YAxis(2).Color = magenta;
    end
      
end
