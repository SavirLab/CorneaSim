function []= EquipotentSim(RLS,short_range_rate,bias_area,dPhi_input,name,div_mat_start,Simulation_length)
%% INPUTS
% 'RLS' :=  replicative life span of progenitor cells i.e maximal number of
%progenitor divisions, progenitor cells that exeeded thier RLS cannot be
%choosen for devision, only for removal

%'short_range_rate' := percentage of divisions that are "coupled" where
%deviding and removed cells can only be choosen from a small area that
%depends on "m"(r_around) so:
%  'short_range_rate' = 1 -> all devisions are coupled short range 
%  'short_range_rate' = 0 -> all devisions are uncoupled long range  

% 'bias_area' :=  Relevant for uncoupled long range interactions - the percentage of
% the area in the center the removed cell is chosen from(pi*rho^2) so :
% small bias_area -> highly biased twardes the center
% large bias_area -> no bias twards the center 

%'dPhi_input' :=  Relevant for coupled short range interactions - the
%degree(theta) of biast twards the center in radians so:
% dPhi_input = 0.2*pi -> 36 degreed centered twards the center  i.e highly biased 
% dPhi_input = 2*pi -> 360 degreed 'centered' twards the center  i.e highly biased 

%'name' :=  file name to save data and video

%''div_mat_start' :=  for the simulation to be accurate we have a
%beforehand calcuated values of steade state nuber of divisions left for
%each cell ( so some of the cells in the center are with lower number of
%divisions/postmitotic) the accompanied Files - "div_mat_start_DATA_short" and
%"div_mat_start_DATA_Long" has this precalculated data, without it the
%simulation will run but will assume all the cells have the same
%replicative capacity RLS in the initial state.

% 'Simulation_length':=  maximal simulation length in terms of simulation
% time-steps where 1 timestep = 1 duplication time. if the cornea is full
% beforehand to save resources the simulation will end after 150 timesteps
% in steady state. (can be modified in the code)

%% OUTPUTS
% matlab mat file with all the variables from the simulation, the most
% important are:
% mov_cnt :=  number of time snapshots the simulation acptured 
% time := time(in termes of division times) for for each snapshot 
% color_mat(i,j,1,k) := liniage marker number in location(i,j) of grid_boundry at time k(1 to
% mov_cnt)
% div_mat(i,j,k) := number of divisions made in location(i,j) of grid_boundry at time k(1 to
% mov_cnt)


%% examples of simulations

%Coupled short range with bias -   EquipotentSim(60,1,1,0.2*pi,'Coupled_with_bias','NaN',200)
%Coupled short range without bias -   EquipotentSim(60,1,1,2*pi,'Coupled_without_bias','NaN',200)
%Uncoupled long range with bias -   EquipotentSim(4,0,0.1,0.2*pi,'Uncoupled_with_bias','NaN',200)
%Uncoupled short range without bias -   EquipotentSim(60,0,0.999,2*pi,'Uncoupled_without_bias','NaN',200)

%% Parametrs
r =100;  % radius in number of cells
gs = 1/10; % Stem cell division rate in 1/days 
gp = 1/1; % pluripotent cell division rate 1/days

TAC_mode=1;% inhibit cells over RLS max ALWAYS ON
push_mode=1;% push cells from division to hole  ALWAYS ON
cent=1;% centripetal movment ALWAYS ON
cent_divide=cent;% centripetal division ALWAYS ON 
rejuvination_mode=0;% ALWAYS ON
Hirr_rate=0.85; % rate of asymetrical devisions
vertical_div=0;%no vertical divisions only horizontal ALWAYS OFF
r_around=5; % 'm' =: the number of cells in the coupled area 

%% Initailize Grid
[ind1,ind2,ind3,grid_boundery] = CircleGrid(r);% generate a circular grid with radius r and space
[ inner_ind ] = GenInCirc(r,grid_boundery,bias_area);% generate inner circle with area of 'bias area' 
grid = grid_boundery;
grid(ind1) = 1;
grid(ind2) = 1;
grid(ind3) = 1;

%% Initialize grid params 
%generate gp grid of division rules for P cells
grid_gp=grid_boundery;
grid_gp(ind1)=0;
grid_gp(ind2)=1;
grid_gp(ind3)=1;
grid_gp_start=grid_gp;

%generate gs grid of division rules for S cells
grid_gs=grid_boundery;
grid_gs(ind1)=1;
grid_gs(ind2)=0;
grid_gs(ind3)=0;
grid_gs=gs.*grid_gs;

%Generate grids for pushing rules
grid_circular = grid_boundery;
grid_circular(ind1) = 1;
grid_circular(ind2) = 0;
grid_circular(ind3) = 0;

grid_inside = grid_boundery; 
grid_inside(ind1) = 0;
grid_inside(ind2) = 1;
grid_inside(ind3) = 1;

%Generate grids for checking and a stoping condition
grid_inside_c=grid;
grid_inside_c(1,1)=1;
grid_inside_c= grid_inside_c==1;

%limbus ind
color = grid_boundery;
color(ind1) =  1:length(ind1); % generate distinct liniage marker for each stem cell 
color(ind2) =0;% randi([1 length(ind1)],length(ind2),1);% add in the center 
color(ind3) =0;% randi([1 length(ind1)],length(ind3),1);% add in the center 

color(1,1)=max(unique(color));%for simulation colors
color_table =  20*round(randi([1 5],length(ind1),1))+20;% 5 color lookup table for generating 5 color image from 200+ colord cells

%initialize age and division matrix
age_mat=grid;

try 
    div_mat=div_mat_start;
    grid_gp=grid_gp_start.*double(div_mat<=RLS);
catch
    disp('WARNING: steady state division matrix INPUT is missing, division matrix will be initialized as div_mat=grid ')
%    pause
    div_mat=grid;
end
%time params 

MC_steps = 4e9;%5e5*50;
dtt=0.5;
time_vec=0:dtt:Simulation_length; 
cnt=1;%step to ending simulation early  need to be <300 

close_mode=1;% holes are close 


%Generate for each point in the limbus places to divide to in the r_around
%area
[ rad_cell_limbus ] = RadiusForEach( r,grid_circular,r_around ); 

%calculate the minimal degree in short mode and
dPhi=max(0.2*pi,dPhi_input); %always on in this version 
[ phi_cell] = PhiForEach( r,grid,dPhi );

rad_cell_best=RadiusForEach2( r,grid,r_around,phi_cell);
save([name,'.mat']); %%% save inital data 
div_mat_time=[];
mov_cnt=1;
%% run
tic
for mc_cnt = 1:MC_steps

% Action rate 
% For P cells one action
% For S cell 2 action;

gs_mat=grid_gs;
gp_mat = gp.*grid_gp;

action_rate(:,:,1)=gp_mat;% P cell division
action_rate(:,:,2)=gs_mat;% S cell division

action_rates_vec = action_rate(:);
    % Find the rate that occur
    ind_rates = find(action_rates_vec);
    action_rates_vec_reduced = action_rates_vec(ind_rates);
    
    Z_rate = sum(action_rates_vec);
    P = cumsum(action_rates_vec_reduced)/Z_rate;
    p_test = rand(1);
    %find who is dividing 
    action_ind = find(P-p_test>=0, 1 );
    action_ind_vec = ind_rates(action_ind);
    [i,j,action_type]  = ind2sub(size(action_rate),action_ind_vec);
    %[i,j]= DIVIDING CELL
    %add age counter
    age_mat(i,j)=age_mat(i,j)+1;
    %add division left counter 
    
    %div from limbus cases
    div_from_limbus=0;
    if grid_boundery(i,j)==1
        div_from_limbus=1;
    end   

    switch action_type
        case 1 % if  P cell is dividing 

            %Add vertical division
            if vertical_div % ALWAYS OFF 
                push_mode=1;
                push_mode=round(rand(1)+0.5-0.05*grid_gp(i,j));
            end
            
            if push_mode % ALWAYS ON (left for testing purposes) 

                %Check if long or short mode
                close_mode=rand(1)-short_range_rate<=0;
                if close_mode 
                    allowed_ind=rad_cell_best{i,j};% hole from the area of the dividing cell in "r_around" radius
                else
                    allowed_ind=inner_ind; %holes only in the center 
                    %allowed_ind=rings_ind{randi_exp( n,mu )};%holes with exp dist declining from the center
                end
                close_mode=1;
                %
                
                G=allowed_ind;
                hole_ind=randi([1,length(G)],1); %generate random hole
                hole=G(hole_ind);
                [i_hole,j_hole]  = ind2sub(size(grid),hole);

                if ~(i_hole==i && j_hole==j)
                    [ trajec_mat ] = Pushgrid(r,grid_inside, i,j,i_hole,j_hole,cent_divide,phi_cell{i,j});%calculate trajectory of the push to th hole
                    [~,t_cols]=size(trajec_mat);

                    switch div_from_limbus %ALWAYS OFF for this type of simulation 
                        case 1 % if dividing from limbus
                            div_i=trajec_mat{2,t_cols}(1);
                            div_j=trajec_mat{2,t_cols}(2);
                            last_traj_is_dev=0;
                            if div_i==i && div_j==j
                                 last_traj_is_dev=1;
                                 div_i=trajec_mat{1,t_cols}(1);
                                 div_j=trajec_mat{1,t_cols}(2);
                            end 
                                
                            for ttt=1:(t_cols-last_traj_is_dev)
                                %check if they are both are 1 cells
                                if  ~(grid_boundery(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))==1 && grid_boundery(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))==2)
                                    color(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=color(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                                    
                                    %chrck if not loimbus cell
                                    if ~(grid_boundery(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))==1)
                                    age_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=age_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                                    div_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=div_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                                    end
                                    
                                end
                            end
                                    
                            if  ~(grid_boundery(div_i,div_j)==1 && grid_boundery(i,j)==2) 
                                    
                                    color(div_i,div_j)=color(i,j);

                                    if ~(grid_boundery(div_i,div_j)==1)
                                        if rejuvination_mode
                                        age_mat(div_i,div_j)=1;
                                        div_mat(div_i,div_j)=1;
                                        else
                                        age_mat(div_i,div_j)=1;
                                        div_mat(div_i,div_j)=1;
                                        end 
                                    end
                            end

                        case 0 %not form the limbus
                            div_i=trajec_mat{2,t_cols}(1);
                            div_j=trajec_mat{2,t_cols}(2);
                            last_traj_is_dev=0;
                            if div_i==i && div_j==j
                                 last_traj_is_dev=1;
                                 div_i=trajec_mat{1,t_cols}(1);
                                 div_j=trajec_mat{1,t_cols}(2);
                            end 
                            
                            for ttt=1:t_cols-last_traj_is_dev
                                %check if they are both are 1 cells color(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))~=0 &&
                                if  ~(grid_boundery(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))==1 && grid_boundery(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))==2) 
                                    color(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=color(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                                    age_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=age_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                                    div_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=div_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                                end
                            end
                                %last move, divide move
                            if ~(grid_boundery(div_i,div_j)==1 && grid_boundery(i,j)==2) 
                                    color(div_i,div_j)=color(i,j);

                                    if rejuvination_mode
                                        age_mat(div_i,div_j)=1;%new cell
                                        div_mat(div_i,div_j)=1;%new cell(rejuvination)
                                    else
                                        age_mat(div_i,div_j)=1;%new cell
                                        div_mat(div_i,div_j)=div_mat(i,j)+1;%new cell with no rejuvination
                                    end
                                    div_mat(i,j)=div_mat(i,j)+1;
                            end

                    end
                    if TAC_mode
                    %cancel division for cells bigger than RLS
                    grid_gp=grid_gp_start.*double(div_mat<=RLS);
                   
                    end
                end
                
            else %vertical devision do nothing horizontally
                
            end
            
        case 2 % S cell
            div_to_p=rand(1)-Hirr_rate<0;
            if div_to_p %if its 1  S-> S+P division

                %Find/make hole in Cornea
                if close_mode
                    allowed_ind=rad_cell_best{i,j};
                else %holes only in the center!
                    allowed_ind=inner_ind; %holes only in the center 
                    %allowed_ind=rings_ind{randi_exp( n,mu )};%holes with exp dist declining from the center
                end
                G=allowed_ind;
                hole_ind=randi([1,length(G)],1); %generate random hole
                hole=G(hole_ind);
                [i_hole,j_hole]  = ind2sub(size(grid),hole);
                
                if ~(i_hole==i && j_hole==j)
                    [ trajec_mat ] = Pushgrid(r,grid_inside, i,j,i_hole,j_hole,cent_divide,phi_cell{i,j});%calculate trajectory of the push to th hole
                    [~,t_cols]=size(trajec_mat);
                    %divide from limbus to 1P and 1S stays
                    div_i=trajec_mat{2,t_cols}(1);
                    div_j=trajec_mat{2,t_cols}(2);
                    last_traj_is_dev=0;
                    if div_i==i && div_j==j
                         last_traj_is_dev=1;
                         div_i=trajec_mat{1,t_cols}(1);
                         div_j=trajec_mat{1,t_cols}(2);
                    end 
                    for ttt=1:(t_cols-last_traj_is_dev)
                        %check if they are both are "2" cells 
                        if  ~(grid_boundery(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))==1 && grid_boundery(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))==2)
                            color(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=color(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));

                            %check if not limbus cell
                            if ~(grid_boundery(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))==1)
                            age_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=age_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                            div_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=div_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                            end

                        end
                    end
                    %last push move== division move 
                    if  ~(grid_boundery(div_i,div_j)==1 && grid_boundery(i,j)==2) 
                            color(div_i,div_j)=color(i,j);
                            if ~(grid_boundery(div_i,div_j)==1)
                                if rejuvination_mode
                                    age_mat(div_i,div_j)=1;
                                    div_mat(div_i,div_j)=1;
                                else
                                    age_mat(div_i,div_j)=1;
                                    div_mat(div_i,div_j)=1;
                                end 
                            end
                    end
                    if TAC_mode
                    %cancel division for cells bigger than TAC
                    grid_gp=grid_gp_start.*double(div_mat<=RLS);
                    end
                 end
    
            else %if it's S to S+S

                %Find/make hole in Limbus
                if close_mode
                    allowed_ind=rad_cell_limbus{i,j};
                else
                    allowed_ind=ind1;
                end 
                G=[allowed_ind];
                hole_ind=randi([1,length(G)],1); %generate random hole
                hole=G(hole_ind);
                [i_hole,j_hole]  = ind2sub(size(grid),hole);
                
                [ trajec_mat ] = Pushgrid(r,grid_circular, i,j,i_hole,j_hole,0,phi_cell{i,j});%calculate circular trajectory of the push to the hole
                [~,t_cols]=size(trajec_mat);
                %push cells circulary to the hole
                if ~(i_hole==i && j_hole==j) && (t_cols<3.5*r)                   
                                        
                    [~,t_cols]=size(trajec_mat);
                    %divide from limbus to 1S and 1S new
                    div_i=trajec_mat{2,t_cols}(1);
                    div_j=trajec_mat{2,t_cols}(2);
                    last_traj_is_dev=0;
                    if div_i==i && div_j==j
                         last_traj_is_dev=1;
                         div_i=trajec_mat{1,t_cols}(1);
                         div_j=trajec_mat{1,t_cols}(2);
                    end 
                    for ttt=1:(t_cols-last_traj_is_dev)
                        if (grid_boundery(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))==1 && grid_boundery(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2))==1)
                        color(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=color(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                        age_mat(trajec_mat{1,ttt}(1),trajec_mat{1,ttt}(2))=age_mat(trajec_mat{2,ttt}(1),trajec_mat{2,ttt}(2));
                        end
                    end
                    %last push move== division move 
                    if  ~(grid_boundery(div_i,div_j)==1 && grid_boundery(i,j)==2)
                        color(div_i,div_j)=color(i,j);
                        age_mat(div_i,div_j)=1;%new cell
                        div_mat(div_i,div_j)=1;%new cell(rejuvination)
                    end
                                        
                end
                if TAC_mode
                    %cancel division for cells bigger than TAC
                    grid_gp=grid_gp_start.*double(div_mat<=RLS);
                end
                
            end  
    end

    %handle time
    p_time = rand(1);
    dt = 1/Z_rate*log(1/p_time);%/action_rates(action_type,action_ind_site);
    
    if mc_cnt==1
        t(mc_cnt) = dt;
    else
        t(mc_cnt) = t(mc_cnt-1)+dt;
    end
    
    
    if t(mc_cnt)>time_vec(length(time_vec))% constrain to time limit 
        break

    elseif t(mc_cnt)>=time_vec(mov_cnt)% if full constrain to extra  100 steps, and then stop. 
        color_mat(:,:,1,mov_cnt) = uint16(color);
        div_mat_time(:,:,1,mov_cnt) = uint16(div_mat);
        time(mov_cnt) = t(mc_cnt);
        mov_cnt = mov_cnt+1;
        disp([num2str(mov_cnt/length(time_vec)*100),'% is done'])% report on progress 
        %if mod(mov_cnt,10)==0;save([name,'.mat']);end
        if 0 &&( (isequal(color>0,grid_inside_c)) || (sum(sum(color>=1))/sum(sum(grid)) >0.998)) 
           cnt=cnt+1;
           if cnt>300;break;end%end early if full after extra 200 steps(50 secs)
        end
    end
end
toc
save([name,'.mat']);

%% record Video
set(gcf,'renderer','OpenGL')
v = VideoWriter([name,'.avi']);
open(v)
figure(7)
cmap = colormap('jet');
cmap(1,:) = [1 1 1];
colormap(cmap);
for i_vid=1:mov_cnt-1
    colormat_mov=color_mat(:,:,1,i_vid);
    colormat_mov2=colormat_mov;
    for iii=1:length(ind1)
        colormat_mov2(colormat_mov==iii)=(iii); % generate 5 color picture from 600+ colors 
    end
    imagesc(colormat_mov2)
    title(['Time:',num2str(time(i_vid))])
    axis equal
    m(i_vid) = getframe(7);
    writeVideo(v,m(i_vid));
end
close(v)
%% analyze and save
save([name,'.mat']);
end
