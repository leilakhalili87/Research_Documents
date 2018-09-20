clear all
close all
%%%%reading file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PathStr=('/home/leila/Leila_sndhard/codes/gb_mobility_test/gb_mob/mom_E0.1_T800/dE0.3/f0.01/zdump_mom_on_type4and5');
S = dir(fullfile(PathStr,'zdump.*.out'));
[~,ndx] = natsortfiles({S.name});
S = S(ndx); % sort structure using indices
%for k = 1:numel(S)
for k=1:1000
    c=[];cc=[];a=[];b=[];
    x_min=-50;
    x_max=100;
    o_f_min=.35;
    o_f_max=.7;
    xx_max=100;
    xx_min=80;
    
    %%read data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(fullfile(PathStr,S(k).name));
    fgetl(fid);
    tline_1 = fgetl(fid);
    t = textscan(tline_1, '%f');
    %fgetl(fid);
    fgetl(fid);
    tline = fgetl(fid);
    a = textscan(tline, '%d');
    atoms = a{1};
    fgetl(fid);
    tline = fgetl(fid);
    a = textscan(tline,'%f %f');
    xlo  = a{1}; xhi = a{2};
    tline = fgetl(fid);
    a = textscan(tline,'%f %f');
    ylo  = a{1}; yhi = a{2};
    tline = fgetl(fid);
    a = textscan(tline,'%f %f');
    zlo  = a{1}; zhi = a{2};
    fgetl(fid);
    B = textscan(fid, '%d %*d %f %f %f %f %f %f', atoms);
    atomID = B{1};
    Xpos = B{2};
    Ypos = B{3};
    Zpos = B{4};
    orien_p=B{5};
    %M=[Xpos,Ypos,Zpos];
    fclose(fid);
    
    m=[Xpos,orien_p];
    
    
    a=find(m(:,1)<x_max & m(:,1)>x_min );
    b=find (m(a,2)<o_f_max & m(a,2)>o_f_min);
    c=m(a(b),:);

    for i=1:15
        qq=find(c(:,1)<xx_max & c(:,1)>xx_min);
        [n_qq,~]=size(qq);
        if n_qq<100
        xx_max=xx_max-10;
        xx_min=xx_min-10;
        else
           break
        end
    end
    cc=c(find(c(:,1)<xx_max & c(:,1)>xx_min),:);
    cc(:,1)=cc(:,1);
    f=fit(cc(:,1),cc(:,2),'poly1');
    ci=confint(f);
    p1=mean(ci(:,1));
    p2=mean(ci(:,2));
    t_1(k)=(t{1});
    x_111(k)=mean(cc(:,1));
    x_1(k)=(.55-p2)/p1;
end

%scatter(t_1,x_1)
% hold on
%scatter(t_1,x_11)
%x_1=x_1-x_1(1);

%   mob=fit(t_1',x_1','poly1');
%     c_mob=confint(mob);
%    
%     p_1_mob=mean(c_mob(:,1))*100
     %ave_p=(max(x_1)-min(x_1))*yhi*zhi*100/ 160.21766208/max(t_1)
figure('DefaultAxesFontSize',18)
 scatter(t_1,x_111,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
ylim([-10 100])
 xlabel('Time (fs)') % x-axis label
ylabel('Position of GB in X direction (A)') % y-axis label