clear all;
close all;

%%
%%%%%structure property%%%%%
a=1;  % lattice parameter
ratio=1.624; % ratio of c/a
%e_i=[1,0,.5;0,1,.5;0,0,.5];  %lattice basis for primitive bcc
%e_i=[a,-.5*a,0;0,a*sqrt(3)/2,0;0,0,ratio];  %lattice basis for primitive hcp
%e_i=[1,0,0;0,1*sqrt(3),0;0,0,ratio]  %lattice basis for multilattice hcp
e_i=[.5,0,.5;.5,.5,0;0,.5,.5];  %lattice basis for primitive fcc


%%
vol=dot(e_i(:,1),cross(e_i(:,2),e_i(:,3)));
%reciprocal lattice basis
e_r(:,1)=cross(e_i(:,2),e_i(:,3))/vol;
e_r(:,2)=cross(e_i(:,3),e_i(:,1))/vol;
e_r(:,3)=cross(e_i(:,1),e_i(:,2))/vol;

kapa=-1;
n = 1; %
[mo,count]=equivalent(n); %generating acceptable equivalent matrices

cc=1;jj=1;tt=1;zz=1;
for i=1:count
    g_i=e_i*mo{i};
    def_grad=kron(g_i(:,1),e_r(:,1)')+kron(g_i(:,2),e_r(:,2)')+kron(g_i(:,3),e_r(:,3)'); %deformation gradient
    H=kron(g_i(:,1),e_r(:,1)')+kron(g_i(:,2),e_r(:,2)')+kron(g_i(:,3),e_r(:,3)');
    C=def_grad'*def_grad; %Right Cauchy Green deformation tensor
    [V,D]=eig(C); 
    if D(1,1)>0 && D(2,2)>0 && D(3,3)>0
   if D(1,1)<1 && abs(D(2,2)-1)<10^-12 && D(3,3)>1
            h{cc}=H;
            ss=sqrt(D(3,3)*(1-D(1,1))/(D(3,3)-D(  1,1)))*V(:,1)+kapa*V(:,3)*sqrt((D(1,1)*(D(3,3)-1))/(D(3,3)-D(1,1)));
            nnn=(sqrt(D(3,3))-sqrt(D(1,1)))/sqrt(D(3,3)-D(1,1))*(-sqrt(1-D(1,1))*V(:,1)+kapa*sqrt(D(3,3)-1)*V(:,3));
            k_1{cc}=nnn/norm(nnn);
            eta_1{cc}=ss*norm(nnn);
            shear(cc)=norm(eta_1{cc});
            Q{cc}=(eye(3)+kron(eta_1{cc},k_1{cc}'))*inv(H);
            % dett(cc)=det(Q{cc});
            if norm(Q{cc}'*Q{cc}-eye(3))>10^-5
                cc
                i
            end
            verom=vrrotmat2vec(det(Q{cc})*Q{cc});
            teta(cc)=real(verom(4)*180/pi);
            cc=cc+1;
        end
    end
end


histogram(shear,100,'FaceColor','r')
xlabel('Shear', 'FontSize', 20);ylabel('Tally', 'FontSize', 20);
set(gca,'FontSize',20); set(gca,'xtick',[0:1:21]);
figure;
histogram(teta,180,'FaceColor','r')
set(gca,'FontSize',20);set(gca,'xtick',[0:45:180]);
xlabel('\theta [\circ]', 'FontSize', 20);ylabel('Tally', 'FontSize', 20);