%tic
clear all, clc, close all;
c = 299792458/sqrt(5.63);
raw_data=zeros(520,37,30);
for r=1:29
    filename=['3D' num2str(r) 'GR_merged.out'];
    
    %iterations = double(h5readatt(filename, '/', 'Iterations'));
    iterations = 520;
    dt = h5readatt(filename, '/', 'dt');
    %time samples
    t_samples = (0:iterations-1)*dt;
    %range samples
    range = t_samples*c/2;
    %spatial step
    dx = 0.02;
    % number of pulses
    np = 37;
    %antenna position in crossF range
    antx = (0:np-1)*dx;
    fieldpath = strcat('/rxs/rx1/Ez');
    fieldc = h5read(filename, fieldpath)';
    field = fieldc(1:520,:);
    %%% Background subtraction 1 %%%
    %figure; imagesc(antx,range,field);
    Escat=(field-repmat(mean(field'),37,1)');
    %figure; imagesc(antx,range,Escat);
    %%% Background subtraction 2 %%%
    %FF = field;
    FF = Escat;
    raw_data(:,:,r)=FF;
%     %%%% SVD Technique %%%%
%     [U,S,V]=svd(FF);
%     s=svd(FF);
%     I=0;
%     %bar(s)
%     k=length(s);
%     %j = input('Number of SVD components to remove: ');
%     for i=3:k
%         I = I + s(i)*U(:,i)*V(:,i)';
%     end
    %figure;imagesc(antx,range,I)

%      Res = I;
%     %============= Match Filtering ===============
%     %image pixels (note size of image and scene are the same with cross range to centre distance equal to zero)
%     pix_x = antx;
%     pix_y = range;
%     pixels = zeros(length(pix_y),length(pix_x));
%     %array to hold time delay to each pixel image
%     delay = zeros(size(pixels));
%     %array to hold reference signal for each pixel
%     refsig = zeros(size(pixels));
%     %frequency
%     f = 2e9;
%     for ii=1:np  % for each pulse
%         
%         fprintf('Pulse %d of %d\n',ii,np);
%         
%         for jj=1:iterations   % for each range bin
%             %Calculate time delays for pixel(jj,ii)
%             for k=1:np
%                 t = sqrt((antx(k)-pix_x(ii))^2 + pix_y(jj)^2)*2/c;
%                 delay(:,k) = t_samples-t;
%             end
%             
%             % Calculate reference signal for 1 amp gaussian wave
%             a = sqrt(2)*pi*f;
%             b = (delay-1/f);
%             d = -2.*sqrt(exp(1)/(2.*a^2)).*a^2.*b;
%             refsig = d.*exp(-a.^2*b.^2);
%             
%             % perform 2D correlation matched filter
%             pixels(jj,ii) = sum(sum(Res.*refsig));
%         end
%         
%     end
%     %Total field
%     subplot(2,2,1);
%     imagesc(antx,range,field)
%     title('Total field')
%     %BG subtraction
%     subplot(2,2,2);
%     imagesc(antx,range,Res);
%     title('Mean Subtraction + SVD');
%     % MF
%     subplot(2,2,3);
%     imagesc(antx,range,pixels);
%     title('Application of MF');
%     % MF + HT
%     HT= abs(hilbert(pixels));
%     subplot(2,2,4);
%     imagesc(antx,range,HT);
%     title('MF + HT');
%     
%     M(:,:,r) = HT;
%     
end
save('raw_data.mat','raw_data');
% M_norm=M./max(max(max(abs(M))));
% 
% [X,Y,Z] = meshgrid(pix_x,pix_y,1:1:30);
% [Xq,Yq,Zq] = meshgrid(min(pix_y):0.001:max(pix_x),pix_y,1:0.1:30);
% 
% M_interp=interp3(X,Y,Z,M_norm,Xq,Yq,Zq,'spline');
% 
% 
% MAX=max(max(max(M_interp)));
% SAR=zeros(size(M_interp));
% 
% for x=100:size(M_interp,1)
%     for y=1:size(M_interp,2)
%         for z=1:size(M_interp,3)
%             if(M_interp(x,y,z)>=0.25)
%                 SAR(x,y,z)=1;
%             end
%         end
%     end
% end
% 
% Chi_Sigma=SAR;
% IsoLevel=0.7;
% 
% axes1 = axes(...
%   'FontName','Times',...
%   'FontSize',14,...
%   'FontWeight','bold',...
%   'LineWidth',1,...
%   'Parent',gcf);
%     N=size(Chi_Sigma);
%     A=isosurface(Chi_Sigma,IsoLevel)
%     p = patch(isosurface(Chi_Sigma,IsoLevel));
%     xlabel('X Axis')
%     ylabel('Y Axis')
%     zlabel('Z Axis')
% 
%     set(p,'FaceColor','red','EdgeColor','none','Parent',axes1);%
%     
%     camlight ;
%     lighting gouraud;
