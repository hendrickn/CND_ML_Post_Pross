clear all; clc; close all;

RX = [0.08:0.015:0.52];

for i=1:30
filename=['3D' num2str(i) 'GR.in'];
fid1 = fopen(filename,'w');

fprintf(fid1,'#domain: 1.00 0.35 0.60 \n');
fprintf(fid1,'#dx_dy_dz: 0.005 0.005 0.005 \n');
fprintf(fid1,'#time_window: 5e-9 \n\n');

fprintf(fid1,'#material: 24 0.00063 1 0 racine \n');
fprintf(fid1,'#waveform: gaussiandotnorm 1 1e9 pulse \n\n');

fprintf(fid1,'#hertzian_dipole: z 0.08 0.30 %f pulse \n',RX(i));
fprintf(fid1,'#rx: 0.18 0.30 %f \n',RX(i));

fprintf(fid1,'#src_steps: 0.02 0 0 \n');
fprintf(fid1,'#rx_steps: 0.02 0 0 \n\n');

fprintf(fid1,'#soil_peplinski: 0.7 0.3 2.0 2.66 0.001 0.10 my_soil \n');
fprintf(fid1,'#fractal_box: 0 0 0 1.00 0.30 0.60 1.5 1 1 1 20 my_soil my_fractal_box 3 \n');
%fprintf(fid1,'#add_surface_roughness: 0 0.30 0 1.0 0.30 0.40 1.5 1 1 0.29 0.31 my_fractal_box 3 \n\n');

fprintf(fid1,'#cylinder: 0.40 0.15 0 0.60 0.15 0.20 0.025 racine \n');
fprintf(fid1,'#cylinder: 0.60 0.15 0.20 0.40 0.15 0.40 0.025 racine \n');
fprintf(fid1,'#cylinder: 0.40 0.15 0.40 0.60 0.15 0.60 0.025 racine \n');

fprintf(fid1,'#geometry_view: 0 0 0 1.00 0.35 0.60 0.005 0.005 0.005 3Dcomproot n \n');

fclose(fid1);
end