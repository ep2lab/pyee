function mf_array = Coil_HPT05M(ss_path)

mf = cell(1,47);
loop_pos = linspace(-35 + 120 ,35 + 120, 48);

for i=1:length(loop_pos)

mf{i} = magnetic_field.loop_2d('RL',30/1000,'Zl',loop_pos(i)/1000,'I',1);

end

mf_array = magnetic_field.array_2d;
mf_array.generators = mf;

%% Check that the magnetic fields are equal
R = double(h5read(ss_path,'/picM/rs')); 
Z = double(h5read(ss_path,'/picM/zs')); 
BrM = double(h5read(ss_path,'/picM/Br')); 
BzM = double(h5read(ss_path,'/picM/Bz')); 

[~,Bz,Br] = mf_array.field_2d(Z,R);
B         = sqrt(Bz.^2+Br.^2);
BM        = sqrt(BzM.^2 + BrM.^2);

scale = mean(BM(:))/mean(B(:));

for i=1:length(loop_pos)

mf_array.generators{i}.I  = mf_array.generators{i}.I * scale;

end

[~,Bz,Br] = mf_array.field_2d(Z,R);
B         = sqrt(Bz.^2+Br.^2);

ampere_turn = scale * length(loop_pos);
B_axis      = B(abs(Z-0.12)<1e-6 & R==0);

figure
surf(Z,R,B); view(2);shading interp;hold on; colorbar;

%% Check single loop
mf2 = magnetic_field.loop_2d('RL',30/1000,'Zl',120/1000,'I',ampere_turn);
[~,Bz,Br] = mf2.field_2d(Z,R);
B         = sqrt(Bz.^2+Br.^2);
B_axis2   = B(abs(Z-0.12)<1e-6 & R==0);

end