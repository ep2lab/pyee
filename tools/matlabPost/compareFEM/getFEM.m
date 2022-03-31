function field = getFEM(fieldTag, Z, R)

mesh = h5read([fieldTag,'.h5'],'/Mesh/mesh/geometry');

ZFEM = mesh(1,:);
RFEM = mesh(2,:);

real = h5read([fieldTag,'.h5'],'/VisualisationVector/real/0');
imag = h5read([fieldTag,'.h5'],'/VisualisationVector/imag/0');

field = real + 1i*imag;

F = scatteredInterpolant(ZFEM', RFEM' ,field');


% Structured mesh
m = size(Z, 1);
n = size(Z, 2);

field = reshape(F(Z(:),R(:)),[m,n]);