c = constants_and_units.constants;

sigma = 27e-20; %Xenon
pc = 1e-3;
T0 = 273 + 20;

nn = pc/(c.kB*T0);

ce = sqrt(8 * c.kB * T0 / pi / c.me);


nn = 3e19;
nu_n = ce * nn * sigma