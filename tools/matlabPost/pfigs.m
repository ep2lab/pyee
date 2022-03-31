Ptot = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, P.*Rr(:,1)*0.01));

Pp = P; Pp(~PointFlag) = 0; Pp(Zr<32.5) = 0;
Pplume = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, Pp.*Rr(:,1)*0.01));

Pc = P; Pc(PointFlag) = 0;
Pchamb = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, Pc.*Rr(:,1)*0.01));

Pt = P; Pt(~PointFlag) = 0; Pt(Zr>32.5) = 0;
Ptube = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, Pt.*Rr(:,1)*0.01));

Po = P; Po(omega_ce>1) = 0;
Pout = 2*pi*trapz(Zr(1,:)*0.01, trapz(Rr(:,1)*0.01, Po.*Rr(:,1)*0.01));


fprintf('\n Total power: %f [W]',Ptot)
fprintf('\n Tube power fraction: %f',Ptube/Ptot)
fprintf('\n Plume power fraction: %f',Pplume/Ptot)
fprintf('\n Chamber power fraction: %f',Pchamb/Ptot)