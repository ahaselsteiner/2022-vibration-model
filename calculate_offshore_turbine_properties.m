m1 = 800 * 10^3;
m2 = 400 * 10^3;
M = m1 + m2;
T = 3;
k1 = M * ((2 * pi) / T)^2;

thickness = 0.023;
r = 3;
IT = 2 * pi * r^3*thickness; % https://www.bau.uni-siegen.de/subdomains/bauinformatik/lehre/tm2/arbeitsblaetter/arbeitsblatt_29_grundformeln_zur_torsion.pdf
GSteel = 80 * 10^9; % https://de.wikipedia.org/wiki/Schubmodul
height = 100;
k2 = GSteel * IT / height; % https://de.wikipedia.org/wiki/Torsion_(Mechanik)


f0_bending = 1 / (2 * pi) * sqrt(k1 / (m1 + m2))
f0_torsion = 1 / (2 * pi) * sqrt(k2 / (I + m2 * d^2))
