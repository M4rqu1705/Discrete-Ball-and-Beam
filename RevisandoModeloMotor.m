s = zpk('s')

T_s = 30E-3;
G_m = 61 / (s * (s + 36));
G_md = c2d(G_m, T_s);
K_m = 7;


lazo_interno_d = minreal(K_m * G_md / (1 + K_m * G_md));

G_bb = 0.4183 / s^2;
G_bbd = c2d(G_bb, T_s, 'zoh');

lazo_abierto = zpk(G_bbd * lazo_interno_d)