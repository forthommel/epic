const size_t c_tcs_t_n = 10;

const double c_tcs_t[c_tcs_t_n][4] = {
	{-1, -0.90001, 6.55601969507633e-05, 2.32781536213997e-06},
	{-0.90001, -0.80002, 9.75939065765482e-05, 3.87175660573723e-06},
	{-0.80002, -0.70003, 0.000116422534975832, 7.07149590056436e-06},
	{-0.70003, -0.60004, 0.000182696072819055, 6.20996872650241e-06},
	{-0.60004, -0.50005, 0.000256425972613311, 2.12721798675027e-05},
	{-0.50005, -0.40006, 0.000429723031786651, 2.43730911383137e-05},
	{-0.40006, -0.30007, 0.000473093194354399, 3.12017578381016e-05},
	{-0.30007, -0.20008, 0.000857078532954347, 4.53934850197056e-05},
	{-0.20008, -0.10009, 0.00112432011570497, 7.35558431494001e-05},
	{-0.10009, -0.0001, 0.00137217284527211, 0.000158880539850254}
};

const size_t c_tcs_QPrim2_n = 10;

const double c_tcs_QPrim2[c_tcs_QPrim2_n][4] = {
	{2, 2.8, 0.00317144926554602, 0.00027046829754705},
	{2.8, 3.6, 0.000915809036294391, 9.34393845518818e-05},
	{3.6, 4.4, 0.000320952208665911, 3.33820358213038e-05},
	{4.4, 5.2, 0.000161221550615195, 1.87459861648839e-05},
	{5.2, 6, 0.000161504423137169, 2.42044332279804e-05},
	{6, 6.8, 5.13735757426015e-05, 6.22012882106756e-06},
	{6.8, 7.6, 3.00299024343722e-05, 2.77396113904273e-06},
	{7.6, 8.4, 1.29688163533164e-05, 1.38815027694446e-06},
	{8.4, 9.2, 1.26530060073274e-05, 1.40747407306789e-06},
	{9.2, 10, 5.98165918163344e-06, 6.12583116719993e-07}
};

const size_t c_tcs_phiL_n = 10;

const double c_tcs_phiL[c_tcs_phiL_n][4] = {
	{0.05, 0.66831853, 0.000715643952916942, 5.81927351337994e-05},
	{0.66831853, 1.28663706, 0.00034189737954421, 2.89342471388395e-05},
	{1.28663706, 1.90495559, 0.000417284039375375, 3.64337569797944e-05},
	{1.90495559, 2.52327412, 0.000269183987201652, 3.29921860256507e-05},
	{2.52327412, 3.14159265, 0.00143216286817972, 0.000203360762312387},
	{3.14159265, 3.75991118, 0.00176181040700952, 0.000211009455583771},
	{3.75991118, 4.37822971, 0.000375923810340841, 4.69453231085265e-05},
	{4.37822971, 4.99654824, 0.000318418783943235, 3.14110726229511e-05},
	{4.99654824, 5.61486677, 0.000350578468983947, 3.64439698646642e-05},
	{5.61486677, 6.2331853, 0.000719342743798803, 9.60824374870119e-05}
};

const size_t c_tcs_thetaL_n = 10;

const double c_tcs_thetaL[c_tcs_thetaL_n][4] = {
	{0.7853, 0.94238, 0.000681579722914401, 7.27992166949909e-05},
	{0.94238, 1.09946, 0.000464235036459225, 5.240502534577e-05},
	{1.09946, 1.25654, 0.000314037853302143, 3.18275086824754e-05},
	{1.25654, 1.41362, 0.000325832711323657, 3.8859801283565e-05},
	{1.41362, 1.5707, 0.000277467768408063, 2.15564411382838e-05},
	{1.5707, 1.72778, 0.000376277897403814, 6.85078027390992e-05},
	{1.72778, 1.88486, 0.000481207166586306, 8.10074769780025e-05},
	{1.88486, 2.04194, 0.000406618541661646, 6.87970063542988e-05},
	{2.04194, 2.19902, 0.000483081758091715, 2.75260880302264e-05},
	{2.19902, 2.3561, 0.000659292687677637, 3.74219715288356e-05}
};

const size_t c_tcs_y_n = 10;

const double c_tcs_y[c_tcs_y_n][4] = {
	{0.05, 0.14, 0, 0},
	{0.14, 0.23, 0.000297558228362713, 2.57633908642573e-05},
	{0.23, 0.32, 0.000519538382297486, 2.64207432904595e-05},
	{0.32, 0.41, 0.00066917622380798, 2.55346228869967e-05},
	{0.41, 0.5, 0.00066943298720003, 2.88774433062688e-05},
	{0.5, 0.59, 0.000662504503318649, 2.06481777656942e-05},
	{0.59, 0.68, 0.000721444421531986, 2.95151992211849e-05},
	{0.68, 0.77, 0.000588675852095423, 3.21439684623914e-05},
	{0.77, 0.86, 0.000561808166671067, 1.96971422114418e-05},
	{0.86, 0.95, 0.000571151154131004, 3.07039136083263e-05}
};

const size_t c_tcs_Q2_n = 10;

const double c_tcs_Q2[c_tcs_Q2_n][4] = {
	{1e-06, 0.0050009, 0.00377950783007817, 0.000293970068927207},
	{0.0050009, 0.0100008, 0.000337976583274448, 2.34271042981231e-05},
	{0.0100008, 0.0150007, 0.000175714406112716, 1.42017413624024e-05},
	{0.0150007, 0.0200006, 0.000117198463393425, 1.34431807704268e-05},
	{0.0200006, 0.0250005, 9.18978567865126e-05, 1.14139494537149e-05},
	{0.0250005, 0.0300004, 0.000102479082623226, 1.0895787720084e-05},
	{0.0300004, 0.0350003, 7.71189644919791e-05, 6.66908479712584e-06},
	{0.0350003, 0.0400002, 6.16053694884647e-05, 6.69234278638629e-06},
	{0.0400002, 0.0450001, 4.66091319083621e-05, 2.85417994500276e-06},
	{0.0450001, 0.05, 5.05328409649155e-05, 4.69375879551902e-06}
};