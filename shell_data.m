% Script run with following version of MATLAB and Chebfun
% MATLAB Release R2016b
% Chebfun Version 5.5.0 01-Jul-2016

% End point stress = +5.749001e-11 (R1 = 7.935607e+00)
% End point stress = -3.482525e-11 (R2 = 1.450776e+01)

R1 =  7.9356072789184839e+00; % Inner Radius
R2 =  1.4507758856012858e+01; % Outer Radius
N = 26; % Number of Chebyshev interpolation points

% Interpolation Points
rn   = [+7.9356072789184839e+00; +7.9615189680144418e+00;
        +8.0388453924863779e+00; +8.1663670684278031e+00;
        +8.3420729028313634e+00; +8.5631919097269407e+00;
        +8.8262369102856955e+00; +9.1270595277133530e+00;
        +9.4609156096303515e+00; +9.8225400461889123e+00;
        +1.0206229804000534e+01; +1.0605933866379377e+01;
        +1.1015348661490094e+01; +1.1428017473441248e+01;
        +1.1837432268551964e+01; +1.2237136330930809e+01;
        +1.2620826088742430e+01; +1.2982450525300987e+01;
        +1.3316306607217991e+01; +1.3617129224645646e+01;
        +1.3880174225204403e+01; +1.4101293232099978e+01;
        +1.4276999066503539e+01; +1.4404520742444964e+01;
        +1.4481847166916900e+01; +1.4507758856012858e+01];

% Barycentric Weights
wn   = [-5.0000000000000000e-01; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +1.0000000000000000e+00;
        -1.0000000000000000e+00; +5.0000000000000000e-01];

% Lame's first parameter
lam  = [+4.5223165114926953e+00; +4.5106135942143997e+00;
        +4.4758691384246738e+00; +4.4191544909522520e+00;
        +4.3421855430725973e+00; +4.2472253122577861e+00;
        +4.1369607586155288e+00; +4.0143648484323942e+00;
        +3.8825553492768892e+00; +3.7446609857608362e+00;
        +3.6037036655273713e+00; +3.4625028906438864e+00;
        +3.3236056303039723e+00; +3.1892422380563561e+00;
        +3.0613067464721202e+00; +2.9413582365841280e+00;
        +2.8306390082074935e+00; +2.7301049193752047e+00;
        +2.6404634005667291e+00; +2.5622151326266529e+00;
        +2.4956960555437537e+00; +2.4411171180185858e+00;
        +2.3985998861341633e+00; +2.3682067395156343e+00;
        +2.3499648634064569e+00; +2.3438835897199599e+00];

% Lame's second parameter / shear modulus
mu   = [+2.7877506026110423e+00; +2.7837289312245326e+00;
        +2.7717960275232603e+00; +2.7523389685820776e+00;
        +2.7259732454929773e+00; +2.6935022943872315e+00;
        +2.6558679909463243e+00; +2.6140975787498646e+00;
        +2.5692521874591141e+00; +2.5223810815109449e+00;
        +2.4744843785801516e+00; +2.4264855125823717e+00;
        +2.3792134408706382e+00; +2.3333936544513287e+00;
        +2.2896464833919552e+00; +2.2484909592631483e+00;
        +2.2103525222221312e+00; +2.1755730506177251e+00;
        +2.1444219636621864e+00; +2.1171074402989900e+00;
        +2.0937870690753737e+00; +2.0745774726065784e+00;
        +2.0595626281086807e+00; +2.0488007335886231e+00;
        +2.0423295537751174e+00; +2.0401702291000281e+00];

% Density
rho  = [+2.2702956758846300e+00; +2.2859339416458995e+00;
        +2.3321353538657470e+00; +2.4065646672183609e+00;
        +2.5047479853879797e+00; +2.6194352992444414e+00;
        +2.7401539191516555e+00; +2.8533155840194984e+00;
        +2.9431818490871375e+00; +2.9937901590993405e+00;
        +2.9916165100693490e+00; +2.9283975266280247e+00;
        +2.8033109498131670e+00; +2.6237618885662091e+00;
        +2.4043825180298932e+00; +2.1644122413808997e+00;
        +1.9241556324412183e+00; +1.7014885561790885e+00;
        +1.5092933554092538e+00; +1.3543203546356015e+00;
        +1.2374866610368698e+00; +1.1552369609547355e+00;
        +1.1014205469502767e+00; +1.0691816209812113e+00;
        +1.0525298428587604e+00; +1.0474483229806277e+00];

% radial displacement at t = 0
ur   = [-1.8141315805162104e+00; -1.8087270254955574e+00;
        -1.7914142047528374e+00; -1.7589219636049498e+00;
        -1.7059323217024207e+00; -1.6254386376388201e+00;
        -1.5095333835610301e+00; -1.3508097780336068e+00;
        -1.1444035094640215e+00; -8.9035853620930017e-01;
        -5.9556347927405595e-01; -2.7424757060932403e-01;
        +5.3707574016091997e-02; +3.6650025558973798e-01;
        +6.4497079701031934e-01; +8.7626244468944514e-01;
        +1.0553211331679795e+00; +1.1841036634337190e+00;
        +1.2693541372590622e+00; +1.3201050903742966e+00;
        +1.3457157444147596e+00; +1.3547011169431924e+00;
        +1.3542030048198126e+00; +1.3498205882244267e+00;
        +1.3455738487067426e+00; +1.3438878054380818e+00];

% radial derivative of radial displacement at t = 0
ur_r = [+2.0476326240114867e-01; +2.1239678281621188e-01;
        +2.3545616438280964e-01; +2.7429351321291018e-01;
        +3.2906815405814444e-01; +3.9908130380533852e-01;
        +4.8187317725049689e-01; +5.7227635719280767e-01;
        +6.6184058531130951e-01; +7.3920599714675173e-01;
        +7.9187667778713977e-01; +8.0926578946446770e-01;
        +7.8601444224394534e-01; +7.2402943295538169e-01;
        +6.3206473469544022e-01; +5.2294117391345130e-01;
        +4.0977729837511934e-01; +3.0294420920687604e-01;
        +2.0873949306477285e-01; +1.2970589183840769e-01;
        +6.5852119913216017e-02; +1.5983371405877356e-02;
        -2.1329223739578795e-02; -4.7249797058461107e-02;
        -6.2538220304036191e-02; -6.7593983321332235e-02];
