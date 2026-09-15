# EXP 3 : IIR-CHEBYSHEV-FITER-DESIGN

## AIM: 

 To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 
```py
clc;
clear;
close;

// ================== USER INPUT ==================
wp     = input('Enter the pass band frequency (Radians)= ');
ws     = input('Enter the stop band frequency (Radians)= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T      = input('Enter the value of sampling time= ');

// ================== PRE-WARPING ==================
omegap = (2/T)*tan(wp/2);
omegas = (2/T)*tan(ws/2);

// ================== FILTER ORDER ==================
Epsilon = sqrt((10^(0.1*alphap))-1);

N_calc = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) ...
         / acosh(omegas/omegap);

N = ceil(N_calc);

// ================== CUTOFF FREQUENCY ==================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));

// ================== ANALOG LOW PASS PROTOTYPE ==================
[pols, gn] = zpch1(N, Epsilon, omegac);
hs = syslin('c', gn, real(poly(pols, 's')));

// ================== BILINEAR TRANSFORMATION ==================
z = poly(0, 'z');
Hz = horner(hs, (2/T)*((1 - z^-1)/(1 + z^-1)));
Hz_sys = syslin('d', Hz);

// ================== NEAT CONSOLE OUTPUT ==================
mprintf("\n====================================================\n");
mprintf("      CHEBYSHEV TYPE-I DIGITAL LOW PASS FILTER\n");
mprintf("====================================================\n");

mprintf("\nInput Specifications:\n");
mprintf("----------------------------------------------------\n");
mprintf("Passband Frequency (wp)      = %8.4f rad/sample\n", wp);
mprintf("Stopband Frequency (ws)      = %8.4f rad/sample\n", ws);
mprintf("Passband Attenuation (Ap)    = %8.2f dB\n", alphap);
mprintf("Stopband Attenuation (As)    = %8.2f dB\n", alphas);
mprintf("Sampling Time (T)            = %8.6f sec\n", T);

mprintf("\nPrewarped Analog Frequencies:\n");
mprintf("----------------------------------------------------\n");
mprintf("Omega_p                      = %10.4f rad/sec\n", omegap);
mprintf("Omega_s                      = %10.4f rad/sec\n", omegas);

mprintf("\nFilter Order:\n");
mprintf("----------------------------------------------------\n");
mprintf("Calculated Order (N)         = %8.4f\n", N_calc);
mprintf("Rounded Filter Order         = %d\n", N);

mprintf("\nDesign Parameters:\n");
mprintf("----------------------------------------------------\n");
mprintf("Ripple Factor (Epsilon)      = %8.6f\n", Epsilon);
mprintf("Cutoff Frequency (Omega_c)   = %10.4f rad/sec\n", omegac);
mprintf("Prototype Gain               = %10.4e\n", gn);

disp("Digital Transfer Function H(z):");
disp(Hz_sys);

// ================== FREQUENCY RESPONSE ==================
[HW, w] = frmag(Hz_sys, 512);

figure();
plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency (w/pi)');
ylabel('Magnitude');
title('Frequency Response of Chebyshev Type-I Low Pass Filter');

```
## PROGRAM (HPF): 
```py
clc;
clear;
close;

// ================== USER INPUT ==================
wp     = input('Enter the pass band frequency (Radians)= ');
ws     = input('Enter the stop band frequency (Radians)= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T      = input('Enter the value of sampling time= ');

// ================== PRE-WARPING ==================
omegap = (2/T)*tan(wp/2);
omegas = (2/T)*tan(ws/2);

// ================== FILTER ORDER ==================
Epsilon = sqrt((10^(0.1*alphap))-1);

N_calc = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1)))/ acosh(omegap/omegas);

N = ceil(N_calc);

// ================== CUTOFF FREQUENCY ==================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));

// ================== ANALOG LOW PASS PROTOTYPE ==================
[pols, gn] = zpch1(N, Epsilon, omegac);
hs = syslin('c', gn, real(poly(pols, 's')));

// ================== LPF → HPF TRANSFORMATION ==================
s = poly(0, 's');
hs_hpf = horner(hs, (omegac^2 / s));

// ================== BILINEAR TRANSFORMATION ==================
z = poly(0, 'z');
Hz = horner(hs_hpf, (2/T)*((1 - z^-1)/(1 + z^-1)));
Hz_sys = syslin('d', Hz);

// ================== NEAT CONSOLE OUTPUT ==================
mprintf("\n====================================================\n");
mprintf("      CHEBYSHEV TYPE-I DIGITAL HIGH PASS FILTER\n");
mprintf("====================================================\n");

mprintf("\nInput Specifications:\n");
mprintf("----------------------------------------------------\n");
mprintf("Passband Frequency (wp)      = %8.4f rad/sample\n", wp);
mprintf("Stopband Frequency (ws)      = %8.4f rad/sample\n", ws);
mprintf("Passband Attenuation (Ap)    = %8.2f dB\n", alphap);
mprintf("Stopband Attenuation (As)    = %8.2f dB\n", alphas);
mprintf("Sampling Time (T)            = %8.6f sec\n", T);

mprintf("\nPrewarped Analog Frequencies:\n");
mprintf("----------------------------------------------------\n");
mprintf("Omega_p                      = %10.4f rad/sec\n", omegap);
mprintf("Omega_s                      = %10.4f rad/sec\n", omegas);

mprintf("\nFilter Order:\n");
mprintf("----------------------------------------------------\n");
mprintf("Calculated Order (N)         = %8.4f\n", N_calc);
mprintf("Rounded Filter Order         = %d\n", N);

mprintf("\nDesign Parameters:\n");
mprintf("----------------------------------------------------\n");
mprintf("Ripple Factor (Epsilon)      = %8.6f\n", Epsilon);
mprintf("Cutoff Frequency (Omega_c)   = %10.4f rad/sec\n", omegac);
mprintf("Prototype Gain               = %10.4e\n", gn);

// Optional display of transfer functions
disp("Analog High Pass Transfer Function:");
disp(hs_hpf);

disp("Digital Transfer Function H(z):");
disp(Hz_sys);

// ================== FREQUENCY RESPONSE ==================
[HW, w] = frmag(Hz_sys, 512);

figure();
plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency (w/pi)');
ylabel('Magnitude');
title('Frequency Response of Chebyshev Type-I High Pass Filter');

```

## PROGRAM (BPF): 
```py
clc;
clear;
close;

// ======================================================
// CHEBYSHEV TYPE-I DIGITAL BAND PASS FILTER
// ======================================================

// ================== USER INPUT ==================
wp1    = input('Enter the lower pass band frequency (Radians)= ');
wp2    = input('Enter the upper pass band frequency (Radians)= ');

ws1    = input('Enter the lower stop band frequency (Radians)= ');
ws2    = input('Enter the upper stop band frequency (Radians)= ');

alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');

T      = input('Enter the value of sampling time= ');


// ======================================================
// PRE-WARPING
// ======================================================

omegap1 = (2/T)*tan(wp1/2);
omegap2 = (2/T)*tan(wp2/2);

omegas1 = (2/T)*tan(ws1/2);
omegas2 = (2/T)*tan(ws2/2);


// ======================================================
// BAND PASS PARAMETERS
// ======================================================

// Bandwidth
B = omegap2 - omegap1;

// Center frequency
omegac = sqrt(omegap1 * omegap2);


// ======================================================
// CONVERT BAND-PASS SPECIFICATION
// TO EQUIVALENT LOW-PASS SPECIFICATION
// ======================================================

Omega_s1 = abs((omegas1^2 - omegac^2) / (B*omegas1));

Omega_s2 = abs((omegas2^2 - omegac^2) / (B*omegas2));

// Take the smaller normalized stopband frequency
Omega_s = min(Omega_s1, Omega_s2);


// ======================================================
// FILTER ORDER
// ======================================================

Epsilon = sqrt((10^(0.1*alphap))-1);

N_calc = acosh( ...
    sqrt(((10^(0.1*alphas))-1) / ...
    ((10^(0.1*alphap))-1)) ...
    ) / acosh(Omega_s);

N = ceil(N_calc);


// ======================================================
// NORMALIZED ANALOG LOW PASS CHEBYSHEV PROTOTYPE
// ======================================================

// Normalized cutoff = 1 rad/sec
[pols, gn] = zpch1(N, Epsilon, 1);

s = poly(0, 's');

Hs_lp = syslin('c', ...
    gn, ...
    real(poly(pols, 's')));


// ======================================================
// LOW PASS --> BAND PASS TRANSFORMATION
// ======================================================

// s --> (s^2 + omegac^2)/(B*s)

Hs_bp = horner(Hs_lp, ...
    (s^2 + omegac^2)/(B*s));


// ======================================================
// BILINEAR TRANSFORMATION
// ======================================================

// s = (2/T) * ((1-z^-1)/(1+z^-1))

z = poly(0, 'z');

Hz = horner(Hs_bp, ...
    (2/T)*((1-z^-1)/(1+z^-1)));

Hz_sys = syslin('d', Hz);


// ======================================================
// NEAT CONSOLE OUTPUT
// ======================================================

mprintf("\n====================================================\n");
mprintf("      CHEBYSHEV TYPE-I DIGITAL BAND PASS FILTER\n");
mprintf("====================================================\n");


mprintf("\nInput Specifications:\n");
mprintf("----------------------------------------------------\n");

mprintf("Lower Passband Frequency (wp1) = %8.4f rad/sample\n", wp1);
mprintf("Upper Passband Frequency (wp2) = %8.4f rad/sample\n", wp2);

mprintf("Lower Stopband Frequency (ws1) = %8.4f rad/sample\n", ws1);
mprintf("Upper Stopband Frequency (ws2) = %8.4f rad/sample\n", ws2);

mprintf("Passband Attenuation (Ap)       = %8.2f dB\n", alphap);
mprintf("Stopband Attenuation (As)       = %8.2f dB\n", alphas);

mprintf("Sampling Time (T)               = %8.6f sec\n", T);


// ======================================================
// PRE-WARPED FREQUENCIES
// ======================================================

mprintf("\nPrewarped Analog Frequencies:\n");
mprintf("----------------------------------------------------\n");

mprintf("Omega_p1 = %10.4f rad/sec\n", omegap1);
mprintf("Omega_p2 = %10.4f rad/sec\n", omegap2);

mprintf("Omega_s1 = %10.4f rad/sec\n", omegas1);
mprintf("Omega_s2 = %10.4f rad/sec\n", omegas2);


// ======================================================
// BAND PASS PARAMETERS
// ======================================================

mprintf("\nBand Pass Parameters:\n");
mprintf("----------------------------------------------------\n");

mprintf("Bandwidth (B)                  = %10.4f rad/sec\n", B);

mprintf("Center Frequency (Omega_c)     = %10.4f rad/sec\n", omegac);

mprintf("Equivalent Stopband Omega_s1   = %10.4f\n", Omega_s1);
mprintf("Equivalent Stopband Omega_s2   = %10.4f\n", Omega_s2);

mprintf("Selected Omega_s               = %10.4f\n", Omega_s);


// ======================================================
// FILTER ORDER
// ======================================================

mprintf("\nFilter Order:\n");
mprintf("----------------------------------------------------\n");

mprintf("Calculated Order (N)           = %8.4f\n", N_calc);

mprintf("Rounded Filter Order            = %d\n", N);


// ======================================================
// DESIGN PARAMETERS
// ======================================================

mprintf("\nDesign Parameters:\n");
mprintf("----------------------------------------------------\n");

mprintf("Ripple Factor (Epsilon)         = %8.6f\n", Epsilon);

mprintf("Bandwidth (B)                   = %10.4f rad/sec\n", B);

mprintf("Center Frequency (Omega_c)      = %10.4f rad/sec\n", omegac);

mprintf("Prototype Gain                  = %10.4e\n", gn);


// ======================================================
// TRANSFER FUNCTIONS
// ======================================================

disp(" ");
disp("Analog Low Pass Prototype H(s):");
disp(Hs_lp);

disp(" ");
disp("Analog Band Pass Transfer Function H(s):");
disp(Hs_bp);

disp(" ");
disp("Digital Band Pass Transfer Function H(z):");
disp(Hz_sys);


// ======================================================
// FREQUENCY RESPONSE
// ======================================================

[HW, w] = frmag(Hz_sys, 1024);

figure();

plot(w/%pi, abs(HW));

xlabel('Normalized Digital Frequency (w/pi)');
ylabel('Magnitude');

title('Frequency Response of Chebyshev Type-I Band Pass Filter');

xgrid();

```

## PROGRAM (BSF): 
```py
clc;
clear;
close;

// ======================================================
// CHEBYSHEV TYPE-I DIGITAL BAND STOP FILTER
// ======================================================

// ================== USER INPUT ==================
wp1    = input('Enter the lower pass band frequency (Radians)= ');
ws1    = input('Enter the lower stop band frequency (Radians)= ');

ws2    = input('Enter the upper stop band frequency (Radians)= ');
wp2    = input('Enter the upper pass band frequency (Radians)= ');

alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');

T      = input('Enter the value of sampling time= ');


// ======================================================
// PRE-WARPING
// ======================================================

omegap1 = (2/T)*tan(wp1/2);
omegas1 = (2/T)*tan(ws1/2);

omegas2 = (2/T)*tan(ws2/2);
omegap2 = (2/T)*tan(wp2/2);


// ======================================================
// BAND STOP PARAMETERS
// ======================================================

// Bandwidth
B = omegas2 - omegas1;

// Center frequency
omegac = sqrt(omegas1 * omegas2);


// ======================================================
// CONVERT BAND-STOP SPECIFICATION
// TO EQUIVALENT LOW-PASS SPECIFICATION
// ======================================================

// For lower stopband edge
Omega_s1 = abs((B*omegas1) / ...
               (omegac^2 - omegas1^2));

// For upper stopband edge
Omega_s2 = abs((B*omegas2) / ...
               (omegas2^2 - omegac^2));


// Passband frequencies transformed to
// equivalent low-pass stopband frequencies

Omega_p1 = abs((B*omegap1) / ...
               (omegac^2 - omegap1^2));

Omega_p2 = abs((B*omegap2) / ...
               (omegap2^2 - omegac^2));

// Select the smallest equivalent passband frequency
Omega_p = min(Omega_p1, Omega_p2);


// ======================================================
// FILTER ORDER
// ======================================================

// Normalize equivalent passband frequency to 1

Omega_s = min(Omega_s1, Omega_s2);

Epsilon = sqrt((10^(0.1*alphap))-1);


// Order calculation
N_calc = acosh( ...
    sqrt(((10^(0.1*alphas))-1) / ...
    ((10^(0.1*alphap))-1)) ...
    ) / acosh(Omega_s/Omega_p);

N = ceil(N_calc);


// ======================================================
// NORMALIZED ANALOG LOW PASS CHEBYSHEV PROTOTYPE
// ======================================================

[pols, gn] = zpch1(N, Epsilon, 1);

s = poly(0, 's');

Hs_lp = syslin('c', ...
    gn, ...
    real(poly(pols, 's')));


// ======================================================
// LOW PASS --> BAND STOP TRANSFORMATION
// ======================================================

// s --> (B*s)/(s^2 + omegac^2)

Hs_bs = horner(Hs_lp, ...
    (B*s)/(s^2 + omegac^2));


// ======================================================
// BILINEAR TRANSFORMATION
// ======================================================

// s = (2/T)*((1-z^-1)/(1+z^-1))

z = poly(0, 'z');

Hz = horner(Hs_bs, ...
    (2/T)*((1-z^-1)/(1+z^-1)));

Hz_sys = syslin('d', Hz);


// ======================================================
// NEAT CONSOLE OUTPUT
// ======================================================

mprintf("\n====================================================\n");
mprintf("      CHEBYSHEV TYPE-I DIGITAL BAND STOP FILTER\n");
mprintf("====================================================\n");

mprintf("\nInput Specifications:\n");
mprintf("----------------------------------------------------\n");

mprintf("Lower Passband Frequency (wp1) = %8.4f rad/sample\n", wp1);
mprintf("Lower Stopband Frequency (ws1) = %8.4f rad/sample\n", ws1);
mprintf("Upper Stopband Frequency (ws2) = %8.4f rad/sample\n", ws2);
mprintf("Upper Passband Frequency (wp2) = %8.4f rad/sample\n", wp2);

mprintf("Passband Attenuation (Ap)       = %8.2f dB\n", alphap);
mprintf("Stopband Attenuation (As)       = %8.2f dB\n", alphas);

mprintf("Sampling Time (T)               = %8.6f sec\n", T);


// ======================================================
// PRE-WARPED FREQUENCIES
// ======================================================

mprintf("\nPrewarped Analog Frequencies:\n");
mprintf("----------------------------------------------------\n");

mprintf("Omega_p1 = %10.4f rad/sec\n", omegap1);
mprintf("Omega_s1 = %10.4f rad/sec\n", omegas1);
mprintf("Omega_s2 = %10.4f rad/sec\n", omegas2);
mprintf("Omega_p2 = %10.4f rad/sec\n", omegap2);


// ======================================================
// BAND STOP PARAMETERS
// ======================================================

mprintf("\nBand Stop Parameters:\n");
mprintf("----------------------------------------------------\n");

mprintf("Bandwidth (B)              = %10.4f rad/sec\n", B);

mprintf("Center Frequency (Omega_c) = %10.4f rad/sec\n", omegac);

mprintf("Equivalent Omega_s1        = %10.4f\n", Omega_s1);
mprintf("Equivalent Omega_s2        = %10.4f\n", Omega_s2);

mprintf("Equivalent Omega_p1        = %10.4f\n", Omega_p1);
mprintf("Equivalent Omega_p2        = %10.4f\n", Omega_p2);


// ======================================================
// FILTER ORDER
// ======================================================

mprintf("\nFilter Order:\n");
mprintf("----------------------------------------------------\n");

mprintf("Calculated Order (N)       = %8.4f\n", N_calc);
mprintf("Rounded Filter Order       = %d\n", N);


// ======================================================
// DESIGN PARAMETERS
// ======================================================

mprintf("\nDesign Parameters:\n");
mprintf("----------------------------------------------------\n");

mprintf("Ripple Factor (Epsilon)    = %8.6f\n", Epsilon);

mprintf("Bandwidth (B)              = %10.4f rad/sec\n", B);

mprintf("Center Frequency           = %10.4f rad/sec\n", omegac);

mprintf("Prototype Gain             = %10.4e\n", gn);


// ======================================================
// TRANSFER FUNCTIONS
// ======================================================

disp(" ");
disp("Analog Low Pass Prototype H(s):");
disp(Hs_lp);

disp(" ");
disp("Analog Band Stop Transfer Function H(s):");
disp(Hs_bs);

disp(" ");
disp("Digital Band Stop Transfer Function H(z):");
disp(Hz_sys);


// ======================================================
// FREQUENCY RESPONSE
// ======================================================

[HW, w] = frmag(Hz_sys, 1024);

figure();

plot(w/%pi, abs(HW));

xlabel('Normalized Digital Frequency (w/pi)');
ylabel('Magnitude');

title('Frequency Response of Chebyshev Type-I Band Stop Filter');

xgrid();

```

## OUTPUT (LPF) : 

<img width="501" height="481" alt="image" src="https://github.com/user-attachments/assets/69de13bf-90f1-4a44-b2c8-dc4254daab73" />

<img width="552" height="454" alt="image" src="https://github.com/user-attachments/assets/5c263a9c-da93-4c67-9307-7517d43bac26" />

## OUTPUT (HPF) : 

<img width="599" height="481" alt="image" src="https://github.com/user-attachments/assets/9e3c138f-7769-4cdd-bdc2-baac75eafefc" />

<img width="505" height="439" alt="image" src="https://github.com/user-attachments/assets/3c6c822a-8d38-42e9-a9bf-8b4c3b972575" />

## OUTPUT (BPF) : 

<img width="503" height="473" alt="image" src="https://github.com/user-attachments/assets/b544763f-eb7b-4427-a698-9812d15973bf" />

<img width="570" height="533" alt="image" src="https://github.com/user-attachments/assets/162463a0-b686-42b9-bd59-fa4ac15db47a" />

<img width="567" height="522" alt="image" src="https://github.com/user-attachments/assets/c21b712f-4f84-4bc6-960f-2823ffe0615a" />

## OUTPUT (BSF) : 

<img width="502" height="578" alt="image" src="https://github.com/user-attachments/assets/8544bbe9-4482-4cb0-9efa-f69350700344" />

<img width="600" height="557" alt="image" src="https://github.com/user-attachments/assets/29bec015-4b1b-4b90-a11e-aefc35fc95f6" />

<img width="663" height="499" alt="image" src="https://github.com/user-attachments/assets/354e7dd1-d986-4ecd-8c5f-eab245e31db2" />

## RESULT: 

Thus we design an IIR Chebyshev filter  using SCILAB. 
