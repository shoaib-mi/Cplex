import math, cplex

# phi = phi_p - phi_c

const 
    beta = 1.0
    step = 0.001 / beta
    tmin = 0.0 / beta
    tmax = 5.0 / beta
    eta = 1.0
    delta_10 = -0.5 * beta
    delta_20 = 0.5 * beta
    omega = 5 * beta
    omega_21 = 1 * beta #5 * 10^14
    phi_c = Pi/2
    phi_p = 0.0
    theta = Pi/4
    length = int( (tmax - tmin) / step ) + 1
    I = complex(0.0,1.0)
    alpha = exp(-I * Pi/4) * sqrt(beta^3) / sqrt(Pi)

var
    h1, h2:  array[length,Complex]
    t:  array[length,float]
    sum_1 = 0.complex
    sum_2 = 0.complex
    f = open("data.csv",fmWrite)
    
t[0] = tmin
h1[0] = sin(theta) * exp(I * phi_p)
h2[0] = cos(theta).complex
f.writeLine(t[0] * beta,",", abs(h1[0])^2,",",abs(h2[0])^2)

for j in 0..length-2:
    sum_1 = 0.complex
    sum_2 = 0.complex
    t[j+1] = t[j] + step
    for i in 0..j-1:
        sum_1 = sum_1 + sqrt(t[j] - t[i]) * ( h1[i+1] - h1[i] + eta * h2[i+1] - eta * h2[i] )
        sum_2 = sum_2 + sqrt(t[j] - t[i]) * ( h2[i+1] - h2[i] + eta * h1[i+1] - eta * h1[i] )
    
    h1[j+1] = h1[j] + step * ( - I * delta_10 * h1[j] + omega * exp( I * (omega_21 * t[j] + phi_c)) * h2[j] - 2.0 * alpha * sqrt(t[j]) * (h1[0] + eta * h2[0]) - 2.0 * alpha * sum_1 )

    h2[j+1] = h2[j] + step * ( - I * delta_20 * h2[j] - omega * exp( - I * (omega_21 * t[j] + phi_c)) * h1[j] - 2.0 * alpha * sqrt(t[j]) * (h2[0] + eta * h1[0]) - 2.0 * alpha * sum_2 )

    f.writeLine(t[j+1] * beta,",",abs(h1[j+1])^2,",",abs(h2[j+1])^2)
    echo "%",round(j/length*100,2)
    
f.close()