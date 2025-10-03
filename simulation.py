from time import time
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import math

U_sup = 3 # V
L = 0.1e-3 # H
R = 1.2 # Ohm
#k_e = 0.02 # Vs/rad
#k_t = 0.02 # Nm/A
k_e = 2e-3
k_t = 2e-3
T_abc = 0.0046720125 # Nm
#B = 0.5e-5 # N*m/(rad/s)
B = 0.5e-7 # N*m/(rad/s)
#J = 0.0002 # kg*m^2
J = 3.8 * 1e-8 # kg*m^2

r = 0.015 # m

# przekładnia
n_1 = 1 # liczba zębów na wale napędzającym
n_2 = 9 # liczba zębów na wale napędzanym
i_g = n_1 / n_2
eta = 0.85 # sprawność

T_abc = T_abc * (i_g * eta) if (i_g != 1) else T_abc

total_time = 0.3 # s
dt = max(L / R * 1e-2, 1e-6) # s

I = [0]
alpha = [0]
omega = [0]
omega_RPM = [0]
omega_load = [0]
omega_RPM_load = [0]
theta = [0]
theta_load = [0]
U = [U_sup]
U_bEMF = [0]
T_m = [0]
T_out = [0]
s = 0

f = 2

if __name__ == '__main__':
    steps = int(total_time / dt)
    start_time = time()
    for i in tqdm(range(1, steps)):
        U_bEMF.append(k_e * omega[i - 1])
        I.append(I[i-1] + (dt * ((U[i-1] - R * I[i-1] - U_bEMF[i])/L)))
        U.append(U_sup)
        omega.append(omega[i-1] + (dt * ((k_t * I[i] - B * omega[i-1] - T_abc) / J)))
        alpha.append((omega[i] - omega[i-1])/dt)
        omega_RPM.append(omega[i] * 60 / (2 * math.pi))
        omega_load.append(omega[i] * i_g)
        omega_RPM_load.append(omega_RPM[i] * i_g)
        theta.append(theta[i-1] + omega[i] * dt)
        theta_load.append(theta_load[i-1] + omega_load[i] * dt)
        T_m.append(k_t * I[i])
        T_out.append(T_m[i] / i_g * eta)
        s = s + r * (theta_load[i] - theta_load[i-1])
    stop_time = time()  
    
    print(f'Simulation is done! It took {stop_time - start_time:.2f}s.')
    print('Last values:')
    print(f'I_n = {I[steps-1]} A')
    print(f'U_n = {U[steps-1]} V')
    print(f'U_bEMF = {U_bEMF[steps-1]} V')
    print(f'T_m = {T_m[steps-1]} Nm')
    print(f'alpha_n = {alpha[steps-1]} rad/s^2')
    print(f'omega_n = {omega[steps-1]} rad/s = {omega_RPM[steps-1]} RPM')
    print(f'theta_n = {theta[steps-1]} rad')
    print(f's = {s} m')
    
    pd.DataFrame.from_dict({
        'I': I,
        'U': U,
        'U_bEMF': U_bEMF,
        'T_m': T_m,
        'alpha': alpha,
        'omega': omega,
        'theta': theta
    }).to_csv('result.csv')
    
    # po zakończonej symulacji (tam gdzie masz print)
    time_axis = [i * dt for i in range(len(I))]

    plt.figure(figsize=(12,8))

    plt.subplot(3,3,1)
    plt.plot(time_axis, I)
    plt.title("Prąd I(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("I [A]")

    plt.subplot(3,3,2)
    plt.plot(time_axis, omega_RPM)
    plt.title("Prędkość kątowa ω(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("ω [RPM]")

    plt.subplot(3,3,3)
    plt.plot(time_axis, theta)
    plt.title("Kąt θ(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("θ [rad]")

    plt.subplot(3,3,4)
    plt.plot(time_axis, alpha)
    plt.title("Przyspieszenie kątowe α(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("α [rad/s²]")
    
    plt.subplot(3,3,5)
    plt.plot(time_axis, U_bEMF)
    plt.title("Back-EMF U_bEMF(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("U_bEMF [V]")
    
    plt.subplot(3,3,6)
    plt.plot(time_axis, U)
    plt.title("Napięcie zasilania U(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("U [V]")
    
    plt.subplot(3,3,7)
    plt.plot(time_axis, T_m)
    plt.title("Moment obrotowy T_m(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("T_m [Nm]")
    
    plt.subplot(3,3,8)
    plt.plot(time_axis, omega_RPM_load)
    plt.title("Wyjściowa prędkość kątowa ω(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("ω [RPM]")
    
    plt.subplot(3,3,9)
    plt.plot(time_axis, T_out)
    plt.title("Wyjściowy moment obrotowy T_out(t)")
    plt.xlabel("Czas [s]")
    plt.ylabel("T_m [Nm]")

    plt.tight_layout()
    plt.show()