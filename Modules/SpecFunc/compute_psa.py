import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import lfilter

def compute_psa(acc, dt, period=1.0, damping=0.05):
    """
    Compute the 5%-damped pseudo spectral acceleration (PSA) at a given period.
    
    Parameters:
    acc : array_like
        Ground acceleration time series (m/s²)
    dt : float
        Time step (s)
    period : float
        Oscillator period (s), default is 1.0
    damping : float
        Damping ratio, default is 0.05 (5%)
    
    Returns:
    psa : float
        Pseudo spectral acceleration (m/s²)
    """
    omega = 2 * np.pi / period
    wn2 = omega**2
    wd = omega * np.sqrt(1 - damping**2)

    # Newmark-beta method coefficients
    beta = 0.25
    gamma = 0.5
    npts = len(acc)

    u = np.zeros(npts)     # displacement
    v = np.zeros(npts)     # velocity
    a = np.zeros(npts)     # acceleration (relative)

    # effective stiffness and constants
    k_eff = wn2 + 2*damping*omega*gamma/dt + 1/(beta*dt**2)
    a1 = 1/(beta*dt)
    a2 = 1/(2*beta) - 1

    for i in range(1, npts):
        dp = - (acc[i] - acc[i-1])  # relative change in force (unit mass)
        du = (dp + a1*v[i-1] + a2*a[i-1]) / k_eff
        dv = gamma * du / (beta*dt) - gamma*v[i-1]/beta + dt*(1 - gamma/(2*beta))*a[i-1]
        da = (du - dt*v[i-1] - dt**2*a[i-1]/2) / (beta*dt**2)

        u[i] = u[i-1] + du
        v[i] = v[i-1] + dv
        a[i] = a[i-1] + da

    psa = omega**2 * np.max(np.abs(u))  # pseudo spectral acceleration
    return psa

