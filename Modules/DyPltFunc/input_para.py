import numpy as np
import matplotlib.pyplot as plt


def plt_omega_func(z1, z2, z3, z4, appendix="test1"):

    ## plot omega (coefficient of deviatory stress vs. confining stress) as a function of depth

    zz = -np.linspace(0.1, 20000, 40)
    zStressDecreaseStart = z3
    zStressDecreaseStop = z4
    zStressIncrease = z2
    zStressIncreaseWidth = z1 - z2
    zStressDecreaseWidth = zStressDecreaseStart - zStressDecreaseStop

    r1 = np.linspace(0, 0, 40)

    for ik, z in enumerate(zz):
        # print(z,ik)

        if z >= zStressIncrease:
            a = (z - zStressIncrease) / zStressIncreaseWidth
            Sx = 3.0 * a * a - 2.0 * a * a * a
            r1[ik] = 1.0 - Sx
        else:
            if z >= zStressDecreaseStart:
                r1[ik] = 1.0
            else:
                if z >= zStressDecreaseStop:
                    a = 1.0 - (z - zStressDecreaseStop) / zStressDecreaseWidth
                    Sx = 3.0 * a * a - 2.0 * a * a * a
                    r1[ik] = 1.0 - Sx
                else:
                    r1[ik] = 0.001

    fig, axe = plt.subplots(1, 2, figsize=(5, 5))
    axe[0].plot(r1, zz, "-k")
    axe[0].set_xlabel("Omega")
    axe[0].set_ylabel("depth")
    plt.savefig("omega-" + appendix + ".png")
    plt.show()

    return fig, axe
