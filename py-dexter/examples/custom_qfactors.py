import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

# ==========================================

psit_values = np.linspace(0, 0.4, 5)
q_array = 2 + 25 * psit_values**2
qfactor = dex.ArrayQfactor(
    flux_array=dex.MagneticFluxArray.Toroidal(psit_values),
    q_array=q_array,
    interp_type="Cubic",
)

dex.plot_qfactor(qfactor, data=True)

# ==========================================

psit_values = np.linspace(0, 0.4, 50)
q_array = 1.3 + 15 * psit_values**4 - 2 * psit_values**3
qfactor = dex.ArrayQfactor(
    flux_array=dex.MagneticFluxArray.Toroidal(psit_values),
    q_array=q_array,
    interp_type="Cubic",
)

dex.plot_qfactor(qfactor, data=True)

# ==========================================

# RFP case
psip_values = [0, 0.1, 0.2, 0.3, 0.4]
q_array = [0.9, 0.7, 0.5, 0.1, -0.2]
qfactor = dex.ArrayQfactor(
    flux_array=dex.MagneticFluxArray.Poloidal(psip_values),
    q_array=q_array,
    interp_type="Cubic",
)

fig, ax = dex.plot_qfactor(qfactor, data=True, show=False)
ax[1].axhline(0, c="k", linestyle="--", label="$q=0$")
ax[1].legend()
plt.show()
