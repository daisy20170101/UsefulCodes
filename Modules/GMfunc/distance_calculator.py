"""
Rupture distance calculations for ground motion models.

This module provides implementations for calculating various rupture distance metrics used in seismic hazard analysis and ground
motion prediction equations.

Author: Sanjay Bora
Version: 0.1.0
Last Modified: 01.10.2025
Organization: GNS Science

References:
    Kaklamanos, J., Baise, L.G., and Boore, D.M. (2011).
    "Estimating Unknown Input Parameters when Implementing the NGA
    Ground-Motion Prediction Equations in Engineering Practice".
    Earthquake Spectra, 27(4), 1219-1235.
"""

import numpy as np


class RuptureDistanceCalculator:
    """
    Corrected class for calculating rupture distances
    """

    def __init__(
        self,
        magnitude: float,
        ztor: float,
        mechanism: str = "SS",
        dip: float = 90.0,
        width: float = None,
        hypo_depth: float = None,
    ):

        self.magnitude = magnitude
        self.mechanism = mechanism.upper()
        self.ztor = ztor
        self.dip = dip

        if width is None:
            self.width = self._calculate_width()
        else:
            self.width = width

        if hypo_depth is None:
            self.hypo_depth = self.ztor + 0.4 * self.width * np.sin(
                np.radians(self.dip)
            )  # See Kaklamanos Paper equation 3 # I have changed the 0.6 to 0.4.
        else:
            self.hypo_depth = hypo_depth

        self._validate_inputs()

    def _validate_inputs(self):
        if self.mechanism not in ["SS", "NM", "RV"]:
            raise ValueError(f"Invalid mechanism: {self.mechanism}")
        if not 0 <= self.dip <= 90:
            raise ValueError(f"Dip must be between 0 and 90 degrees")
        if self.ztor < 0:
            raise ValueError(f"Ztor must be non-negative")

    def _calculate_width(self) -> float:
        """
        Calculate down-dip rupture width from magnitude using empirical scaling.

        Estimates the rupture width using magnitude-dependent relationships
        from Wells & Coppersmith (1994). Different coefficients are applied
        based on the fault mechanism type, as rupture dimensions vary
        systematically with faulting style.

        Returns:
            float: Estimated down-dip rupture width in kilometers.

        References:
            Wells, D.L. and Coppersmith, K.J. (1994).
            "New Empirical Relationships among Magnitude, Rupture Length,
            Rupture Width, Rupture Area, and Surface Displacement".
            Bulletin of the Seismological Society of America, 84(4), 974-1002.
            https://doi.org/10.1785/BSSA0840040974"""
        m = self.magnitude

        if self.mechanism == "SS":
            width = 10 ** (-1.0 + 0.3 * m)
        elif self.mechanism == "NM":
            width = 10 ** (-1.14 + 0.35 * m)
        elif self.mechanism == "RV":
            width = 10 ** (-1.61 + 0.41 * m)
        else:
            width = 10 ** (-1.01 + 0.32 * m)

        return width

    def get_rjb_rrup(self, rx):
        """
        Calculate Joyner-Boore (Rjb) and rupture (Rrup) distances.

        Computes two fundamental distance metrics used in ground motion prediction
        equations based on the horizontal distance (rx) from the surface projection
        of the top edge of the rupture plane.

        This is a corrected implementation that properly handles dipping faults
        and follows the Kaklamanos et al. (2011) methodology, with verification
        against Peter Stafford's reference implementation.

        Args:
            rx: Horizontal distance(s) perpendicular to fault strike in km.
                Can be scalar or array. Negative values indicate hanging wall side
                (updip from surface projection), positive values indicate footwall
                side (downdip).

        Returns:
            tuple: A tuple containing:
                - rrup (np.ndarray): Closest distance to rupture plane in km.
                - rjb (np.ndarray): Joyner-Boore distance (closest horizontal
                  distance to surface projection of rupture) in km.
        Notes:
        Distance Definitions:
            - **Rjb (Joyner-Boore)**: Shortest horizontal distance from site
                to surface projection of rupture plane. Always >= 0.
            - **Rrup**: Shortest 3D distance from site to rupture plane.

        Calculation Method:
            For vertical faults (dip = 90°):
                Rrup = √(Rjb² + Ztor²)

                For dipping faults, three regions are considered:
                    1. Updip (rx < x₁): Site above top of rupture
                       Rrup = √(rx² + Ztor²)

                    2. Above fault plane (x₁ ≤ rx ≤ x₂): Site between surface
                       projections of top and bottom of rupture
                       Rrup = rx·sin(δ) + Ztor·cos(δ)

                    3. Downdip (rx > x₂): Site beyond bottom of rupture
                       Rrup = √((rx - W·cos(δ))² + (Ztor + W·sin(δ))²)

                Where:
                    x₁ = Ztor·tan(δ)  (transition to fault plane)
                    x₂ = x₁ + W/cos(δ)  (transition to downdip)
                    W = rupture width, δ = dip angle, Ztor = depth to top

            Corrections Applied:
                - Changed Rrup calculation in "above fault plane" region from
                  (rx - x₁)·sin(δ) to rx·sin(δ) based on Kaklamanos paper
                  and verification with Peter Stafford's code
                - Proper handling of Rjb for negative rx values (hanging wall)

                Raises:
            ValueError: If rx contains NaN or infinite values (implicit via numpy).

        See Also:
            get_rhyp : Calculate hypocentral distance
            _validate_inputs : Ensures fault geometry is valid

        References:
            .. [1] Kaklamanos, J., Baise, L.G., and Boore, D.M. (2011).
                   "Estimating Unknown Input Parameters when Implementing the NGA
                   Ground-Motion Prediction Equations in Engineering Practice".
                   Earthquake Spectra, 27(4), 1219-1235.
            .. [2] Joyner, W.B. and Boore, D.M. (1981).
                   "Peak horizontal acceleration and velocity from strong-motion
                   records including records from the 1979 Imperial Valley,
                   California, earthquake". Bull. Seismol. Soc. Am., 71(6), 2011-2038.

        Warnings:
            - Results are sensitive to fault dip angle for sites near the rupture
            - For very shallow ruptures (Ztor ≈ 0), numerical precision may be limited
            - Distances are computed assuming a planar rupture surface
        """
        rx = np.atleast_1d(rx)

        # Convert dip to radians
        dip_rad = np.radians(self.dip)

        # Calculate projections
        w_cos = self.width * np.cos(dip_rad)  # Horizontal width
        w_sin = self.width * np.sin(dip_rad)  # Vertical component

        # Initialize arrays
        rjb = np.zeros_like(rx, dtype=float)
        rrup = np.zeros_like(rx, dtype=float)

        # Calculate Rjb - this should be correct
        rjb = np.where(rx < 0, np.abs(rx), np.maximum(0.0, rx - w_cos))

        # Calculate Rrup
        if self.dip == 90.0:  # Vertical fault
            # For vertical fault, Rrup is hypotenuse of Rjb and Ztor
            rrup = np.sqrt(rjb**2 + self.ztor**2)
        else:  # Dipping fault
            # Calculate transition points
            x1 = self.ztor * np.tan(dip_rad) if self.ztor > 0 else 0
            x2 = x1 + self.width / np.cos(dip_rad)

            # Three regions for dipping fault
            for i, x in enumerate(rx):
                if x < x1:  # Updip region
                    rrup[i] = np.sqrt(x**2 + self.ztor**2)
                elif x <= x2:  # Above fault plane
                    rrup[i] = (x) * np.sin(dip_rad) + self.ztor * np.cos(
                        dip_rad
                    )  # AI code gave (x-x1). I changed it to the paper and checked Peter's code.
                else:  # Downdip region
                    rrup[i] = np.sqrt((x - w_cos) ** 2 + (self.ztor + w_sin) ** 2)

        return np.squeeze(rrup), np.squeeze(rjb)

    def calculate_ry(self, rx, rjb, alpha: float = 90.0):

        rx = np.atleast_1d(rx)
        rjb = np.atleast_1d(rjb)

        azimuth_rad = np.radians(alpha)

        if alpha in [0.0, 180.0, -180.0]:
            ry = rjb
        elif alpha == 90.0 or alpha == -90.0:
            ry = np.zeros_like(rjb)
        else:
            ry = np.abs(rx / np.tan(azimuth_rad))
        return np.squeeze(ry)

    def calculate_rhypo(self, rx: float, ry: float):
        """
        Calculate hypocentral distance. This code is not checked.

        Parameters
        ----------
        rx : float or array-like
            Distance perpendicular to strike
        ry : float or array-like
            Distance along strike from hypocenter

        Returns
        -------
        rhypo : array
            Hypocentral distance
        """
        rx = np.atleast_1d(rx)
        ry = np.atleast_1d(ry)

        # Get Rrup first
        rrup, _ = self.get_rjb_rrup(rx)

        # Simple approximation - more complex calculation would consider
        # actual hypocenter location on fault plane
        rhypo = np.sqrt(rrup**2 + ry**2)

        return rhypo

    def check_calculations(self, rx_test=None):
        """
        Method to verify calculations are working
        """
        if rx_test is None:
            rx_test = np.array([0, 10, 20, 30])

        rrup, rjb = self.get_rjb_rrup(rx_test)

        print(f"\nFault Parameters:")
        print(f"  Magnitude: {self.magnitude}")
        print(f"  Mechanism: {self.mechanism}")
        print(f"  Dip: {self.dip}°")
        print(f"  Ztor: {self.ztor} km")
        print(f"  Width: {self.width:.2f} km")
        print(f"  Horizontal width: {self.width * np.cos(np.radians(self.dip)):.2f} km")

        print(f"\nDistance Calculations:")
        for i in range(len(rx_test)):
            print(
                f"  Rx={rx_test[i]:6.1f} km -> Rjb={rjb[i]:6.2f} km, Rrup={rrup[i]:6.2f} km"
            )

        # Check if all values are the same (which would indicate a problem)
        if np.allclose(rx_test, rjb) and np.allclose(rjb, rrup):
            print("\n⚠️ WARNING: All values are identical - there may be an issue!")
        else:
            print("\n✓ Values are different as expected")


# Simple test function
def simple_test():
    """Run a simple test to verify calculations"""
    print("SIMPLE TEST - M7.0 Strike-Slip")
    print("=" * 60)

    # Create calculator
    calc = RuptureDistanceCalculator(7.0, "SS", ztor=0, dip=90)

    # Test at specific distances
    rx_values = [0, 10, 20]
    rrup, rjb = calc.get_rjb_rrup(rx_values)

    print(f"Fault width: {calc.width:.2f} km")
    print(f"\nResults:")
    print(f"{'Rx (input)':>12} | {'Rjb':>10} | {'Rrup':>10}")
    print("-" * 40)
    for i in range(len(rx_values)):
        print(f"{rx_values[i]:>12.1f} | {rjb[i]:>10.2f} | {rrup[i]:>10.2f}")

    # Expected values for verification
    print(f"\nExpected Rjb values:")
    print(f"  At Rx=0: should be 0")
    print(f"  At Rx=10: should be 0 if width > 10, else Rx-width")
    print(f"  At Rx=20: should be max(0, 20-width)")

    print(f"\nFor vertical fault with Ztor=0:")
    print(f"  Rrup should equal Rjb when Ztor=0")
    print(f"  With Ztor>0, Rrup = sqrt(Rjb² + Ztor²)")
