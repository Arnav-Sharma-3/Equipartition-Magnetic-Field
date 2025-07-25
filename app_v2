import streamlit as st
import pandas as pd
import math
from math import *

# --------------------------------------------------
# Constants
# --------------------------------------------------
CGS_KPC  = 3.08567758128e21    # cm per kiloparsec
CGS_MPC  = 3.08567758128e24    # cm per Megaparsec
C1       = 6.266e18
C3       = 2.368e-3
M_E      = 9.1093837139e-28
C_LIGHT  = 2.99792458e10
X_FACTOR = 0.0

# --------------------------------------------------
# Cosmology Calculator (from redshift)
# --------------------------------------------------
def run_cosmology_calculator(z, H0, WM, WV):
    h = H0 / 100.
    WR = 4.165E-5 / (h * h)
    WK = 1 - WM - WR - WV
    az = 1.0 / (1.0 + z)
    Tyr = 977.8
    c = 299792.458

    n = 1000
    age = 0.0
    for i in range(n):
        a = az * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        age += 1.0 / adot
    zage = az * age / n

    DTT = 0.0
    DCMR = 0.0
    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
        DTT += 1.0 / adot
        DCMR += 1.0 / (a * adot)

    DTT = (1.0 - az) * DTT / n
    DCMR = (1.0 - az) * DCMR / n

    x = sqrt(abs(WK)) * DCMR
    if x > 0.1:
        if WK > 0:
            ratio = 0.5 * (exp(x) - exp(-x)) / x
        else:
            ratio = sin(x) / x
    else:
        y = x * x
        if WK < 0:
            y = -y
        ratio = 1. + y / 6. + y * y / 120.

    DCMT = ratio * DCMR
    DA = az * DCMT
    DL = DA / (az * az)

    DA_Mpc = (c / H0) * DA
    DL_Mpc = (c / H0) * DL

    return DL_Mpc, DA_Mpc, az  # Luminosity Distance, Angular Diameter, Scale Factor

# --------------------------------------------------
# Magnetic Field Calculator
# --------------------------------------------------
def compute_fields(alpha, g1, g2, v0, s_v0, l, b, w, D_l, Sf):
    l_cm = l * Sf * CGS_KPC
    b_cm = b * Sf * CGS_KPC
    w_cm = w * Sf * CGS_KPC
    D_l_cm = D_l * CGS_MPC
    v0_hz = v0 * 1e6
    s_v0_cgs = s_v0 * 1e-23

    p = 2 * alpha + 1
    V = (4 / 3) * math.pi * l_cm * b_cm * w_cm * 0.125
    L1 = 4 * math.pi * D_l_cm**2 * s_v0_cgs * v0_hz**alpha

    T3 = (g2 - 1)**(2 - p) - (g1 - 1)**(2 - p)
    T4 = (g2 - 1)**(2 * (1 - alpha)) - (g1 - 1)**(2 * (1 - alpha))
    T5 = (g2 - 1)**(3 - p) - (g1 - 1)**(3 - p)
    T6 = T3 * T4 / T5

    T1 = 3 * L1 / (2 * C3 * (M_E * C_LIGHT**2)**(2 * alpha - 1))
    T2 = (1 + X_FACTOR) / (1 - alpha) * (3 - p) / (2 - p) * (math.sqrt(2/3) * C1)**(1 - alpha)
    A = T1 * T2 * T6
    L = L1 / (1 - alpha) * (math.sqrt(2/3) * C1 * (M_E * C_LIGHT**2)**2)**(1 - alpha) * T4

    B_min = ((4 * math.pi * (1 + alpha) * A) / V)**(1 / (3 + alpha))
    B_eq = (2 / (1 + alpha))**(1 / (3 + alpha)) * B_min

    u_b = B_min**2 / (8 * math.pi)
    u_p = A / V * B_min**(-1 + alpha)
    u_tot = u_p + u_b

    return alpha, B_min * 1e6, B_eq * 1e6, D_l, L, u_p, u_b, u_tot

# --------------------------------------------------
# Streamlit App
# --------------------------------------------------
st.set_page_config(page_title="CosmoMagnetic Field Calculator", layout="wide")
st.title("🌀 Magnetic Field Calculator with Redshift-based Cosmology")

st.markdown(
    """
    Upload a CSV with columns:  
    `Source, z, alpha, gamma1, gamma2, v0, s_v0, l, b, w`  
    — where **l, b, w** are in **kpc**, **z** is redshift, **v0** in **MHz**, **s_v0** in **Jy**.
    """
)

# Cosmology settings
st.sidebar.header("Cosmology Parameters")
H0 = st.sidebar.number_input("Hubble Constant (H₀)", value=70.0)
WM = st.sidebar.slider("Ω Matter (Ωₘ)", 0.0, 1.5, 0.3)
WV = st.sidebar.slider("Ω Vacuum (ΩΛ)", 0.0, 1.5, 0.7)

uploaded_file = st.file_uploader("Upload your data file", type=["csv", "tsv", "txt"])
if uploaded_file:
    sep = "\t" if uploaded_file.name.endswith((".tsv", ".txt")) else ","
    try:
        df = pd.read_csv(uploaded_file, sep=sep, comment="#")
    except Exception as e:
        st.error(f"❌ Could not read file: {e}")
    else:
        required = ["Source", "z", "alpha", "gamma1", "gamma2", "v0", "s_v0", "l", "b", "w"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            st.error(f"❌ Missing columns: {', '.join(missing)}")
        else:
            # Compute cosmology-based values and magnetic fields
            results = []
            for _, row in df.iterrows():
                DL, DA, Sf = run_cosmology_calculator(row["z"], H0, WM, WV)
                output = compute_fields(
                    row["alpha"], row["gamma1"], row["gamma2"],
                    row["v0"], row["s_v0"],
                    row["l"], row["b"], row["w"],
                    DL, Sf
                )
                results.append({
                    "Source": row["Source"],
                    "z": row["z"],
                    "Alpha": output[0],
                    "B_min (μG)": output[1],
                    "B_eq (μG)": output[2],
                    "D_L (Mpc)": f"{output[3]:.2f}",
                    "L (erg/s)": f"{output[4]:.2e}",
                    "u_p (erg/cm³)": f"{output[5]:.2e}",
                    "u_B (erg/cm³)": f"{output[6]:.2e}",
                    "u_total (erg/cm³)": f"{output[7]:.2e}",
                    "D_A (Mpc)": f"{DA:.2f}",
                    "Scale Factor (a)": f"{Sf:.4f}"
                })

            df_out = pd.DataFrame(results)
            st.success("✅ Calculation complete!")
            st.dataframe(df_out)

            csv_data = df_out.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📥 Download Results as CSV",
                data=csv_data,
                file_name="magnetic_field_with_cosmology.csv",
                mime="text/csv"
            )

st.markdown(
    """
    <hr style="margin-top: 3rem;">
    <div style='text-align: center; font-size: 0.9rem; color: gray;'>
        Created by <b>Arnav Sharma</b><br>
        Under the Guidance of <b>Dr. Chiranjib Konar</b>
    </div>
    """,
    unsafe_allow_html=True
)
