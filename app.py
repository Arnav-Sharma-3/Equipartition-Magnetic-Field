import streamlit as st
import pandas as pd
import math

# --------------------------------------------------
# Constants for CGS conversions and synchrotron math
# --------------------------------------------------
CGS_KPC  = 3.08567758128e21    # cm per kiloparsec
CGS_MPC  = 3.08567758128e24    # cm per Megaparsec
C1       = 6.266e18            # synchrotron constant
C3       = 2.368e-3            # synchrotron constant
M_E      = 9.1093837139e-28    # electron mass (g)
C_LIGHT  = 2.99792458e10       # speed of light (cm/s)
X_FACTOR = 0.0                 # proton/electron energy ratio

def compute_fields(alpha, g1, g2, v0, s_v0, l, b, w, D_l, Sf, x=X_FACTOR):
    """
    Compute B_min and B_eq in µG for one row of input.
      - l, b, w should be in kpc
      - D_l in Mpc
      - v0 in MHz
      - s_v0 in Jy
    """
    # Convert kpc → cm, Mpc → cm
    l   *= Sf * CGS_KPC
    b   *= Sf * CGS_KPC
    w   *= Sf * CGS_KPC
    D_l *= CGS_MPC

    # Convert MHz → Hz, Jy → erg/s/cm²/Hz
    v0   *= 1e6
    s_v0 *= 1e-23

    # Derived quantities
    p   = 2 * alpha + 1
    V   = (4/3) * math.pi * l * b * w * 0.125
    L1  = 4 * math.pi * D_l**2 * s_v0 * v0**alpha

    # Synchrotron integration terms
    T3 = (g2 - 1)**(2 - p) - (g1 - 1)**(2 - p)
    T4 = (g2 - 1)**(2 * (1 - alpha)) - (g1 - 1)**(2 * (1 - alpha))
    T5 = (g2 - 1)**(3 - p) - (g1 - 1)**(3 - p)
    T6 = T3 * T4 / T5

    # Field-strength prefactors
    T1 = 3 * L1 / (2 * C3 * (M_E * C_LIGHT**2)**(2 * alpha - 1))
    T2 = (1 + x) / (1 - alpha) * (3 - p) / (2 - p) * (math.sqrt(2/3) * C1)**(1 - alpha)
    A  = T1 * T2 * T6
    L  = L1/ (1 - alpha) * (math.sqrt(2/3) * C1 * (M_E * C_LIGHT**2)**2 )**(1 - alpha) * T4
    
    # Compute B_min, B_eq in Gauss, then convert to µG
    B_min = ((4 * math.pi * (1 + alpha) * A) / V)**(1 / (3 + alpha))
    B_eq  = (2 / (1 + alpha))**(1 / (3 + alpha)) * B_min
    return B_min * 1e6, B_eq * 1e6, L.2f # µG

# -----------------------
# Streamlit App Layout
# -----------------------
st.set_page_config(page_title="Galaxy Magnetic Field Calculator", layout="centered")
st.title("🌀 Lobe Magnetic Field Estimator")
st.markdown(
    """
    Upload a CSV/TSV with columns:  
    `Source, alpha, gamma1, gamma2, v0, s_v0, l, b, w, D_l, Sf`  
    — where **l, b, w** are in **kpc**, **D_l** in **Mpc**, **v0** in **MHz**, **s_v0** in **Jy**.  
    The app will compute **B_min** & **B_eq** in µG for each row.
    """
)

uploaded_file = st.file_uploader("Upload your data file", type=["csv", "tsv", "txt"])
if uploaded_file:
    # Auto-detect delimiter
    sep = "\t" if uploaded_file.name.endswith((".tsv", ".txt")) else ","
    try:
        df = pd.read_csv(uploaded_file, sep=sep, comment="#")
    except Exception as e:
        st.error(f"📂 Could not read file: {e}")
    else:
        # Validate columns
        required = ["Source","alpha","gamma1","gamma2","v0","s_v0","l","b","w","D_l","Sf"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            st.error(f"❌ Missing columns: {', '.join(missing)}")
        else:
            # Vectorized compute
            results = df.apply(
                lambda r: compute_fields(
                    r["alpha"], r["gamma1"], r["gamma2"],
                    r["v0"],    r["s_v0"],
                    r["l"],     r["b"],     r["w"],
                    r["D_l"],   r["Sf"]
                ), axis=1, result_type="expand"
            )
            df["B_min (µG)"] = results[0].round(3)
            df["B_eq (µG)"]  = results[1].round(3)
            df["L (W/m^2)"]  = results[2].round(3)
            
            # Only show the three columns
            df_out = df[["Source", "B_min (µG)", "B_eq (µG)", "L (W/m^2)"]]
            st.success("✅ Calculation complete!")
            st.dataframe(df_out)

            # CSV download button
            csv_data = df_out.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📥 Download Results (CSV)",
                data=csv_data,
                file_name="magnetic_fields_results.csv",
                mime="text/csv"
            )

# ✅ Credits: Always visible at the bottom
st.markdown(
    """
    <hr style="margin-top: 3rem; margin-bottom: 1rem;">
    <div style='text-align: center; font-size: 0.9rem; color: gray;'>
        Created by <b>Arnav Sharma</b><br>
        Under the Guidance of <b>Dr. Chiranjib Konar</b>
    </div>
    """,
    unsafe_allow_html=True
)
