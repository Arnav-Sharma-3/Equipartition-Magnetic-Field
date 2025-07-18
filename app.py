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
    # Convert kpc → cm, Mpc → cm
    l_cm = l * Sf * CGS_KPC
    b_cm = b * Sf * CGS_KPC
    w_cm = w * Sf * CGS_KPC
    D_l_cm = D_l * CGS_MPC

    # Convert MHz → Hz, Jy → erg/s/cm²/Hz
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
    T2 = (1 + x) / (1 - alpha) * (3 - p) / (2 - p) * (math.sqrt(2/3) * C1)**(1 - alpha)
    A = T1 * T2 * T6
    L = L1 / (1 - alpha) * (math.sqrt(2/3) * C1 * (M_E * C_LIGHT**2)**2)**(1 - alpha) * T4

    B_min = ((4 * math.pi * (1 + alpha) * A) / V)**(1 / (3 + alpha))
    B_eq = (2 / (1 + alpha))**(1 / (3 + alpha)) * B_min

    u_b = B_min**2 / (8 * math.pi)
    u_p = A / V * B_min**(-1 + alpha)
    u_tot = u_p + u_b

    return alpha, B_min * 1e6, B_eq * 1e6, D_l_cm, L, u_p, u_b, u_tot


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
    The app will compute **B_min**, **B_eq**, **Energy Densities**, and **Luminosity** in CGS units.
    """
)

uploaded_file = st.file_uploader("Upload your data file", type=["csv", "tsv", "txt"])
if uploaded_file:
    sep = "\t" if uploaded_file.name.endswith((".tsv", ".txt")) else ","
    try:
        df = pd.read_csv(uploaded_file, sep=sep, comment="#")
    except Exception as e:
        st.error(f"📂 Could not read file: {e}")
    else:
        required = ["Source","alpha","gamma1","gamma2","v0","s_v0","l","b","w","D_l","Sf"]
        missing = [c for c in required if c not in df.columns]
        if missing:
            st.error(f"❌ Missing columns: {', '.join(missing)}")
        else:
            results = df.apply(
                lambda r: compute_fields(
                    r["alpha"], r["gamma1"], r["gamma2"],
                    r["v0"],    r["s_v0"],
                    r["l"],     r["b"],     r["w"],
                    r["D_l"],   r["Sf"]
                ), axis=1, result_type="expand"
            )
            df_out = pd.DataFrame({
                "Source": df["Source"],
                "Alpha": results[0],
                "B_min (\u00b5G)": results[1].round(3),
                "B_eq (\u00b5G)": results[2].round(3),
                "D_L (cm)": results[3].apply(lambda x: f"{x:.2e}"),
                "L (erg/s)": results[4].apply(lambda x: f"{x:.2e}"),
                "u_p (erg/cm³)": results[5].apply(lambda x: f"{x:.2e}"),
                "u_B (erg/cm³)": results[6].apply(lambda x: f"{x:.2e}"),
                "u_total (erg/cm³)": results[7].apply(lambda x: f"{x:.2e}")
            })

            st.success("✅ Calculation complete!")
            st.dataframe(df_out)

            csv_data = df_out.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="📅 Download Results (CSV)",
                data=csv_data,
                file_name="magnetic_fields_results.csv",
                mime="text/csv"
            )

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
