import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="NGS Yield Calculator", page_icon="🧬", layout="wide")

st.title("🧬 TapeStation & Qubit Mass Yield Calculator")
st.write("Upload your data logs to calculate the total regional mass ($ng$) present in a $50\\mu L$ sample tube volume.")

# Constants
TOTAL_VOLUME_UL = 50.0

col1, col2 = st.columns(2)

with col1:
    ts_file = st.file_uploader("1. Upload TapeStation Region CSV File", type=["csv"])
with col2:
    qb_file = st.file_uploader("2. Upload Qubit Run CSV File", type=["csv"])

if ts_file is not None and qb_file is not None:
    try:
        # Load and clean TapeStation Data
        ts_df = pd.read_csv(ts_file)
        
        # Isolate rows marked with '%cfDNA' comment flags (50-700 bp rows)
        ts_cfdna = ts_df[ts_df['Region Comment'] == '%cfDNA'].copy()
        if ts_cfdna.empty:
            # Fallback if comment flags differ: grab every alternate row
            ts_cfdna = ts_df.iloc[::2].copy()
            
        ts_cfdna = ts_cfdna.reset_index(drop=True)

        # Load and clean Qubit Data
        qb_df = pd.read_csv(qb_file)
        qb_df = qb_df.dropna(subset=['Original Sample Conc.']).reset_index(drop=True)

        # Cross-verify row pairing counts
        min_rows = min(len(ts_cfdna), len(qb_df))
        
        if len(ts_cfdna) != len(qb_df):
            st.warning(f"⚠️ Row mismatch! TapeStation regions: {len(ts_cfdna)} | Qubit samples: {len(qb_df)}. Processing first {min_rows} records.")
        
        # Build unified execution dataset
        merged_records = []
        for i in range(min_rows):
            well_id = ts_cfdna.loc[i, 'WellId']
            sample_desc = ts_cfdna.loc[i, 'Sample Description']
            pct_of_total = float(ts_cfdna.loc[i, '% of Total'])
            
            qubit_name = qb_df.loc[i, 'Sample Name']
            raw_qubit_conc = float(qb_df.loc[i, 'Original Sample Conc.'])
            
            # 1. Calculate the raw concentration for this specific region (ng/uL)
            calculated_ng_ul = raw_qubit_conc * (pct_of_total / 100.0)
            
            # 2. Multiply by 50uL total volume to find total mass (ng) in the tube
            total_mass_ng = calculated_ng_ul * TOTAL_VOLUME_UL
            
            merged_records.append({
                "Well ID": well_id,
                "Sample Description": sample_desc,
                "TapeStation % of Total": f"{pct_of_total}%",
                "Qubit Tube Name": qubit_name,
                "Raw Qubit (ng/µL)": raw_qubit_conc,
                "Calculated Region (ng/µL)": round(calculated_ng_ul, 4),
                "Total Regional Mass (ng in 50µL)": round(total_mass_ng, 2)
            })
            
        result_df = pd.DataFrame(merged_records)
        
        # Summary Metrics Panel
        st.subheader("📊 Yield Metrics Summary")
        m1, m2, m3 = st.columns(3)
        m1.metric("Total Samples Processed", len(result_df))
        m2.metric("Avg Region Concentration", f"{result_df['Calculated Region (ng/µL)'].mean():.2f} ng/µL")
        m3.metric("Avg Tube Mass Yield", f"{result_df['Total Regional Mass (ng in 50µL)'].mean():.2f} ng")
        
        # Interactive Grid View
        st.subheader("📋 Final Dataset View")
        st.dataframe(result_df, use_container_width=True)
        
        # Generate export CSV string buffer 
        csv_buffer = io.StringIO()
        result_df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()
        
        st.download_button(
            label="📥 Download Total Mass Calculation CSV",
            data=csv_data,
            file_name="total_tube_mass_yields.csv",
            mime="text/csv"
        )
        st.success("✅ Dataset and mass calculations generated successfully.")
        
    except Exception as e:
        st.error(f"Execution Error occurred processing data layouts: {e}")

