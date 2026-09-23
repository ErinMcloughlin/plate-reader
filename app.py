import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="Direct ID NGS Calculator", page_icon="🧬", layout="wide")

st.title("🧬 Alphanumeric Sample ID Matching Calculator")
st.write("This app joins data logs using strict, direct text matching on your sample identity columns.")

# Volume factor constant definition
TOTAL_VOLUME_UL = 50.0

col1, col2 = st.columns(2)

with col1:
    ts_file = st.file_uploader("1. Upload TapeStation Region CSV File", type=["csv"])
with col2:
    qb_file = st.file_uploader("2. Upload Qubit Run CSV File", type=["csv"])

if ts_file is not None and qb_file is not None:
    try:
        # 1. Read and clean TapeStation Data
        ts_df = pd.read_csv(ts_file, encoding='latin1')
        ts_df.columns = ts_df.columns.str.strip()
        
        # Standardize matching key formats to clean strings
        if 'Sample Description' in ts_df.columns:
            ts_df['Sample Description'] = ts_df['Sample Description'].astype(str).str.strip()
        else:
            st.error("❌ TapeStation file is missing the 'Sample Description' column.")
            st.stop()
            
        # Isolate rows where tracking start window is exactly 50
        ts_df['From [bp]'] = pd.to_numeric(ts_df['From [bp]'], errors='coerce')
        ts_filtered = ts_df[ts_df['From [bp]'] == 50].copy()

        # 2. Read and clean Qubit Data
        qb_df = pd.read_csv(qb_file, encoding='latin1')
        qb_df.columns = qb_df.columns.str.strip()
        
        # Flexibly locate the Qubit identifier column (accepts 'Sample Description', 'Sample Name', or 'Sample ID')
        qb_id_col = None
        for col_name in ['Sample Description', 'Sample Name', 'Sample ID']:
            if col_name in qb_df.columns:
                qb_id_col = col_name
                break
                
        if qb_id_col:
            qb_df[qb_id_col] = qb_df[qb_id_col].astype(str).str.strip()
            # Rename temporarily to create a perfect mirror key for pandas merge functionality
            qb_df = qb_df.rename(columns={qb_id_col: 'Sample Description'})
        else:
            st.error("❌ Qubit file is missing an identifier column like 'Sample Description' or 'Sample Name'.")
            st.stop()
            
        qb_df = qb_df.dropna(subset=['Original Sample Conc.'])

        # 3. Perform a strict inner merge matching the exact string keys
        merged_df = pd.merge(ts_filtered, qb_df, on='Sample Description', how='inner')

        if merged_df.empty:
            st.error("❌ Failed to pair samples. No exact alphanumeric text matches were found between the files.")
            st.info("💡 Ensure both files use identical strings (e.g., both contain 'exDNA26012358').")
            st.stop()

        # 4. Perform the mass balance calculations
        merged_df['% of Total'] = pd.to_numeric(merged_df['% of Total'], errors='coerce')
        merged_df['Original Sample Conc.'] = pd.to_numeric(merged_df['Original Sample Conc.'], errors='coerce')
        
        merged_df['Calculated Region (ng/µL)'] = merged_df['Original Sample Conc.'] * (merged_df['% of Total'] / 100.0)
        merged_df['Total Regional Mass (ng in 50µL)'] = merged_df['Calculated Region (ng/µL)'] * TOTAL_VOLUME_UL

        # Clean display presentation formatting arrays
        display_df = pd.DataFrame({
            "Well ID": merged_df['WellId'],
            "Sample Description": merged_df['Sample Description'],
            "Region Window": "50-" + merged_df['To [bp]'].astype(str) + " bp",
            "TapeStation % of Total": merged_df['% of Total'].astype(str) + "%",
            "Raw Qubit (ng/µL)": merged_df['Original Sample Conc.'],
            "Calculated Region (ng/µL)": merged_df['Calculated Region (ng/µL)'].round(4),
            "Total Regional Mass (ng in 50µL)": merged_df['Total Regional Mass (ng in 50µL)'].round(2)
        })

        # Summary Metrics Panels
        st.subheader("📊 Cross-Matched Run Analysis Summary")
        m1, m2, m3 = st.columns(3)
        m1.metric("Successfully Paired Records", len(display_df))
        m2.metric("Avg Target Pool Concentration", f"{display_df['Calculated Region (ng/µL)'].mean():.2f} ng/µL")
        m3.metric("Avg Calculated Yield Weight", f"{display_df['Total Regional Mass (ng in 50µL)'].mean():.2f} ng")

        # Display Data Spreadsheet Grid Output Panel
        st.subheader("📋 Matched Alphanumeric Output Matrix")
        st.dataframe(display_df, use_container_width=True)

        # Generate download export system logic triggers
        csv_buffer = io.StringIO()
        display_df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()

        st.download_button(
            label="📥 Download Verified Matching Yield Matrix CSV",
            data=csv_data,
            file_name="matched_ngs_yield_report.csv",
            mime="text/csv"
        )
        st.success("✅ Calculations executed safely on matched records data targets.")

    except Exception as e:
        st.error(f"Processing Error: {e}")

