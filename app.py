import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="NGS QC Yield Calculator", page_icon="🧬", layout="wide")

st.title("🧬 Alphanumeric Sample ID Matching & QC Filter")
st.write("This app joins data logs using strict direct text matching, checks failure metrics, and flags low-quality samples.")

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
        
        qb_id_col = None
        for col_name in ['Sample Description', 'Sample Name', 'Sample ID']:
            if col_name in qb_df.columns:
                qb_id_col = col_name
                break
                
        if qb_id_col:
            qb_df[qb_id_col] = qb_df[qb_id_col].astype(str).str.strip()
            qb_df = qb_df.rename(columns={qb_id_col: 'Sample Description'})
        else:
            st.error("❌ Qubit file is missing an identifier column like 'Sample Description' or 'Sample Name'.")
            st.stop()
            
        qb_df = qb_df.dropna(subset=['Original Sample Conc.'])

        # 3. Perform direct text inner merge matching keys
        merged_df = pd.merge(ts_filtered, qb_df, on='Sample Description', how='inner')

        if merged_df.empty:
            st.error("❌ Failed to pair samples. No exact alphanumeric text matches were found between the files.")
            st.stop()

        # 4. Perform the mass balance calculations
        merged_df['% of Total'] = pd.to_numeric(merged_df['% of Total'], errors='coerce')
        merged_df['Original Sample Conc.'] = pd.to_numeric(merged_df['Original Sample Conc.'], errors='coerce')
        
        merged_df['Calculated Region (ng/µL)'] = merged_df['Original Sample Conc.'] * (merged_df['% of Total'] / 100.0)
        merged_df['Total Regional Mass (ng in 50µL)'] = merged_df['Calculated Region (ng/µL)'] * TOTAL_VOLUME_UL

        # 5. Apply Strict QC Metrics Logic Check
        def check_qc_status(row):
            # Fail if % of Total <= 60 OR Total Regional Mass <= 10
            if row['% of Total'] <= 60.0 or row['Total Regional Mass (ng in 50µL)'] <= 10.0:
                return "FAIL"
            return "PASS"

        merged_df['QC Status'] = merged_df.apply(check_qc_status, axis=1)

        # Clean display columns layout array mapping
        display_df = pd.DataFrame({
            "Well ID": merged_df['WellId'],
            "Sample Description": merged_df['Sample Description'],
            "Region Window": "50-" + merged_df['To [bp]'].astype(str) + " bp",
            "TapeStation % of Total": merged_df['% of Total'].round(2),
            "Raw Qubit (ng/µL)": merged_df['Original Sample Conc.'],
            "Calculated Region (ng/µL)": merged_df['Calculated Region (ng/µL)'].round(4),
            "Total Regional Mass (ng in 50µL)": merged_df['Total Regional Mass (ng in 50µL)'].round(2),
            "QC Status": merged_df['QC Status']
        })

        # Calculate tracking dashboard count variables
        total_count = len(display_df)
        pass_count = len(display_df[display_df['QC Status'] == "PASS"])
        fail_count = len(display_df[display_df['QC Status'] == "FAIL"])

        # Summary Metrics Panels
        st.subheader("📊 Cross-Matched Run Analysis Summary")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Successfully Paired Records", total_count)
        m2.metric("✅ Passed QC Check", pass_count)
        m3.metric("❌ Failed QC Check", fail_count, delta=f"-{fail_count}" if fail_count > 0 else None, delta_color="inverse")
        m4.metric("Avg Tube Mass Yield", f"{display_df['Total Regional Mass (ng in 50µL)'].mean():.2f} ng")

        # Visual filtering control system options
        st.subheader("📋 Output Matrix Data Viewer")
        status_filter = st.selectbox("Filter table view display parameters:", ["Show All Samples", "Show Only PASS Samples", "Show Only FAIL Samples"])
        
        if status_filter == "Show Only PASS Samples":
            filtered_display = display_df[display_df['QC Status'] == "PASS"]
        elif status_filter == "Show Only FAIL Samples":
            filtered_display = display_df[display_df['QC Status'] == "FAIL"]
        else:
            filtered_display = display_df

        # Apply colorful background highlights to pass/fail status cells in the interactive grid UI panel
        def color_qc(val):
            color = 'background-color: #ffcccc; color: #cc0000; font-weight: bold' if val == 'FAIL' else 'background-color: #ccffcc; color: #006600; font-weight: bold'
            return color

        st.dataframe(
            filtered_display.style.map(color_qc, subset=['QC Status']), 
            use_container_width=True
        )

        # Generate download export report files buffers triggers layout elements string mapping configurations
        csv_buffer = io.StringIO()
        display_df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()

        st.download_button(
            label="📥 Download Complete Report Matrix with QC Flags CSV",
            data=csv_data,
            file_name="ngs_yield_qc_report.csv",
            mime="text/csv"
        )
        st.success("✅ Calculations executed safely with verified matching quality controls data metrics.")

    except Exception as e:
        st.error(f"Processing Error: {e}")
