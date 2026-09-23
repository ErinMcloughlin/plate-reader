import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="NGS Library Multi-Run Hub", page_icon="🧬", layout="wide")

st.title("🧬 Alphanumeric Sample ID Matching & Smart Rerun Filter")
st.write("Upload your data logs below. If a rerun is needed, uploading all 4 files will automatically replace initial failures with updated data.")

# Constants
TOTAL_VOLUME_UL = 50.0

# ----------------------------------------------------
# 1. UI UPLOADER LAYOUT
# ----------------------------------------------------
st.subheader("📁 Data Log Inputs")
col1, col2 = st.columns(2)

with col1:
    ts_file_1 = st.file_uploader("TapeStation File (Initial Run)", type=["csv"], key="ts1")
    ts_file_2 = st.file_uploader("TapeStation File (Rerun Only - Optional)", type=["csv"], key="ts2")

with col2:
    qb_file_1 = st.file_uploader("Qubit File (Initial Run)", type=["csv"], key="qb1")
    qb_file_2 = st.file_uploader("Qubit File (Rerun Only - Optional)", type=["csv"], key="qb2")

# Helper function to parse, join, and compute standard parameters
def process_data(ts_in, qb_in):
    if ts_in is None or qb_in is None:
        return None
    
    # Read TapeStation
    ts_df = pd.read_csv(ts_in, encoding='latin1')
    ts_df.columns = ts_df.columns.str.strip()
    ts_df['Sample Description'] = ts_df['Sample Description'].astype(str).str.strip()
    
    # Filter for exact 50bp regions
    ts_df['From [bp]'] = pd.to_numeric(ts_df['From [bp]'], errors='coerce')
    ts_filtered = ts_df[ts_df['From [bp]'] == 50].copy()

    # Read Qubit
    qb_df = pd.read_csv(qb_in, encoding='latin1')
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
        st.error("❌ Qubit file is missing an identifier column.")
        st.stop()
        
    qb_df = qb_df.dropna(subset=['Original Sample Conc.'])

    # Merge data tables
    merged = pd.merge(ts_filtered, qb_df, on='Sample Description', how='inner')
    if merged.empty:
        return pd.DataFrame()

    # Core Math Formulas
    merged['% of Total'] = pd.to_numeric(merged['% of Total'], errors='coerce')
    merged['Original Sample Conc.'] = pd.to_numeric(merged['Original Sample Conc.'], errors='coerce')
    
    merged['Calculated Region (ng/µL)'] = merged['Original Sample Conc.'] * (merged['% of Total'] / 100.0)
    merged['Total Regional Mass (ng in 50µL)'] = merged['Calculated Region (ng/µL)'] * TOTAL_VOLUME_UL

    # Apply Quality Control Flag Checks
    merged['QC Status'] = merged.apply(
        lambda r: "FAIL" if (r['% of Total'] <= 60.0 or r['Total Regional Mass (ng in 50µL)'] <= 10.0) else "PASS", 
        axis=1
    )

    return pd.DataFrame({
        "Well ID": merged['WellId'],
        "Sample Description": merged['Sample Description'],
        "Region Window": "50-" + merged['To [bp]'].astype(str) + " bp",
        "TapeStation % of Total": merged['% of Total'].round(2),
        "Raw Qubit (ng/µL)": merged['Original Sample Conc.'],
        "Calculated Region (ng/µL)": merged['Calculated Region (ng/µL)'].round(4),
        "Total Regional Mass (ng in 50µL)": merged['Total Regional Mass (ng in 50µL)'].round(2),
        "QC Status": merged['QC Status']
    })

# ----------------------------------------------------
# 2. SMART PIPELINE EXECUTION
# ----------------------------------------------------
if ts_file_1 and qb_file_1:
    try:
        # Run Initial Calculations
        final_df = process_data(ts_file_1, qb_file_1)
        
        if final_df.empty:
            st.error("❌ No exact sample ID matches found in the initial data logs.")
            st.stop()

        is_rerun_mode = False

        # If user provides rerun files, execute smart overwrite matching rules
        if ts_file_2 and qb_file_2:
            rerun_df = process_data(ts_file_2, qb_file_2)
            
            if not rerun_df.empty:
                is_rerun_mode = True
                
                # Filter out the initial failures that have an updated entry in the rerun data
                rerun_sample_ids = rerun_df['Sample Description'].tolist()
                
                # Drop old records for these specific samples from our baseline run
                final_df = final_df[~final_df['Sample Description'].isin(rerun_sample_ids)]
                
                # Append the fresh new rerun data records straight onto our dataset matrix
                final_df = pd.concat([final_df, rerun_df], ignore_index=True)
            else:
                st.warning("⚠️ Rerun files were uploaded but no matching sample data points could be paired.")

        # ----------------------------------------------------
        # 3. DASHBOARD SUMMARY DISPLAY PANELS
        # ----------------------------------------------------
        st.write("---")
        if is_rerun_mode:
            st.subheader("📊 Combined Multi-Run Analysis Summary (Reruns Merged)")
        else:
            st.subheader("📊 Initial Run Analysis Summary")

        t_count = len(final_df)
        p_count = len(final_df[final_df['QC Status'] == "PASS"])
        f_count = len(final_df[final_df['QC Status'] == "FAIL"])

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Reported Samples", t_count)
        m2.metric("✅ Passed QC Check", p_count)
        m3.metric("❌ Failed QC Check", f_count, delta=f"-{f_count}" if f_count > 0 else None, delta_color="inverse")
        m4.metric("Avg Tube Mass Yield", f"{final_df['Total Regional Mass (ng in 50µL)'].mean():.2f} ng")

        # ----------------------------------------------------
        # 4. INTERACTIVE VIEW DROPDOWN FILTER
        # ----------------------------------------------------
        st.subheader("📋 Output Matrix Data Viewer")
        status_filter = st.selectbox(
            "Filter table view display parameters:", 
            ["Show All Samples", "Show Only PASS Samples", "Show Only FAIL Samples"]
        )
        
        if status_filter == "Show Only PASS Samples":
            filtered_display = final_df[final_df['QC Status'] == "PASS"]
        elif status_filter == "Show Only FAIL Samples":
            filtered_display = final_df[final_df['QC Status'] == "FAIL"]
        else:
            filtered_display = final_df

        # Apply colorful background highlights to pass/fail status cells
        def color_qc(val):
            return 'background-color: #ffcccc; color: #cc0000; font-weight: bold' if val == 'FAIL' else 'background-color: #ccffcc; color: #006600; font-weight: bold'

        st.dataframe(
            filtered_display.style.map(color_qc, subset=['QC Status']), 
            use_container_width=True
        )

        # Generate download export system csv string layout buffer options
        # Note: The download always exports the full sheet (including ruruns), regardless of the dropdown filter view
        csv_buffer = io.StringIO()
        final_df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()

        st.download_button(
            label="📥 Download Consolidated Report CSV",
            data=csv_data,
            file_name="consolidated_ngs_yield_report.csv",
            mime="text/csv"
        )
        st.success("✅ Calculations executed safely across matching data streams successfully.")

    except Exception as e:
        st.error(f"Processing Error: {e}")
