import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="NGS Library Multi-Run Hub", page_icon="🧬", layout="wide")
# ----------------------------------------------------
# 1. USER INPUTS & FILE UPLOADER LAYOUT
# ----------------------------------------------------
# Wrap the constraint in a clean card container block
with st.container(border=True):
    st.markdown("### 📋 Sample Number Check")
    
    # Use columns to keep the width concise and aligned
    input_col, info_col = st.columns([1, 2])
    
    with input_col:
        expected_samples_count = st.number_input(
            "Expected Unique Samples:",
            min_value=1,
            value=1,  # Setting default to 10 prevents immediate 0 errors on load
            step=1,
            help="The analytical pipeline will gatekeep processing until your verified file rows match this value exactly."
        )
        
    with info_col:
        st.markdown("<br>", unsafe_allow_html=True) # Simple vertical spacer alignment
        st.caption(
            "💡 **Quality Gate:** This value checks unique entries in your `Sample Description` log metrics. "
        )
st.title("🧬 Plate Upload")
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

# Helper function to find the flexible Qubit ID column name
def find_qubit_id_col(df):
    for col_name in ['Sample Description', 'Sample Name', 'Sample ID']:
        if col_name in df.columns:
            return col_name
    stripped_cols = {c.strip(): c for c in df.columns}
    for col_name in ['Sample Description', 'Sample Name', 'Sample ID']:
        if col_name in stripped_cols:
            return stripped_cols[col_name]
    return None

# Helper function to process data strictly for dashboard view computations
def process_data(ts_df_in, qb_df_in):
    # Standardize column headers for reliable merge math without breaking raw copies
    ts_calc = ts_df_in.copy()
    ts_calc.columns = ts_calc.columns.str.strip()
    ts_calc['Sample Description'] = ts_calc['Sample Description'].astype(str).str.strip()
    ts_calc['From [bp]'] = pd.to_numeric(ts_calc['From [bp]'], errors='coerce')
    
    # Identify ALL unique samples present in the TapeStation file to check for completeness
    all_ts_samples = ts_calc['Sample Description'].dropna().unique()

    # Read Qubit
    qb_calc = qb_df_in.copy()
    qb_calc.columns = qb_calc.columns.str.strip()
    
    qb_id_col = find_qubit_id_col(qb_calc)
    if qb_id_col:
        qb_calc[qb_id_col] = qb_calc[qb_id_col].astype(str).str.strip()
        qb_calc = qb_calc.rename(columns={qb_id_col: 'Sample Description'})
    else:
        st.error("❌ Qubit file is missing an identifier column.")
        st.stop()
        
    qb_calc = qb_calc.dropna(subset=['Original Sample Conc.'])

    # Build the tracking array
    processed_records = []
    
    for sample_id in all_ts_samples:
        if not sample_id or sample_id == 'nan' or sample_id == '':
            continue
            
        # Isolate rows for this specific sample
        sample_ts_rows = ts_calc[ts_calc['Sample Description'] == sample_id]
        # Look for the exact 100bp region row
        region_100_row = sample_ts_rows[sample_ts_rows['From [bp]'] == 100]
        # Look for matching Qubit entry
        qubit_row = qb_calc[qb_calc['Sample Description'] == sample_id]
        
        # Grab well ID safely from values array
        well_id = sample_ts_rows['WellId'].values[0] if 'WellId' in sample_ts_rows.columns and not sample_ts_rows.empty else "N/A"
        
        # FLAG CONDITIONAL: If the sample exists but lacks a 100bp row entry
        if region_100_row.empty:
            raw_qubit = float(qubit_row['Original Sample Conc.'].values[0]) if not qubit_row.empty else 0.0
            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": "No 100bp Region Found",
                "TapeStation % of Total": 0.0,
                "Raw Qubit (ng/µL)": raw_qubit,
                "Calculated Region (ng/µL)": 0.0,
                "Total Regional Mass (ng in 100µL)": 0.0,
                "QC Status": "MISSING 100BP"
            })
            continue

        # If it has the 100bp row, check if it also matches a Qubit record
        if not qubit_row.empty:
            pct_of_total = float(region_100_row['% of Total'].values[0])
            raw_qubit_conc = float(qubit_row['Original Sample Conc.'].values[0])
            to_bp = region_100_row['To [bp]'].values[0]
            
            # Core Math Formulas
            calculated_ng_ul = raw_qubit_conc * (pct_of_total / 100.0)
            total_mass_ng = calculated_ng_ul * TOTAL_VOLUME_UL
            
            qc_status = "FAIL" if (pct_of_total <= 60.0 or total_mass_ng <= 10.0) else "PASS"
            
            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": f"100-{to_bp} bp",
                "TapeStation % of Total": round(pct_of_total, 2),
                "Raw Qubit (ng/µL)": raw_qubit_conc,
                "Calculated Region (ng/µL)": round(calculated_ng_ul, 4),
                "Total Regional Mass (ng in 50µL)": round(total_mass_ng, 2),
                "QC Status": qc_status
            })

    return pd.DataFrame(processed_records)

# ----------------------------------------------------
# 2. DATA MERGING PIPELINE
# ----------------------------------------------------
if ts_file_1 and qb_file_1:
    try:
        # Load raw files into memory preserving formatting exactly (Latin-1 preserves micro symbol 'µ')
        master_ts_df = pd.read_csv(ts_file_1, encoding='latin1')
        master_qb_df = pd.read_csv(qb_file_1, encoding='latin1')

        # Identify raw column indices for strict structural manipulation
        ts_clean_cols = master_ts_df.columns.str.strip()
        ts_desc_idx = list(ts_clean_cols).index('Sample Description') if 'Sample Description' in ts_clean_cols else None
        ts_from_idx = list(ts_clean_cols).index('From [bp]') if 'From [bp]' in ts_clean_cols else None
        ts_pct_idx = list(ts_clean_cols).index('% of Total') if '% of Total' in ts_clean_cols else None

        qb_clean_cols = master_qb_df.columns.str.strip()
        qb_id_col_raw = find_qubit_id_col(master_qb_df)
        qb_id_idx = list(master_qb_df.columns).index(qb_id_col_raw) if qb_id_col_raw else None
        qb_conc_idx = list(qb_clean_cols).index('Original Sample Conc.') if 'Original Sample Conc.' in qb_clean_cols else None

        is_rerun_mode = False

        # In-Place replacement block if rerun pairs are uploaded
        if ts_file_2 and qb_file_2:
            is_rerun_mode = True
            raw_ts_rerun = pd.read_csv(ts_file_2, encoding='latin1')
            raw_qb_rerun = pd.read_csv(qb_file_2, encoding='latin1')
            
            raw_ts_rerun.columns = raw_ts_rerun.columns.str.strip()
            raw_qb_rerun.columns = raw_qb_rerun.columns.str.strip()

            # Overwrite TapeStation Master Sheet records where From [bp] == 100
            raw_ts_rerun['From [bp]'] = pd.to_numeric(raw_ts_rerun['From [bp]'], errors='coerce')
            ts_rerun_filtered = raw_ts_rerun[raw_ts_rerun['From [bp]'] == 100]
            
            for _, rerun_row in ts_rerun_filtered.iterrows():
                sample_id = str(rerun_row['Sample Description']).strip()
                new_pct = rerun_row['% of Total']
                
                ts_mask = (master_ts_df.iloc[:, ts_desc_idx].astype(str).str.strip() == sample_id) & \
                          (pd.to_numeric(master_ts_df.iloc[:, ts_from_idx], errors='coerce') == 100)
                if ts_mask.any():
                    master_ts_df.iloc[ts_mask, ts_pct_idx] = new_pct

            # Overwrite Qubit Master Sheet records
            for _, rerun_row in raw_qb_rerun.dropna(subset=['Original Sample Conc.']).iterrows():
                qb_rerun_id_col = find_qubit_id_col(raw_qb_rerun)
                sample_id = str(rerun_row[qb_rerun_id_col]).strip()
                new_conc = rerun_row['Original Sample Conc.']
                
                qb_mask = (master_qb_df.iloc[:, qb_id_idx].astype(str).str.strip() == sample_id)
                if qb_mask.any():
                    master_qb_df.iloc[qb_mask, qb_conc_idx] = new_conc

        # Compute data strictly for screen visualization grid parameters
        final_df = process_data(master_ts_df, master_qb_df)

        if final_df.empty:
            st.error("❌ No exact sample ID matches found in the data log parameters.")
            st.stop()
        # Compute data strictly for screen visualization grid parameters
        final_df = process_data(master_ts_df, master_qb_df)

        if final_df.empty:
            st.error("❌ No exact sample ID matches found in the data log parameters.")
            st.stop()

        # ----------------------------------------------------
        # CRITICAL VALIDATION CHECK (INSERT THIS BLOCK HERE)
        # ----------------------------------------------------
        total_reported_samples = len(final_df)
        
        if total_reported_samples != expected_samples_count:
            st.error(f"🚨 **Sample Count Mismatch! Processing Blocked.**")
            st.error(f"Expected: **{expected_samples_count}** unique samples | Detected in files: **{total_reported_samples}** unique samples.")
            st.info("💡 Please verify your input log files or update the expected sample input number above to match your run sequence layout.")
            st.stop() # Stops execution instantly, blocking downstream dashboard elements

        # ----------------------------------------------------
        # 3. DASHBOARD SUMMARY PANEL
        # ----------------------------------------------------
        st.write("---")
        if is_rerun_mode:
            st.subheader("📊 Combined Multi-Run Analysis Summary (Reruns Merged)")
        else:
            st.subheader("📊 Initial Run Analysis Summary")

        t_count = len(final_df)
        p_count = len(final_df[final_df['QC Status'] == "PASS"])
        f_count = len(final_df[final_df['QC Status'] == "FAIL"])
        m_count = len(final_df[final_df['QC Status'] == "MISSING 100BP"])

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Total Reported Samples", t_count)
        m2.metric("✅ Passed QC Check", p_count)
        m3.metric("❌ Failed QC Check", f_count, delta=f"-{f_count}" if f_count > 0 else None, delta_color="inverse")
        m4.metric("⚠️ Missing 100bp Regions", m_count)

        # ----------------------------------------------------
        # 4. INTERACTIVE VIEW DROPDOWN FILTER
        # ----------------------------------------------------
        st.subheader("📋 Output Matrix Data Viewer")
        status_filter = st.selectbox(
            "Filter table view display parameters:", 
            ["Show All Samples", "Show Only PASS Samples", "Show Only FAIL Samples", "Show Only MISSING 100BP Samples"]
        )
        
        if status_filter == "Show Only PASS Samples":
            filtered_display = final_df[final_df['QC Status'] == "PASS"]
        elif status_filter == "Show Only FAIL Samples":
            filtered_display = final_df[final_df['QC Status'] == "FAIL"]
        elif status_filter == "Show Only MISSING 100BP Samples":
            filtered_display = final_df[final_df['QC Status'] == "MISSING 100BP"]
        else:
            filtered_display = final_df

        # Apply colorful background highlights to cells dynamically
        def color_qc(val):
            if val == 'FAIL':
                return 'background-color: #ffcccc; color: #cc0000; font-weight: bold'
            elif val == 'PASS':
                return 'background-color: #ccffcc; color: #006600; font-weight: bold'
            elif val == 'MISSING 100BP':
                return 'background-color: #ffe6cc; color: #cc6600; font-weight: bold'
            return ''

        st.dataframe(
            filtered_display.style.map(color_qc, subset=['QC Status']),
            use_container_width=True
        )

        # ----------------------------------------------------
        # 5. EXPORT FORMAT GENERATION SYSTEM
        # ----------------------------------------------------
        st.write("---")
        st.subheader("📥 Download Modified Instrument Files")
        st.info("These files maintain the exact headers and layout rows of your first uploaded files. If reruns were provided, the data is integrated directly into them.")

        ts_buffer = io.StringIO()
        master_ts_df.to_csv(ts_buffer, index=False)
        ts_csv_bytes = ts_buffer.getvalue()

        qb_buffer = io.StringIO()
        master_qb_df.to_csv(qb_buffer, index=False)
        qb_csv_bytes = qb_buffer.getvalue()

        dl_col1, dl_col2 = st.columns(2)
        
        with dl_col1:
            st.download_button(
                label="📥 Download Updated TapeStation File",
                data=ts_csv_bytes,
                file_name="updated_tapestation_report.csv",
                mime="text/csv"
            )
            
        with dl_col2:
            st.download_button(
                label="📥 Download Updated Qubit File",
                data=qb_csv_bytes,
                file_name="updated_qubit_report.csv",
                mime="text/csv"
            )
            
        st.success("✅ Instrument files exported successfully.")

    except Exception as e:
        st.error(f"Processing Error: {e}")
