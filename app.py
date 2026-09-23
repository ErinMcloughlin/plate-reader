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
def process_data(ts_df_in, qb_df_in, is_rerun_run=False):
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
        region_100_row = sample_ts_rows[sample_ts_rows['From [bp]'] == 100]
        qubit_row = qb_calc[qb_calc['Sample Description'] == sample_id]
        
        well_id = sample_ts_rows['WellId'].values[0] if 'WellId' in sample_ts_rows.columns and not sample_ts_rows.empty else "N/A"
        
        # Determine baseline failure flags (used to track recovery)
        is_baseline_missing = region_100_row.empty
        baseline_failed = False
        
        if not is_baseline_missing and not qubit_row.empty:
            b_pct = float(region_100_row['% of Total'].values[0])
            b_q_conc = float(qubit_row['Original Sample Conc.'].values[0])
            b_mass = b_q_conc * (b_pct / 100.0) * TOTAL_VOLUME_UL
            if b_pct <= 60.0 or b_mass <= 10.0:
                baseline_failed = True

        # Process active current metrics (including any overwritten files)
        if region_100_row.empty:
            raw_qubit = float(qubit_row['Original Sample Conc.'].values[0]) if not qubit_row.empty else 0.0
            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": "No 100bp Region Found",
                "TapeStation % of Total": 0.0,
                "Raw Qubit (ng/µL)": raw_qubit,
                "Calculated Region (ng/µL)": 0.0,
                "Total Regional Mass (ng in 50µL)": 0.0,
                "QC Status": "MISSING 100BP"
            })
            continue

        if not qubit_row.empty:
            pct_of_total = float(region_100_row['% of Total'].values[0])
            raw_qubit_conc = float(qubit_row['Original Sample Conc.'].values[0])
            to_bp = region_100_row['To [bp]'].values[0]
            
            calculated_ng_ul = raw_qubit_conc * (pct_of_total / 100.0)
            total_mass_ng = calculated_ng_ul * TOTAL_VOLUME_UL
            
            # Determine active evaluation status
            if pct_of_total <= 60.0 or total_mass_ng <= 10.0:
                qc_status = "FAIL"
            else:
                # ✅ CORRECTED LINE: Uses function parameter flag directly
                if (is_baseline_missing or baseline_failed) and is_rerun_run:
                    qc_status = "RECOVERED"
                else:
                    qc_status = "PASS"
            
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
# 2. DATA MERGING & VERIFICATION PIPELINE
# ----------------------------------------------------
if ts_file_1 and qb_file_1:
    try:
        # 1. Load original raw files into memory
        master_ts_df = pd.read_csv(ts_file_1, encoding='latin1')
        master_qb_df = pd.read_csv(qb_file_1, encoding='latin1')

        # 2. Calculate baseline values first to find original dropouts
        baseline_df = process_data(master_ts_df, master_qb_df, is_rerun_run=False)
        
        # Save a list of sample IDs that originally failed or were missing
        original_failures = []
        if not baseline_df.empty:
            original_failures = baseline_df[
                baseline_df['QC Status'].isin(['FAIL', 'MISSING 100BP'])
            ]['Sample Description'].tolist()

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
        is_rerun_mode = False
        # Initialize an empty list to capture our audit trail logs
        audit_trail_log = []

        # 3. Apply Rerun Merging if files are present
        if ts_file_2 and qb_file_2:
            is_rerun_mode = True
            raw_ts_rerun = pd.read_csv(ts_file_2, encoding='latin1')
            raw_qb_rerun = pd.read_csv(qb_file_2, encoding='latin1')
            
            raw_ts_rerun.columns = raw_ts_rerun.columns.str.strip()
            raw_qb_rerun.columns = raw_qb_rerun.columns.str.strip()

            raw_ts_rerun['From [bp]'] = pd.to_numeric(raw_ts_rerun['From [bp]'], errors='coerce')
            ts_rerun_filtered = raw_ts_rerun[raw_ts_rerun['From [bp]'] == 100]
            
            # Trace TapeStation updates
            for _, rerun_row in ts_rerun_filtered.iterrows():
                sample_id = str(rerun_row['Sample Description']).strip()
                new_pct = rerun_row['% of Total']
                
                ts_mask = (master_ts_df.iloc[:, ts_desc_idx].astype(str).str.strip() == sample_id) & \
                          (pd.to_numeric(master_ts_df.iloc[:, ts_from_idx], errors='coerce') == 100)
                
                if ts_mask.any():
                    # Record old value before overwrite
                    old_pct = master_ts_df.iloc[ts_mask, ts_pct_idx].values[0]
                    master_ts_df.iloc[ts_mask, ts_pct_idx] = new_pct
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File Type": "TapeStation",
                        "Data Parameter Updated": "% of Total (100bp Region)",
                        "Original Baseline Value": old_pct,
                        "New Overwritten Value": new_pct
                    })

            # Trace Qubit updates
            for _, rerun_row in raw_qb_rerun.dropna(subset=['Original Sample Conc.']).iterrows():
                qb_rerun_id_col = find_qubit_id_col(raw_qb_rerun)
                sample_id = str(rerun_row[qb_rerun_id_col]).strip()
                new_conc = rerun_row['Original Sample Conc.']
                
                qb_mask = (master_qb_df.iloc[:, qb_id_idx].astype(str).str.strip() == sample_id)
                
                if qb_mask.any():
                    # Record old value before overwrite
                    old_conc = master_qb_df.iloc[qb_mask, qb_conc_idx].values[0]
                    master_qb_df.iloc[qb_mask, qb_conc_idx] = new_conc
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File Type": "Qubit",
                        "Data Parameter Updated": "Original Sample Conc. (ng/µL)",
                        "Original Baseline Value": old_conc,
                        "New Overwritten Value": new_conc
                    })

        # Convert the audit list into a structural tracking DataFrame
        audit_df = pd.DataFrame(audit_trail_log)

        # ----------------------------------------------------
        # CRITICAL VALIDATION CHECK
        # ----------------------------------------------------
        total_reported_samples = len(final_df)
        if total_reported_samples != expected_samples_count:
            st.error(f"🚨 **Sample Count Mismatch! Processing Blocked.**")
            st.error(f"Expected: **{expected_samples_count}** unique samples | Detected in files: **{total_reported_samples}** unique samples.")
            st.info("💡 Please verify your input log files or update the expected sample input number above.")
            st.stop()

        # ----------------------------------------------------
        # 3. DASHBOARD SUMMARY PANEL
        # ----------------------------------------------------
        st.write("---")
        st.subheader("📊 Combined Run Analysis Summary" if is_rerun_mode else "📊 Initial Run Analysis Summary")

        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Total Reported Samples", len(final_df))
        m2.metric("✅ Passed QC Check", len(final_df[final_df['QC Status'] == "PASS"]))
        m3.metric("🚀 Recovered Status", len(final_df[final_df['QC Status'] == "RECOVERED"]))
        m4.metric("❌ Failed QC Check", len(final_df[final_df['QC Status'] == "FAIL"]))
        m5.metric("⚠️ Missing 100bp Regions", len(final_df[final_df['QC Status'] == "MISSING 100BP"]))

        # ----------------------------------------------------
        # 4. INTERACTIVE VIEW DROPDOWN FILTER
        # ----------------------------------------------------
        st.subheader("📋 Output Matrix Data Viewer")
        status_filter = st.selectbox(
            "Filter table view display parameters:", 
            ["Show All Samples", "Show Only PASS Samples", "Show Only RECOVERED Samples", "Show Only FAIL Samples", "Show Only MISSING 100BP Samples"]
        )
        
        if status_filter == "Show Only PASS Samples":
            filtered_display = final_df[final_df['QC Status'] == "PASS"]
        elif status_filter == "Show Only RECOVERED Samples":
            filtered_display = final_df[final_df['QC Status'] == "RECOVERED"]
        elif status_filter == "Show Only FAIL Samples":
            filtered_display = final_df[final_df['QC Status'] == "FAIL"]
        elif status_filter == "Show Only MISSING 100BP Samples":
            filtered_display = final_df[final_df['QC Status'] == "MISSING 100BP"]
        else:
            filtered_display = final_df

        def color_qc(val):
            if val == 'FAIL':
                return 'background-color: #ffcccc; color: #cc0000; font-weight: bold'
            elif val == 'PASS':
                return 'background-color: #ccffcc; color: #006600; font-weight: bold'
            elif val == 'RECOVERED':
                return 'background-color: #e6f7ff; color: #0050b3; font-weight: bold' 
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
        st.subheader("📥 Download Modified Instrument Files & Audit Trail Logs")
        st.info("These files maintain the exact headers and layout rows of your first uploaded files. The change log traces your history entries.")

        # Re-pack instrument files
        ts_buffer = io.StringIO()
        master_ts_df.to_csv(ts_buffer, index=False)
        ts_csv_bytes = ts_buffer.getvalue()

        qb_buffer = io.StringIO()
        master_qb_df.to_csv(qb_buffer, index=False)
        qb_csv_bytes = qb_buffer.getvalue()

        # Re-pack Audit Change Log
        audit_csv_bytes = ""
        if is_rerun_mode and not audit_df.empty:
            audit_buffer = io.StringIO()
            audit_df.to_csv(audit_buffer, index=False)
            audit_csv_bytes = audit_buffer.getvalue()

        # Render 3 columns side-by-side for neat alignment
        dl_col1, dl_col2, dl_col3 = st.columns(3)
        
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
            
        with dl_col3:
            if is_rerun_mode and audit_csv_bytes != "":
                st.download_button(
                    label="📜 Download Rerun Modification Trace Log",
                    data=audit_csv_bytes,
                    file_name="rerun_modification_audit_log.csv",
                    mime="text/csv"
                )
            else:
                st.button("📜 Download Rerun Modification Trace Log", disabled=True, help="This option unlocks only when 4 files are parsed.")
                
        st.success("✅ Output matrices and validation trace files generated successfully.")


