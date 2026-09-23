import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="NGS Library Multi-Run Hub", page_icon="🧬", layout="wide")

# ====================================================
# 1. USER VALIDATION CONSTRAINTS & THRESHOLDS SETUP
# ====================================================
with st.container(border=True):
    st.markdown("### 📋 Run Parameters & Study Quality Gates")
    col_input, col_study, col_info = st.columns([1.2, 1.2, 2])
    
    with col_input:
        expected_samples_count = st.number_input(
            "Expected Unique Samples:",
            min_value=1,
            value=84,  
            step=1,
            help="The analytical pipeline will gatekeep processing until verified file rows match this value exactly."
        )
        
    with col_study:
        selected_study = st.selectbox(
            "Select Associated Study Framework:",
            ["Two Tube Kit (HALE)", "Four Tube Kits (Procares)"],
            help="Choosing a study sets the explicit upper limits. Samples breaching these points will flag warning alerts."
        )
        
    with col_info:
        if "HALE" in selected_study:
            q_lim, ts_lim = "1.318 ng/µL", "89.69%"
        else:
            q_lim, ts_lim = "3.68 ng/µL", "92.89%"
            
        st.markdown(f"""
        <div style="background-color: #f8f9fa; padding: 12px; border-radius: 5px; border-left: 4px solid #0288d1; margin-top: 5px;">
            <p style="margin: 0; font-size: 13px; font-weight: bold; color: #0288d1;">🛡️ Active Upper Limit Warning Matrices</p>
            <p style="margin: 5px 0 0 0; font-size: 12px; color: #333;"><b>Qubit Concentration Limit:</b> Above {q_lim}</p>
            <p style="margin: 2px 0 0 0; font-size: 12px; color: #333;"><b>TapeStation % of Total Limit:</b> Above {ts_lim}</p>
            <p style="margin: 2px 0 0 0; font-size: 12px; color: #333;"><b>Average Size [bp] Limit:</b> Above 350 bp</p>
        </div>
        """, unsafe_allow_html=True)

st.title("🧬 Plate Upload & Analysis Dashboard")
st.write("Upload your data logs below. Providing rerun files for Round 2 or Round 3 will automatically overwrite older failures with the latest data.")

# Constants
TOTAL_VOLUME_UL = 50.0

# ----------------------------------------------------
# 1. UI UPLOADER LAYOUT (EXPANDED TO 6 FILES)
# ----------------------------------------------------
st.subheader("📁 Data Log Inputs")

tab1, tab2, tab3 = st.tabs(["▶️ Round 1: Initial Baseline", "🔄 Round 2: First Rerun", "🔁 Round 3: Second Rerun"])

with tab1:
    col1, col2 = st.columns(2)
    with col1:
        ts_file_1 = st.file_uploader("TapeStation File (Initial Run)", type=["csv"], key="ts1")
    with col2:
        qb_file_1 = st.file_uploader("Qubit File (Initial Run)", type=["csv"], key="qb1")

with tab2:
    col3, col4 = st.columns(2)
    with col3:
        ts_file_2 = st.file_uploader("TapeStation File (Rerun 1 - Optional)", type=["csv"], key="ts2")
    with col4:
        qb_file_2 = st.file_uploader("Qubit File (Rerun 1 - Optional)", type=["csv"], key="qb2")

with tab3:
    col5, col6 = st.columns(2)
    with col5:
        ts_file_3 = st.file_uploader("TapeStation File (Rerun 2 - Optional)", type=["csv"], key="ts3")
    with col6:
        qb_file_3 = st.file_uploader("Qubit File (Rerun 2 - Optional)", type=["csv"], key="qb3")
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
    ts_calc = ts_df_in.copy()
    ts_calc.columns = ts_calc.columns.str.strip()
    ts_calc['Sample Description'] = ts_calc['Sample Description'].astype(str).str.strip()
    ts_calc['From [bp]'] = pd.to_numeric(ts_calc['From [bp]'], errors='coerce')
    
    all_ts_samples = ts_calc['Sample Description'].dropna().unique()

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
    processed_records = []
    
    for sample_id in all_ts_samples:
        if not sample_id or sample_id in ['nan', '']:
            continue
            
        sample_ts_rows = ts_calc[ts_calc['Sample Description'] == sample_id]
        region_100_row = sample_ts_rows[sample_ts_rows['From [bp]'] == 100]
        qubit_row = qb_calc[qb_calc['Sample Description'] == sample_id]
        
        well_id = sample_ts_rows['WellId'].values if 'WellId' in sample_ts_rows.columns and not sample_ts_rows.empty else "N/A"
        
        if region_100_row.empty:
            raw_qubit = float(qubit_row['Original Sample Conc.'].values) if not qubit_row.empty else 0.0
            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": "No 100bp Region Found",
                "TapeStation % of Total": 0.0,
                "Raw Qubit (ng/µL)": raw_qubit,
                "Calculated Region (ng/µL)": 0.0,
                "Total Regional Mass (ng in 50µL)": 0.0,
                "QC Status": "MISSING 100BP",
                "Average Size [bp]": 0.0
            })
            continue

        if not qubit_row.empty:
            try:
                pct_val = region_100_row['% of Total'].values
                pct_of_total = float(pct_val) if pd.notna(pct_val) and str(pct_val).strip() != "" else 0.0
            except:
                pct_of_total = 0.0
                
            raw_qubit_conc = float(qubit_row['Original Sample Conc.'].values)
            to_bp = region_100_row['To [bp]'].values
            
            calculated_ng_ul = raw_qubit_conc * (pct_of_total / 100.0)
            total_mass_ng = calculated_ng_ul * TOTAL_VOLUME_UL
            
            # --- STUDY UPPER LIMIT CONTROL EVALUATIONS ---
            if "HALE" in selected_study:
                qubit_limit = 1.318
                tapestation_limit = 89.69
            else:
                qubit_limit = 3.68
                tapestation_limit = 92.89

            # Safely extract Average Size scalar value from the array
            try:
                avg_size_val = region_100_row['Average Size [bp]'].values
                avg_size = float(avg_size_val) if pd.notna(avg_size_val) and str(avg_size_val).strip() != "" else 0.0
            except:
                avg_size = 0.0

            # Check if values breach upper limits (including >350 bp limit)
            is_above_qubit = raw_qubit_conc > qubit_limit
            is_above_tapestation = pct_of_total > tapestation_limit
            is_above_size = avg_size > 350.0

            # Determine Active QC Status Hierarchy
            if pct_of_total <= 60.0 or total_mass_ng <= 10.0:
                qc_status = "FAIL"
            elif is_above_qubit or is_above_tapestation or is_above_size:
                qc_status = "ABOVE UPPER LIMIT"
            else:
                qc_status = "PASS"
            
            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": f"100-{to_bp} bp",
                "Average Size [bp]": round(avg_size, 1),  
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
        master_ts_df = pd.read_csv(ts_file_1, encoding='latin1')
        master_qb_df = pd.read_csv(qb_file_1, encoding='latin1')

        # Compute initial run baseline failures to register historic dropouts
        baseline_df = process_data(master_ts_df, master_qb_df)
        
        original_failures = []
        if not baseline_df.empty:
            original_failures = baseline_df[
                baseline_df['QC Status'].isin(['FAIL', 'MISSING 100BP'])
            ]['Sample Description'].tolist()

        ts_clean_cols = master_ts_df.columns.str.strip()
        ts_desc_idx = list(ts_clean_cols).index('Sample Description') if 'Sample Description' in ts_clean_cols else None
        ts_from_idx = list(ts_clean_cols).index('From [bp]') if 'From [bp]' in ts_clean_cols else None
        ts_pct_idx = list(ts_clean_cols).index('% of Total') if '% of Total' in ts_clean_cols else None

        qb_clean_cols = master_qb_df.columns.str.strip()
        qb_id_col_raw = find_qubit_id_col(master_qb_df)
        qb_id_idx = list(master_qb_df.columns).index(qb_id_col_raw) if qb_id_col_raw else None
        qb_conc_idx = list(qb_clean_cols).index('Original Sample Conc.') if 'Original Sample Conc.' in qb_clean_cols else None

        is_rerun_mode = False
        audit_trail_log = []

        # ----------------------------------------------------
        # ROUND 2 PROCESSING (RERUN 1)
        # ----------------------------------------------------
        if ts_file_2 and qb_file_2:
            is_rerun_mode = True
            raw_ts_rerun_1 = pd.read_csv(ts_file_2, encoding='latin1')
            raw_qb_rerun_1 = pd.read_csv(qb_file_2, encoding='latin1')
            
            raw_ts_rerun_1.columns = raw_ts_rerun_1.columns.str.strip()
            raw_qb_rerun_1.columns = raw_qb_rerun_1.columns.str.strip()

            raw_ts_rerun_1['From [bp]'] = pd.to_numeric(raw_ts_rerun_1['From [bp]'], errors='coerce')
            ts_rerun_1_filtered = raw_ts_rerun_1[raw_ts_rerun_1['From [bp]'] == 100]
            
            for _, rerun_row in ts_rerun_1_filtered.iterrows():
                sample_id = str(rerun_row['Sample Description']).strip()
                new_pct = rerun_row['% of Total']
                
                ts_mask = (master_ts_df.iloc[:, ts_desc_idx].astype(str).str.strip() == sample_id) & \
                          (pd.to_numeric(master_ts_df.iloc[:, ts_from_idx], errors='coerce') == 100)
                if ts_mask.any():
                    old_pct = master_ts_df.iloc[ts_mask, ts_pct_idx].values
                    master_ts_df.iloc[ts_mask, ts_pct_idx] = new_pct
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File": "TapeStation",
                        "Rerun Round": "Round 2 (Rerun 1)",
                        "Parameter Updated": "% of Total (100bp Region)",
                        "Old Value": f"{old_pct}%" if pd.notna(old_pct) else "BLANK",
                        "New Value": f"{new_pct}%"
                    })

            for _, rerun_row in raw_qb_rerun_1.dropna(subset=['Original Sample Conc.']).iterrows():
                qb_rerun_id_col = find_qubit_id_col(raw_qb_rerun_1)
                sample_id = str(rerun_row[qb_rerun_id_col]).strip()
                new_conc = rerun_row['Original Sample Conc.']
                
                qb_mask = (master_qb_df.iloc[:, qb_id_idx].astype(str).str.strip() == sample_id)
                if qb_mask.any():
                    old_conc = master_qb_df.iloc[qb_mask, qb_conc_idx].values
                    master_qb_df.iloc[qb_mask, qb_conc_idx] = new_conc
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File": "Qubit",
                        "Rerun Round": "Round 2 (Rerun 1)",
                        "Parameter Updated": "Original Sample Conc. (ng/µL)",
                        "Old Value": f"{old_conc} ng/µL" if pd.notna(old_conc) else "BLANK",
                        "New Value": f"{new_conc} ng/µL"
                    })

        # ----------------------------------------------------
        # ROUND 3 PROCESSING (RERUN 2)
        # ----------------------------------------------------
        if ts_file_3 and qb_file_3:
            is_rerun_mode = True
            raw_ts_rerun_2 = pd.read_csv(ts_file_3, encoding='latin1')
            raw_qb_rerun_2 = pd.read_csv(qb_file_3, encoding='latin1')
            
            raw_ts_rerun_2.columns = raw_ts_rerun_2.columns.str.strip()
            raw_qb_rerun_2.columns = raw_qb_rerun_2.columns.str.strip()

            raw_ts_rerun_2['From [bp]'] = pd.to_numeric(raw_ts_rerun_2['From [bp]'], errors='coerce')
            ts_rerun_2_filtered = raw_ts_rerun_2[raw_ts_rerun_2['From [bp]'] == 100]
            
            for _, rerun_row in ts_rerun_2_filtered.iterrows():
                sample_id = str(rerun_row['Sample Description']).strip()
                new_pct = rerun_row['% of Total']
                
                ts_mask = (master_ts_df.iloc[:, ts_desc_idx].astype(str).str.strip() == sample_id) & \
                          (pd.to_numeric(master_ts_df.iloc[:, ts_from_idx], errors='coerce') == 100)
                if ts_mask.any():
                    old_pct = master_ts_df.iloc[ts_mask, ts_pct_idx].values
                    master_ts_df.iloc[ts_mask, ts_pct_idx] = new_pct
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File": "TapeStation",
                        "Rerun Round": "Round 3 (Rerun 2)",
                        "Parameter Updated": "% of Total (100bp Region)",
                        "Old Value": f"{old_pct}%" if pd.notna(old_pct) else "BLANK",
                        "New Value": f"{new_pct}%"
                    })

            for _, rerun_row in raw_qb_rerun_2.dropna(subset=['Original Sample Conc.']).iterrows():
                qb_rerun_id_col = find_qubit_id_col(raw_qb_rerun_2)
                sample_id = str(rerun_row[qb_rerun_id_col]).strip()
                new_conc = rerun_row['Original Sample Conc.']
                
                qb_mask = (master_qb_df.iloc[:, qb_id_idx].astype(str).str.strip() == sample_id)
                if qb_mask.any():
                    old_conc = master_qb_df.iloc[qb_mask, qb_conc_idx].values
                    master_qb_df.iloc[qb_mask, qb_conc_idx] = new_conc
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File": "Qubit",
                        "Rerun Round": "Round 3 (Rerun 2)",
                        "Parameter Updated": "Original Sample Conc. (ng/µL)",
                        "Old Value": f"{old_conc} ng/µL" if pd.notna(old_conc) else "BLANK",
                        "New Value": f"{new_conc} ng/µL"
                    })

        audit_df = pd.DataFrame(audit_trail_log)
        final_df = process_data(master_ts_df, master_qb_df)

        if final_df.empty:
            st.error("❌ No exact sample ID matches found in the data log parameters.")
            st.stop()

        # Dynamic Recovery Assessment
        if is_rerun_mode and original_failures:
            def adjust_for_recovery(row):
                if row['Sample Description'] in original_failures and row['QC Status'] == 'PASS':
                    return 'RECOVERED'
                return row['QC Status']
            final_df['QC Status'] = final_df.apply(adjust_for_recovery, axis=1)
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
        m5.metric("⚠️ Above Upper Limit", len(final_df[final_df['QC Status'] == "ABOVE UPPER LIMIT"]))
        # ----------------------------------------------------
        # 4. INTERACTIVE VIEW DROPDOWN FILTER
        # ----------------------------------------------------
        st.subheader("📋 Output Matrix Data Viewer")
        status_filter = st.selectbox(
            "Filter table view display parameters:", 
            ["Show All Samples", "Show Only PASS Samples", "Show Only RECOVERED Samples", "Show Only FAIL Samples", "Show Only MISSING 100BP Samples", "Show Only ABOVE UPPER LIMIT Samples"]
        )
        
        if status_filter == "Show Only PASS Samples":
            filtered_display = final_df[final_df['QC Status'] == "PASS"]
        elif status_filter == "Show Only RECOVERED Samples":
            filtered_display = final_df[final_df['QC Status'] == "RECOVERED"]
        elif status_filter == "Show Only FAIL Samples":
            filtered_display = final_df[final_df['QC Status'] == "FAIL"]
        elif status_filter == "Show Only MISSING 100BP Samples":
            filtered_display = final_df[final_df['QC Status'] == "MISSING 100BP"]
        elif status_filter == "Show Only ABOVE UPPER LIMIT Samples":
            filtered_display = final_df[final_df['QC Status'] == "ABOVE UPPER LIMIT"]
        else:
            filtered_display = final_df

        # Apply colorful background highlights to cells dynamically based on study criteria
        def color_qc_row(row):
            styles = [''] * len(row)
            qc_status = row['QC Status']
            
            cols = list(row.index)
            status_idx = cols.index('QC Status') if 'QC Status' in cols else -1
            qubit_idx = cols.index('Raw Qubit (ng/µL)') if 'Raw Qubit (ng/µL)' in cols else -1
            tapestation_idx = cols.index('TapeStation % of Total') if 'TapeStation % of Total' in cols else -1
            size_idx = cols.index('Average Size [bp]') if 'Average Size [bp]' in cols else -1

            if status_idx != -1:
                if qc_status == 'FAIL':
                    styles[status_idx] = 'background-color: #ffcccc; color: #cc0000; font-weight: bold;'
                elif qc_status == 'PASS':
                    styles[status_idx] = 'background-color: #ccffcc; color: #006600; font-weight: bold;'
                elif qc_status == 'RECOVERED':
                    styles[status_idx] = 'background-color: #e6f7ff; color: #0050b3; font-weight: bold;'
                elif qc_status == 'MISSING 100BP':
                    styles[status_idx] = 'background-color: #ffe6cc; color: #cc6600; font-weight: bold;'
                elif qc_status == 'ABOVE UPPER LIMIT':
                    styles[status_idx] = 'background-color: #fff2cc; color: #d68100; font-weight: bold;'

            if "HALE" in selected_study:
                qubit_limit, ts_limit = 1.318, 89.69
            else:
                qubit_limit, ts_limit = 3.68, 92.89

            # Highlight specific cell matrices if breaching thresholds
            if qubit_idx != -1 and float(row['Raw Qubit (ng/µL)']) > qubit_limit:
                styles[qubit_idx] = 'background-color: #fff2cc; color: #d68100; font-weight: bold;'

            if tapestation_idx != -1 and float(row['TapeStation % of Total']) > ts_limit:
                styles[tapestation_idx] = 'background-color: #fff2cc; color: #d68100; font-weight: bold;'

            if size_idx != -1 and float(row['Average Size [bp]']) > 350.0:
                styles[size_idx] = 'background-color: #fff2cc; color: #d68100; font-weight: bold;'

            return styles

        st.dataframe(
            filtered_display.style.apply(color_qc_row, axis=1), 
            use_container_width=True
        )

        # ----------------------------------------------------
        # 5. EXPORT FORMAT GENERATION SYSTEM
        # ----------------------------------------------------
        st.write("---")
        st.subheader("📥 Download Modified Instrument Files & Audit Trail Logs")
        st.info("These files maintain the exact headers and layout rows of your first uploaded files.")

        ts_buffer = io.StringIO()
        master_ts_df.to_csv(ts_buffer, index=False)
        ts_csv_text = ts_buffer.getvalue()
        ts_csv_text = ts_csv_text.replace("Conc. [pg/Âµl]", "Conc. [pg/µl]")
        ts_csv_bytes = ts_csv_text.encode('latin1', errors='ignore')

        qb_buffer = io.StringIO()
        master_qb_df.to_csv(qb_buffer, index=False)
        qb_csv_text = qb_buffer.getvalue()
        qb_csv_bytes = qb_csv_text.encode('latin1', errors='ignore')

        audit_csv_bytes = b""
        if is_rerun_mode and not audit_df.empty:
            audit_buffer = io.StringIO()
            audit_df.to_csv(audit_buffer, index=False)
            audit_csv_bytes = audit_buffer.getvalue().encode('utf-8')

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
            if is_rerun_mode and audit_csv_bytes != b"":
                st.download_button(
                    label="📜 Download Modification Trace Log",
                    data=audit_csv_bytes,
                    file_name="rerun_modification_audit_log.csv",
                    mime="text/csv"
                )
            else:
                st.button("📜 Download Modification Trace Log", disabled=True, help="Upload rerun files to log modifications.")

        st.summary = st.success("✅ Output matrices and validation trace files generated successfully.")

    except Exception as e:
        st.error(f"Processing Error: {e}")

