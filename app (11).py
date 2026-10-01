import streamlit as st
import pandas as pd
import io
import numpy as np

st.set_page_config(page_title="TapeStation_Qubit_Analysis", page_icon="🧬", layout="wide")

# ====================================================
# 0. GLOBAL QUBIT OUTLIER / RUN REVIEW CONSTANTS
# ====================================================
# 3SD upper limit for Qubit concentration, from the raw (untransformed) historical data:
#   n = 2,164 samples | mean = 1.659 ng/µL | SD = 3.509 ng/µL | mean + 3SD = 12.18 ng/µL
QUBIT_3SD_LIMIT = 12.18          # Samples with Raw Qubit above this are flagged for rerun
# Qubit level used for the run-level review check depends on the study selected:
#   Two Tube Kit (HALE) -> 1.318 ng/µL | Four Tube Kits (Procares) -> 3.68 ng/µL
RUN_REVIEW_QUBIT_LIMITS = {"HALE": 1.318, "Procares": 3.68}
RUN_REVIEW_FRACTION = 0.10       # Run flagged if MORE than 10% of samples exceed RUN_REVIEW_QUBIT_LIMIT

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
            value=1,  
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
            RUN_REVIEW_QUBIT_LIMIT = RUN_REVIEW_QUBIT_LIMITS["HALE"]
        else:
            q_lim, ts_lim = "3.68 ng/µL", "92.89%"
            RUN_REVIEW_QUBIT_LIMIT = RUN_REVIEW_QUBIT_LIMITS["Procares"]
            
        st.markdown(f"""
        <div style="background-color: #f8f9fa; padding: 12px; border-radius: 5px; border-left: 4px solid #0288d1; margin-top: 5px;">
            <p style="margin: 0; font-size: 13px; font-weight: bold; color: #0288d1;">🛡️ Active Upper Limit Warning Matrices</p>
            <p style="margin: 5px 0 0 0; font-size: 12px; color: #333;"><b>Qubit Concentration Limit:</b> Above {q_lim}</p>
            <p style="margin: 2px 0 0 0; font-size: 12px; color: #333;"><b>TapeStation % of Total Limit:</b> Above {ts_lim}</p>
            <p style="margin: 2px 0 0 0; font-size: 12px; color: #333;"><b>Average Size [bp] Limit:</b> Above 350 bp</p>
            <p style="margin: 2px 0 0 0; font-size: 12px; color: #cc0000;"><b>Qubit 3SD Outlier (Rerun):</b> Above {QUBIT_3SD_LIMIT} ng/µL</p>
            <p style="margin: 2px 0 0 0; font-size: 12px; color: #cc0000;"><b>Run Review:</b> More than {RUN_REVIEW_FRACTION:.0%} of samples above {RUN_REVIEW_QUBIT_LIMIT} ng/µL</p>
        </div>
        """, unsafe_allow_html=True)

st.title("🧬 Plate Upload & Analysis Dashboard")
st.write("Upload your data logs below. Providing rerun files for Round 2 or Round 3 will automatically trigger %CV evaluations and outlier analysis.")

# Constants
TOTAL_VOLUME_UL = 50.0
MIN_TOTAL_MASS_NG = 10.0    # Mass pass: Total Regional Mass (Qubit x TapeStation % of Total x 50 µL) must be >= 10 ng
REGION_FROM_BP = 100        # TapeStation region row to use: the one whose "From [bp]" equals this value

# ----------------------------------------------------
# 1. UI UPLOADER LAYOUT (EXPANDED TO 6 FILES)
# ----------------------------------------------------
st.subheader("📁 Data Log Inputs")

tab1, tab2, tab3 = st.tabs(["▶️ Round 1: Initial Baseline", "🔄 Round 2: First Rerun", "🔁 Round 3: Second Rerun"])

with tab1:
    col1, col2 = st.columns(2)
    with col1:
        ts_file_1 = st.file_uploader("TapeStation File (Initial Run)", type=["csv", "xlsx", "xls"], key="ts1")
    with col2:
        qb_file_1 = st.file_uploader("Qubit File (Initial Run)", type=["csv", "xlsx", "xls"], key="qb1")

with tab2:
    col3, col4 = st.columns(2)
    with col3:
        ts_file_2 = st.file_uploader("TapeStation File (Rerun 1 - Optional)", type=["csv", "xlsx", "xls"], key="ts2")
    with col4:
        qb_file_2 = st.file_uploader("Qubit File (Rerun 1 - Optional)", type=["csv", "xlsx", "xls"], key="qb2")

with tab3:
    col5, col6 = st.columns(2)
    with col5:
        ts_file_3 = st.file_uploader("TapeStation File (Rerun 2 - Optional)", type=["csv", "xlsx", "xls"], key="ts3")
    with col6:
        qb_file_3 = st.file_uploader("Qubit File (Rerun 2 - Optional)", type=["csv", "xlsx", "xls"], key="qb3")

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

# Helper: read an uploaded CSV or Excel file
def read_data_file(uploaded_file, encoding='latin1'):
    name = getattr(uploaded_file, 'name', '') or ''
    if name.lower().endswith(('.xlsx', '.xls')):
        uploaded_file.seek(0)
        return pd.read_excel(uploaded_file)
    return pd.read_csv(uploaded_file, encoding=encoding)

# Helper: %CV from a list of measurements (ignores blanks/NaN). Returns None if not computable.
def compute_cv(values):
    clean = [float(v) for v in values if v is not None and pd.notna(v)]
    if len(clean) < 2:
        return None
    mean_val = np.mean(clean)
    if mean_val <= 0:
        return None
    return (np.std(clean, ddof=1) / mean_val) * 100

# Helper: format %CV text the same way for every sample ("TS: x% | QB: y%")
def format_cv(ts_cv, qb_cv, suffix=""):
    ts_txt = f"{ts_cv:.1f}%" if ts_cv is not None else "N/A"
    qb_txt = f"{qb_cv:.1f}%" if qb_cv is not None else "N/A"
    return f"TS: {ts_txt} | QB: {qb_txt}{suffix}"

# Helper: marks, for each sample, the TapeStation region row whose "From [bp]" equals REGION_FROM_BP (100 bp)
def first_region_mask(df, desc_idx, from_idx):
    desc = df.iloc[:, desc_idx].astype(str).str.strip()
    frm = pd.to_numeric(df.iloc[:, from_idx], errors='coerce')
    mask = pd.Series(False, index=df.index)
    first_idx = desc[frm == REGION_FROM_BP].drop_duplicates(keep='first').index
    mask.loc[first_idx] = True
    return mask

# %CV acceptance threshold used for NO REPEAT NEEDED decisions
CV_LIMIT = 20.0

# Helper function to process data strictly for dashboard view computations
def process_data(ts_df_in, qb_df_in, status_overrides=None, calculated_cv_map=None):
    if status_overrides is None:
        status_overrides = {}
    if calculated_cv_map is None:
        calculated_cv_map = {}

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
        # Use the sample's region row that starts at REGION_FROM_BP (100 bp)
        region_100_row = sample_ts_rows[sample_ts_rows['From [bp]'] == REGION_FROM_BP].head(1)
        qubit_row = qb_calc[qb_calc['Sample Description'] == sample_id]
        
        # Protected index lookup for scalar conversions
        well_id = sample_ts_rows['WellId'].iloc[0] if 'WellId' in sample_ts_rows.columns and not sample_ts_rows.empty else "N/A"
        cv_display_val = calculated_cv_map.get(sample_id, "")

        if region_100_row.empty:
            raw_qubit = float(qubit_row['Original Sample Conc.'].iloc[0]) if not qubit_row.empty else 0.0
            current_qc = status_overrides.get(sample_id, "FAIL")

            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": "No Region Found",
                "TapeStation % of Total": 0.0,
                "Raw Qubit (ng/µL)": raw_qubit,
                "Calculated Region (ng/µL)": 0.0,
                "Qubit Total (ng in 50µL)": round(raw_qubit * TOTAL_VOLUME_UL, 2),
                "Total Regional Mass (ng in 50µL)": 0.0,
                "QC Status": current_qc,
                "Average Size [bp]": 0.0,
                "Calculated %CV (Reruns)": cv_display_val
            })
            continue

        if not qubit_row.empty:
            try:
                pct_val = region_100_row['% of Total'].iloc[0]
                pct_of_total = float(pct_val) if pd.notna(pct_val) and str(pct_val).strip() != "" else 0.0
            except:
                pct_of_total = 0.0
                
            raw_qubit_conc = float(qubit_row['Original Sample Conc.'].iloc[0])
            to_bp = region_100_row['To [bp]'].iloc[0] if 'To [bp]' in region_100_row.columns else ""
            from_bp = region_100_row['From [bp]'].iloc[0]
            from_bp = int(from_bp) if float(from_bp).is_integer() else from_bp
            try:
                to_bp = int(float(to_bp)) if float(to_bp).is_integer() else to_bp
            except (TypeError, ValueError):
                pass
            
            calculated_ng_ul = raw_qubit_conc * (pct_of_total / 100.0)
            total_mass_ng = calculated_ng_ul * TOTAL_VOLUME_UL   # used for pass/fail (must be >= MIN_TOTAL_MASS_NG)
            qubit_total_ng = raw_qubit_conc * TOTAL_VOLUME_UL
            
            # --- STUDY UPPER LIMIT CONTROL EVALUATIONS ---
            if "HALE" in selected_study:
                qubit_limit = 1.318
                tapestation_limit = 89.69
            else:
                qubit_limit = 3.68
                tapestation_limit = 92.89

            try:
                avg_size_val = region_100_row['Average Size [bp]'].iloc[0]
                avg_size = float(avg_size_val) if pd.notna(avg_size_val) and str(avg_size_val).strip() != "" else 0.0
            except:
                avg_size = 0.0

            is_above_qubit = raw_qubit_conc > qubit_limit
            is_above_tapestation = pct_of_total > tapestation_limit
            is_above_size = avg_size > 350.0
            is_above_3sd = raw_qubit_conc > QUBIT_3SD_LIMIT

            if sample_id in status_overrides:
                qc_status = status_overrides[sample_id]
            elif pct_of_total <= 60.0 or total_mass_ng < MIN_TOTAL_MASS_NG:
                qc_status = "FAIL"
            elif is_above_3sd:
                # Qubit reading is a 3SD outlier vs. historical data -> flag for rerun
                qc_status = "FAIL"
            elif is_above_qubit or is_above_tapestation or is_above_size:
                qc_status = "ABOVE UPPER LIMIT"
            else:
                qc_status = "PASS"
            
            processed_records.append({
                "Well ID": well_id,
                "Sample Description": sample_id,
                "Region Window": f"{from_bp}-{to_bp} bp",
                "Average Size [bp]": round(avg_size, 1),  
                "TapeStation % of Total": round(pct_of_total, 2),
                "Raw Qubit (ng/µL)": raw_qubit_conc,
                "Calculated Region (ng/µL)": round(calculated_ng_ul, 4),
                "Qubit Total (ng in 50µL)": round(qubit_total_ng, 2),
                "Total Regional Mass (ng in 50µL)": round(total_mass_ng, 2),
                "QC Status": qc_status,
                "Calculated %CV (Reruns)": cv_display_val
            })

    return pd.DataFrame(processed_records)

# ----------------------------------------------------
# 2. DATA MERGING & VERIFICATION PIPELINE
# ----------------------------------------------------
if ts_file_1 and qb_file_1:
    try:
        master_ts_df = read_data_file(ts_file_1, encoding='latin1')
        master_qb_df = read_data_file(qb_file_1, encoding='latin1')

        pipeline_status_overrides = {}
        pipeline_cv_reporting = {}

        baseline_df = process_data(master_ts_df, master_qb_df)
        
        original_failures = []
        if not baseline_df.empty:
            original_failures = baseline_df[
                baseline_df['QC Status'] == 'FAIL'
            ]['Sample Description'].tolist()

        # --- RUN-LEVEL REVIEW CHECK (evaluated on the initial run input) ---
        run_review_required = False
        run_review_count = 0
        run_review_total = len(baseline_df)
        run_review_pct = 0.0
        if run_review_total > 0:
            run_review_count = int((baseline_df['Raw Qubit (ng/µL)'] > RUN_REVIEW_QUBIT_LIMIT).sum())
            run_review_pct = run_review_count / run_review_total
            run_review_required = run_review_pct > RUN_REVIEW_FRACTION

        ts_clean_cols = master_ts_df.columns.str.strip()
        ts_desc_idx = list(ts_clean_cols).index('Sample Description') if 'Sample Description' in ts_clean_cols else None
        ts_from_idx = list(ts_clean_cols).index('From [bp]') if 'From [bp]' in ts_clean_cols else None
        ts_pct_idx = list(ts_clean_cols).index('% of Total') if '% of Total' in ts_clean_cols else None

        qb_clean_cols = master_qb_df.columns.str.strip()
        qb_id_col_raw = find_qubit_id_col(master_qb_df)
        qb_id_idx = list(master_qb_df.columns).index(qb_id_col_raw) if qb_id_col_raw else None
        qb_conc_idx = list(qb_clean_cols).index('Original Sample Conc.') if 'Original Sample Conc.' in qb_clean_cols else None

        # Store measurement columns as decimals so rerun values and Round 3 averages can be written back
        if ts_pct_idx is not None:
            master_ts_df.isetitem(ts_pct_idx, pd.to_numeric(master_ts_df.iloc[:, ts_pct_idx], errors='coerce').astype(float))
        if qb_conc_idx is not None:
            master_qb_df.isetitem(qb_conc_idx, pd.to_numeric(master_qb_df.iloc[:, qb_conc_idx], errors='coerce').astype(float))

        is_rerun_mode = False
        audit_trail_log = []
        # Track every originally failed sample so %CV can be reported for all repeats
        sample_history = {s_id: {'ts': [], 'qb': []} for s_id in original_failures}
        auto_repeat_samples = set()   # Round 2 samples repeated because %CV >= 20%
        round3_samples = set()        # Samples evaluated with 3 runs (final decision)

        # Safe tracking collection without scalar extraction error risks
        master_region_mask = first_region_mask(master_ts_df, ts_desc_idx, ts_from_idx)
        for idx in master_ts_df.index[master_region_mask]:
            s_id = str(master_ts_df.loc[idx].iloc[ts_desc_idx]).strip()
            if s_id in sample_history:
                sample_history[s_id]['ts'].append(pd.to_numeric(master_ts_df.loc[idx].iloc[ts_pct_idx], errors='coerce'))

        for idx, row in master_qb_df.iterrows():
            s_id = str(row[qb_id_col_raw]).strip()
            if s_id in sample_history:
                sample_history[s_id]['qb'].append(pd.to_numeric(row['Original Sample Conc.'], errors='coerce'))

        # Ensure multi-entry logs choose only the first index safely
        for s_id in sample_history:
            if len(sample_history[s_id]['ts']) > 1:
                sample_history[s_id]['ts'] = [sample_history[s_id]['ts'][0]]
            if len(sample_history[s_id]['qb']) > 1:
                sample_history[s_id]['qb'] = [sample_history[s_id]['qb'][0]]

        # ----------------------------------------------------
        # ROUND 2 PROCESSING (RERUN 1)
        # ----------------------------------------------------
        if ts_file_2 and qb_file_2:
            is_rerun_mode = True
            raw_ts_rerun_1 = read_data_file(ts_file_2, encoding='latin1')
            raw_qb_rerun_1 = read_data_file(qb_file_2, encoding='latin1')
            
            raw_ts_rerun_1.columns = raw_ts_rerun_1.columns.str.strip()
            raw_qb_rerun_1.columns = raw_qb_rerun_1.columns.str.strip()

            raw_ts_rerun_1['From [bp]'] = pd.to_numeric(raw_ts_rerun_1['From [bp]'], errors='coerce')
            ts_rerun_1_filtered = raw_ts_rerun_1[first_region_mask(raw_ts_rerun_1, raw_ts_rerun_1.columns.get_loc('Sample Description'), raw_ts_rerun_1.columns.get_loc('From [bp]'))]
            
            for _, rerun_row in ts_rerun_1_filtered.iterrows():
                sample_id = str(rerun_row['Sample Description']).strip()
                new_pct = float(rerun_row['% of Total'])
                
                if sample_id in sample_history:
                    sample_history[sample_id]['ts'].append(new_pct)

                ts_mask = (master_ts_df.iloc[:, ts_desc_idx].astype(str).str.strip() == sample_id) &                           master_region_mask
                if ts_mask.any():
                    old_pct = master_ts_df.iloc[ts_mask, ts_pct_idx].iloc[0]
                    master_ts_df.iloc[ts_mask, ts_pct_idx] = new_pct
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File": "TapeStation",
                        "Rerun Round": "Round 2 (Rerun 1)",
                        "Parameter Updated": "% of Total (Region)",
                        "Old Value": f"{old_pct}%" if pd.notna(old_pct) else "BLANK",
                        "New Value": f"{new_pct}%"
                    })

            for _, rerun_row in raw_qb_rerun_1.dropna(subset=['Original Sample Conc.']).iterrows():
                qb_rerun_id_col = find_qubit_id_col(raw_qb_rerun_1)
                sample_id = str(rerun_row[qb_rerun_id_col]).strip()
                new_conc = float(rerun_row['Original Sample Conc.'])
                
                if sample_id in sample_history:
                    sample_history[sample_id]['qb'].append(new_conc)
                
                qb_mask = (master_qb_df.iloc[:, qb_id_idx].astype(str).str.strip() == sample_id)
                if qb_mask.any():
                    old_conc = master_qb_df.iloc[qb_mask, qb_conc_idx].iloc[0]
                    master_qb_df.iloc[qb_mask, qb_conc_idx] = new_conc
                    
                    audit_trail_log.append({
                        "Sample ID": sample_id,
                        "Instrument File": "Qubit",
                        "Rerun Round": "Round 2 (Rerun 1)",
                        "Parameter Updated": "Original Sample Conc. (ng/µL)",
                        "Old Value": f"{old_conc} ng/µL" if pd.notna(old_conc) else "BLANK",
                        "New Value": f"{new_conc} ng/µL"
                    })

            for s_id, metrics in sample_history.items():
                # Report %CV for EVERY rerun sample (repeat or no-repeat) whenever at least one instrument has 2 readings
                ts_cv = compute_cv(metrics['ts'][:2])
                qb_cv = compute_cv(metrics['qb'][:2])
                if ts_cv is None and qb_cv is None:
                    continue

                pipeline_cv_reporting[s_id] = format_cv(ts_cv, qb_cv)

                if ts_cv is not None and qb_cv is not None:
                    if ts_cv < CV_LIMIT and qb_cv < CV_LIMIT:
                        pipeline_status_overrides[s_id] = "NO REPEAT NEEDED"
                        audit_trail_log.append({
                            "Sample ID": s_id,
                            "Instrument File": "Pipeline Logic",
                            "Rerun Round": "Round 2 Analysis",
                            "Parameter Updated": "QC Status Designation",
                            "Old Value": "FAIL",
                            "New Value": f"NO REPEAT NEEDED (%CV TS: {ts_cv:.1f}%, QB: {qb_cv:.1f}%)"
                        })

                # Any %CV at or above 20% -> automatically repeated, even if the rerun values pass QC
                cv_breaches = [c for c in (ts_cv, qb_cv) if c is not None and c >= CV_LIMIT]
                if cv_breaches:
                    pipeline_status_overrides[s_id] = "FAIL"
                    auto_repeat_samples.add(s_id)
                    audit_trail_log.append({
                        "Sample ID": s_id,
                        "Instrument File": "Pipeline Logic",
                        "Rerun Round": "Round 2 Analysis",
                        "Parameter Updated": "QC Status Designation",
                        "Old Value": "FAIL",
                        "New Value": f"REPEAT - %CV at or above {CV_LIMIT:.0f}% ({format_cv(ts_cv, qb_cv)})"
                    })

        # ----------------------------------------------------
        # ROUND 3 PROCESSING (RERUN 2)
        # ----------------------------------------------------
        if ts_file_3 and qb_file_3:
            is_rerun_mode = True
            raw_ts_rerun_2 = read_data_file(ts_file_3, encoding='latin1')
            raw_qb_rerun_2 = read_data_file(qb_file_3, encoding='latin1')
            
            raw_ts_rerun_2.columns = raw_ts_rerun_2.columns.str.strip()
            raw_qb_rerun_2.columns = raw_qb_rerun_2.columns.str.strip()

            raw_ts_rerun_2['From [bp]'] = pd.to_numeric(raw_ts_rerun_2['From [bp]'], errors='coerce')
            ts_rerun_2_filtered = raw_ts_rerun_2[first_region_mask(raw_ts_rerun_2, raw_ts_rerun_2.columns.get_loc('Sample Description'), raw_ts_rerun_2.columns.get_loc('From [bp]'))]
            
            for _, rerun_row in ts_rerun_2_filtered.iterrows():
                sample_id = str(rerun_row['Sample Description']).strip()
                new_pct = float(rerun_row['% of Total'])
                
                if sample_id in sample_history and len(sample_history[sample_id]['ts']) < 3:
                    sample_history[sample_id]['ts'].append(new_pct)

            for _, rerun_row in raw_qb_rerun_2.dropna(subset=['Original Sample Conc.']).iterrows():
                qb_rerun_id_col = find_qubit_id_col(raw_qb_rerun_2)
                sample_id = str(rerun_row[qb_rerun_id_col]).strip()
                new_conc = float(rerun_row['Original Sample Conc.'])
                
                if sample_id in sample_history and len(sample_history[sample_id]['qb']) < 3:
                    sample_history[sample_id]['qb'].append(new_conc)

            for s_id, metrics in sample_history.items():
                if len(metrics['ts']) == 3 and len(metrics['qb']) == 3:
                    
                    def filter_outlier_and_average(values_list):
                        arr = np.array(values_list)
                        median = np.median(arr)
                        distances = np.abs(arr - median)
                        outlier_idx = np.argmax(distances)
                        remaining_values = np.delete(arr, outlier_idx)
                        return np.mean(remaining_values), arr[outlier_idx], remaining_values

                    final_ts_avg, ts_outlier, rem_ts = filter_outlier_and_average(metrics['ts'])
                    final_qb_avg, qb_outlier, rem_qb = filter_outlier_and_average(metrics['qb'])

                    ts_cv = (np.std(rem_ts, ddof=1) / final_ts_avg) * 100 if final_ts_avg > 0 else 0
                    qb_cv = (np.std(rem_qb, ddof=1) / final_qb_avg) * 100 if final_qb_avg > 0 else 0
                    pipeline_cv_reporting[s_id] = f"TS: {ts_cv:.1f}% | QB: {qb_cv:.1f}% (Outliers Dropped)"

                    # Final run: judge the averaged values, not the Round 2 %CV repeat
                    if pipeline_status_overrides.get(s_id) == "FAIL":
                        pipeline_status_overrides.pop(s_id)
                    round3_samples.add(s_id)

                    ts_mask = (master_ts_df.iloc[:, ts_desc_idx].astype(str).str.strip() == s_id) &                               master_region_mask
                    if ts_mask.any():
                        master_ts_df.iloc[ts_mask, ts_pct_idx] = final_ts_avg

                    qb_mask = (master_qb_df.iloc[:, qb_id_idx].astype(str).str.strip() == s_id)
                    if qb_mask.any():
                        master_qb_df.iloc[qb_mask, qb_conc_idx] = final_qb_avg

                    audit_trail_log.append({
                        "Sample ID": s_id,
                        "Instrument File": "Multi-Run Optimization",
                        "Rerun Round": "Round 3 (Outlier Purge)",
                        "Parameter Updated": "Averaged Quant Matrix",
                        "Old Value": f"TS Outlier dropped: {ts_outlier}% | QB Outlier dropped: {qb_outlier} ng/µL",
                        "New Value": f"Scrubbed TS Mean: {final_ts_avg:.2f}% | Scrubbed QB Mean: {final_qb_avg:.2f} ng/µL"
                    })

        audit_df = pd.DataFrame(audit_trail_log)
        final_df = process_data(master_ts_df, master_qb_df, pipeline_status_overrides, pipeline_cv_reporting)

        if final_df.empty:
            st.error("❌ No exact sample ID matches found in the data log parameters.")
            st.stop()

        final_df['Is Originally Failed'] = final_df['Sample Description'].isin(original_failures)
        
        def calculate_recovery_flag(row):
            if row['Is Originally Failed'] and row['QC Status'] in ['PASS', 'ABOVE UPPER LIMIT', 'NO REPEAT NEEDED']:
                return True
            return False
            
        final_df['Is Recovered'] = final_df.apply(calculate_recovery_flag, axis=1)

        def resolve_status_hierarchy(row):
            if row['Is Originally Failed'] and row['QC Status'] == 'PASS':
                return 'RECOVERED'
            return row['QC Status']
            
        final_df['QC Status'] = final_df.apply(resolve_status_hierarchy, axis=1)

        # After the third run a failing sample is a final FAIL, not a sample to repeat
        final_fail_mask = final_df['Sample Description'].isin(round3_samples) & (final_df['QC Status'] == 'FAIL')
        final_df.loc[final_fail_mask, 'QC Status'] = 'FINAL FAIL'

        # Repeat samples that have no rerun data yet still get a %CV entry so the column is never blank
        if is_rerun_mode:
            repeat_no_cv_mask = (
                (final_df['QC Status'] == 'FAIL')
                & (final_df['Calculated %CV (Reruns)'].astype(str).str.strip() == "")
            )
            final_df.loc[repeat_no_cv_mask, 'Calculated %CV (Reruns)'] = "TS: N/A | QB: N/A (Not in rerun files)"

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
        # RUN-LEVEL TECHNICAL SUPERVISOR REVIEW FLAG
        # (shown as a red bar like the sample count flag, but does NOT block processing)
        # ----------------------------------------------------
        if run_review_required:
            st.error("🚨 **Run Flagged for Technical Supervisor Review!**")
            st.error(
                f"**{run_review_count}** of **{run_review_total}** samples "
                f"(**{run_review_pct:.1%}**) have a Raw Qubit above **{RUN_REVIEW_QUBIT_LIMIT} ng/µL** "
                f"| Allowed: {RUN_REVIEW_FRACTION:.0%} or fewer."
            )
            st.info("💡 A technical supervisor must review this run before results are released.")

        # ----------------------------------------------------
        # 3. DASHBOARD SUMMARY PANEL
        # ----------------------------------------------------
        st.write("---")
        st.subheader("📊 Combined Run Analysis Summary" if is_rerun_mode else "📊 Initial Run Analysis Summary")

        samples_to_repeat_count = len(final_df[final_df['QC Status'] == 'FAIL'])
        above_3sd_count = int(
            ((final_df['QC Status'] == 'FAIL') & (final_df['Raw Qubit (ng/µL)'] > QUBIT_3SD_LIMIT)).sum()
        )

        if is_rerun_mode:
            # Rerun dashboard: only total samples and samples to repeat
            m1, m2, m3 = st.columns(3)
            m1.metric("Total Reported Samples", len(final_df))
            m2.metric("🛑 Samples to Repeat", samples_to_repeat_count)
            m3.metric("📈 Qubit 3SD Outliers (Repeat)", above_3sd_count)
        else:
            # Initial run dashboard
            m1, m2, m3, m4, m5 = st.columns(5)
            m1.metric("Total Reported Samples", len(final_df))
            m2.metric("✅ Passed QC Check", len(final_df[final_df['QC Status'] == "PASS"]))
            m3.metric("⚠️ Above Upper Limit (Total)", len(final_df[final_df['QC Status'] == "ABOVE UPPER LIMIT"]))
            m4.metric("🛑 Samples to Repeat", samples_to_repeat_count)
            m5.metric("📈 Qubit 3SD Outliers (Repeat)", above_3sd_count)

        # ----------------------------------------------------
        # 4. INTERACTIVE VIEW DROPDOWN FILTER
        # ----------------------------------------------------
        st.subheader("📋 Output Matrix Data Viewer")
        
        base_filters = ["Show All Samples", "Show Only PASS Samples", "Show Only ABOVE UPPER LIMIT Samples", "Show Only RECOVERED Samples", "Show Only Samples to Repeat"]
        if is_rerun_mode:
            base_filters.insert(4, "Show Only NO REPEAT NEEDED Samples")
            if (final_df['QC Status'] == 'FINAL FAIL').any():
                base_filters.append("Show Only FINAL FAIL Samples")
            
        status_filter = st.selectbox("Filter table view display parameters:", base_filters)
        
        if status_filter == "Show Only PASS Samples":
            filtered_display = final_df[final_df['QC Status'] == "PASS"]
        elif status_filter == "Show Only ABOVE UPPER LIMIT Samples":
            filtered_display = final_df[final_df['QC Status'] == "ABOVE UPPER LIMIT"]
        elif status_filter == "Show Only RECOVERED Samples":
            filtered_display = final_df[(final_df['QC Status'] == "RECOVERED") | (final_df['Is Recovered'] == True)]
        elif status_filter == "Show Only NO REPEAT NEEDED Samples":
            filtered_display = final_df[final_df['QC Status'] == "NO REPEAT NEEDED"]
        elif status_filter == "Show Only Samples to Repeat":
            filtered_display = final_df[final_df['QC Status'] == 'FAIL']
        elif status_filter == "Show Only FINAL FAIL Samples":
            filtered_display = final_df[final_df['QC Status'] == 'FINAL FAIL']
        else:
            filtered_display = final_df

        if not is_rerun_mode:
            filtered_display = filtered_display.drop(columns=["Calculated %CV (Reruns)"])

        def color_qc_row(row):
            styles = [''] * len(row)
            cols = list(row.index)
            
            if 'QC Status' not in cols:
                return styles
                
            qc_status = row['QC Status']
            is_recovered_flag = row.get('Is Recovered', False)
            
            status_idx = cols.index('QC Status')
            qubit_idx = cols.index('Raw Qubit (ng/µL)') if 'Raw Qubit (ng/µL)' in cols else -1
            tapestation_idx = cols.index('TapeStation % of Total') if 'TapeStation % of Total' in cols else -1
            size_idx = cols.index('Average Size [bp]') if 'Average Size [bp]' in cols else -1

            if qc_status == 'FAIL':
                styles[status_idx] = 'background-color: #ffcccc; color: #cc0000; font-weight: bold;'
            elif qc_status == 'FINAL FAIL':
                styles[status_idx] = 'background-color: #cc0000; color: #ffffff; font-weight: bold;'
            elif qc_status == 'PASS':
                styles[status_idx] = 'background-color: #ccffcc; color: #006600; font-weight: bold;'
            elif qc_status == 'RECOVERED':
                styles[status_idx] = 'background-color: #e6f7ff; color: #0050b3; font-weight: bold;'
            elif qc_status == 'NO REPEAT NEEDED':
                styles[status_idx] = 'background-color: #eaf2ff; color: #106ba3; font-weight: bold; border: 1px dashed #106ba3;'
            elif qc_status == 'ABOVE UPPER LIMIT':
                if is_recovered_flag:
                    styles[status_idx] = 'background-color: #fff2cc; color: #4a148c; font-weight: bold; border: 2px solid #4a148c;'
                else:
                    styles[status_idx] = 'background-color: #fff2cc; color: #d68100; font-weight: bold;'

            if "HALE" in selected_study:
                qubit_limit, ts_limit = 1.318, 89.69
            else:
                qubit_limit, ts_limit = 3.68, 92.89

            fail_cell_style = 'background-color: #ffcccc; color: #cc0000; font-weight: bold;'
            upper_cell_style = 'background-color: #fff2cc; color: #d68100; font-weight: bold;'

            if qc_status in ['FAIL', 'FINAL FAIL']:
                # ---- SAMPLES TO REPEAT / FINAL FAILS: highlight only the value(s) causing the failure ----
                mass_idx = cols.index('Total Regional Mass (ng in 50µL)') if 'Total Regional Mass (ng in 50µL)' in cols else -1
                region_idx = cols.index('Region Window') if 'Region Window' in cols else -1
                cv_idx = cols.index('Calculated %CV (Reruns)') if 'Calculated %CV (Reruns)' in cols else -1

                if row['Region Window'] == 'No Region Found':
                    if region_idx != -1:
                        styles[region_idx] = fail_cell_style
                else:
                    if tapestation_idx != -1 and float(row['TapeStation % of Total']) <= 60.0:
                        styles[tapestation_idx] = fail_cell_style
                    if mass_idx != -1 and float(row['Total Regional Mass (ng in 50µL)']) < MIN_TOTAL_MASS_NG:
                        styles[mass_idx] = fail_cell_style
                    # Qubit 3SD outlier -> highlight the Qubit value as the reason for rerun
                    if qubit_idx != -1 and float(row['Raw Qubit (ng/µL)']) > QUBIT_3SD_LIMIT:
                        styles[qubit_idx] = fail_cell_style

                # Highlight the %CV in red for every sample that failed again after a rerun
                if cv_idx != -1:
                    cv_text = str(row['Calculated %CV (Reruns)']).strip()
                    if cv_text and "Not in rerun files" not in cv_text:
                        styles[cv_idx] = fail_cell_style
            else:
                # ---- ALL OTHER SAMPLES: keep the existing upper-limit highlighting ----
                if qubit_idx != -1 and float(row['Raw Qubit (ng/µL)']) > qubit_limit:
                    styles[qubit_idx] = upper_cell_style
                if tapestation_idx != -1 and float(row['TapeStation % of Total']) > ts_limit:
                    styles[tapestation_idx] = upper_cell_style
                if size_idx != -1 and float(row['Average Size [bp]']) > 350.0:
                    styles[size_idx] = upper_cell_style

            return styles

        st.dataframe(
            filtered_display.style.apply(color_qc_row, axis=1), 
            use_container_width=True,
            height=500
        )

        # ----------------------------------------------------
        # 5. EXPORT FORMAT GENERATION SYSTEM
        # ----------------------------------------------------
        st.write("---")
        st.subheader("📥 Download Modified Instrument Files & Audit Trail Logs")
        
        ts_buffer = io.StringIO()
        master_ts_df.to_csv(ts_buffer, index=False)
        ts_csv_text = ts_buffer.getvalue().replace("Conc. [pg/Âµl]", "Conc. [pg/µl]")
        ts_csv_bytes = ts_csv_text.encode('latin1', errors='ignore')

        qb_buffer = io.StringIO()
        master_qb_df.to_csv(qb_buffer, index=False)
        qb_csv_bytes = qb_buffer.getvalue().encode('latin1', errors='ignore')

        audit_csv_bytes = b""
        if is_rerun_mode and not audit_df.empty:
            audit_buffer = io.StringIO()
            audit_df.to_csv(audit_buffer, index=False)
            audit_csv_bytes = audit_buffer.getvalue().encode('utf-8')

        dl_col1, dl_col2, dl_col3 = st.columns(3)
        with dl_col1:
            st.download_button(label="📥 Download Updated TapeStation File", data=ts_csv_bytes, file_name="updated_tapestation_report.csv", mime="text/csv")
        with dl_col2:
            st.download_button(label="📥 Download Updated Qubit File", data=qb_csv_bytes, file_name="updated_qubit_report.csv", mime="text/csv")
        with dl_col3:
            if is_rerun_mode and audit_csv_bytes != b"":
                st.download_button(label="📜 Download Modification Trace Log", data=audit_csv_bytes, file_name="rerun_modification_audit_log.csv", mime="text/csv")
            else:
                st.button("📜 Download Modification Trace Log", disabled=True)

        st.success("✅ Output matrices and validation trace files generated successfully.")

    except Exception as e:
        st.error(f"Processing Error: {e}")
