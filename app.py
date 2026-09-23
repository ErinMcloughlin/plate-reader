import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="NGS Library Multi-Run QC", page_icon="🧬", layout="wide")

st.title("🧬 Multi-Run NGS Yield & Rerun Validation Hub")
st.write("Process initial plate parameters, identify structural dropouts, and track failure recovery runs.")

# Constants
TOTAL_VOLUME_UL = 50.0

def process_and_qc(ts_file, qb_file):
    """Helper function to parse, join, calculate, and QC a pair of TS and Qubit files."""
    # 1. Read and clean TapeStation Data
    ts_df = pd.read_csv(ts_file, encoding='latin1')
    ts_df.columns = ts_df.columns.str.strip()
    ts_df['Sample Description'] = ts_df['Sample Description'].astype(str).str.strip()
    
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
        return None, "Qubit file missing an identifier column like 'Sample Description' or 'Sample Name'."
        
    qb_df = qb_df.dropna(subset=['Original Sample Conc.'])

    # 3. Perform direct text inner merge matching keys
    merged_df = pd.merge(ts_filtered, qb_df, on='Sample Description', how='inner')
    if merged_df.empty:
        return None, "No exact alphanumeric text matches found between these two files."

    # 4. Perform the mass balance calculations
    merged_df['% of Total'] = pd.to_numeric(merged_df['% of Total'], errors='coerce')
    merged_df['Original Sample Conc.'] = pd.to_numeric(merged_df['Original Sample Conc.'], errors='coerce')
    
    merged_df['Calculated Region (ng/µL)'] = merged_df['Original Sample Conc.'] * (merged_df['% of Total'] / 100.0)
    merged_df['Total Regional Mass (ng in 50µL)'] = merged_df['Calculated Region (ng/µL)'] * TOTAL_VOLUME_UL

    # 5. Apply Strict QC Metrics Logic Check
    merged_df['QC Status'] = merged_df.apply(
        lambda r: "FAIL" if (r['% of Total'] <= 60.0 or r['Total Regional Mass (ng in 50µL)'] <= 10.0) else "PASS", 
        axis=1
    )

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
    return display_df, None

# Formatting helper function for table display cells
def color_qc(val):
    if val == 'FAIL':
        return 'background-color: #ffcccc; color: #cc0000; font-weight: bold'
    elif val == 'PASS':
        return 'background-color: #ccffcc; color: #006600; font-weight: bold'
    return ''

# ==========================================
# STAGE 1: INITIAL RUN PROCESSING
# ==========================================
st.header("1️⃣ Step 1: Process Initial Plate Baseline Run")
c1, c2 = st.columns(2)
with c1:
    ts_file_1 = st.file_uploader("Upload Initial TapeStation Region CSV", type=["csv"], key="ts1")
with c2:
    qb_file_1 = st.file_uploader("Upload Initial Qubit Run CSV", type=["csv"], key="qb1")

initial_df = None
initial_failures = set()

if ts_file_1 and qb_file_1:
    df1, err1 = process_and_qc(ts_file_1, qb_file_1)
    if err1:
        st.error(f"Error processing Step 1: {err1}")
    else:
        initial_df = df1
        initial_failures = set(df1[df1['QC Status'] == "FAIL"]['Sample Description'].tolist())
        
        # Summary Dashboard Display Panel
        t_count = len(df1)
        p_count = len(df1[df1['QC Status'] == "PASS"])
        f_count = len(initial_failures)
        
        m1, m2, m3 = st.columns(3)
        m1.metric("Total Baseline Samples", t_count)
        m2.metric("✅ Passed Baseline", p_count)
        m3.metric("❌ Failed Baseline (Need Rerun)", f_count, delta=f"-{f_count}" if f_count > 0 else None, delta_color="inverse")
        
        with st.expander("👁️ View Full Initial Baseline Data Sheet Grid"):
            st.dataframe(df1.style.map(color_qc, subset=['QC Status']), use_container_width=True)

# ==========================================
# STAGE 2: FAILURE RERUN PROCESSING
# ==========================================
st.write("---")
st.header("2️⃣ Step 2: Evaluate Failure Rerun Recovery Data")

if not initial_df:
    st.info("💡 Please complete Step 1 first to identify baseline processing failures.")
elif len(initial_failures) == 0:
    st.success("🎉 Outstanding! Your initial run had 0 failures. No second stage testing required.")
else:
    st.warning(f"📋 **Action Required:** Upload rerun data matrices specifically assessing these {len(initial_failures)} failed sample IDs: `{list(initial_failures)}`")
    
    c3, c4 = st.columns(2)
    with c3:
        ts_file_2 = st.file_uploader("Upload Rerun TapeStation Region CSV", type=["csv"], key="ts2")
    with c4:
        qb_file_2 = st.file_uploader("Upload Rerun Qubit Run CSV", type=["csv"], key="qb2")
        
    if ts_file_2 and qb_file_2:
        df2, err2 = process_and_qc(ts_file_2, qb_file_2)
        if err2:
            st.error(f"Error processing Step 2: {err2}")
        else:
            # Build comparative analysis dashboard tracking metrics shifts
            tracking_records = []
            
            for idx, rerun_row in df2.iterrows():
                s_id = rerun_row['Sample Description']
                
                # Retrieve matching original run history entry metrics
                orig_row = initial_df[initial_df['Sample Description'] == s_id]
                
                if not orig_row.empty:
                    orig_pct = orig_row.iloc[0]['TapeStation % of Total']
                    orig_mass = orig_row.iloc[0]['Total Regional Mass (ng in 50µL)']
                    
                    tracking_records.append({
                        "Sample Description": s_id,
                        "Orig % of Total": f"{orig_pct}%",
                        "Rerun % of Total": f"{rerun_row['TapeStation % of Total']}%",
                        "Orig Mass (ng)": orig_mass,
                        "Rerun Mass (ng)": rerun_row['Total Regional Mass (ng in 50µL)'],
                        "Updated QC Status": rerun_row['QC Status']
                    })
                else:
                    # Sample uploaded in rerun wasn't flagged or present in initial failure array scope
                    tracking_records.append({
                        "Sample Description": s_id,
                        "Orig % of Total": "N/A",
                        "Rerun % of Total": f"{rerun_row['TapeStation % of Total']}%",
                        "Orig Mass (ng)": "N/A",
                        "Rerun Mass (ng)": rerun_row['Total Regional Mass (ng in 50µL)'],
                        "Updated QC Status": rerun_row['QC Status']
                    })
                    
            comparison_df = pd.DataFrame(tracking_records)
            
            recovered_count = len(comparison_df[comparison_df['Updated QC Status'] == "PASS"])
            still_failing = len(comparison_df[comparison_df['Updated QC Status'] == "FAIL"])
            
            st.subheader("📊 Rerun Recovery Progress Analysis Dashboard")
            rc1, rc2 = st.columns(2)
            rc1.metric("🚀 Successfully Recovered (Now PASS)", recovered_count)
            rc2.metric("⚠️ Remaining Failures (Still FAIL)", still_failing, delta=f"+{still_failing}" if still_failing > 0 else None, delta_color="inverse")
            
            # Render interactive recovery grid table layout panel
            st.subheader("📋 Detailed Rerun Tracking Progress Log Table")
            st.dataframe(comparison_df.style.map(color_qc, subset=['Updated QC Status']), use_container_width=True)
            
            # Generate export data package downloads buffers triggers
            csv_buffer = io.StringIO()
            comparison_df.to_csv(csv_buffer, index=False)
            csv_data = csv_buffer.getvalue()
            
            st.download_button(
                label="📥 Download Failure Recovery Progress Report CSV",
                data=csv_data,
                file_name="rerun_recovery_validation_report.csv",
                mime="text/csv"
            )

