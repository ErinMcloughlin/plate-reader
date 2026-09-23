import streamlit as st
import pandas as pd

st.set_page_config(page_title="Plate Log Reader", page_icon="📊", layout="wide")

st.title("📊 License Plate File Reader")
st.write("Upload your structural log files (CSV, TXT, or Excel) containing extracted license plate data.")

# File uploader widget accepting tabular/text data files
uploaded_file = st.file_uploader(
    "Choose a data file...", 
    type=["csv", "txt", "xlsx"]
)

if uploaded_file is not None:
    try:
        # 1. Parse file based on extension
        file_details = {"FileName": uploaded_file.name, "FileType": uploaded_file.type}
        
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        elif uploaded_file.name.endswith('.xlsx'):
            df = pd.read_excel(uploaded_file)
        else:
            # For raw .txt files, read line by line into a DataFrame
            lines = [line.decode("utf-8").strip() for line in uploaded_file.readlines()]
            df = pd.DataFrame(lines, columns=["Raw Plate Record"])

        st.success(f"Successfully loaded: **{file_details['FileName']}**")
        
        # 2. Add an interactive search/filter widget
        st.subheader("🔍 Filter & Analyze Logs")
        search_query = st.text_input("Search for a specific license plate string:")
        
        # Apply filter if user types something
        if search_query:
            # Case-insensitive partial matching across all text columns
            mask = df.astype(str).apply(lambda x: x.str.contains(search_query, case=False)).any(axis=1)
            filtered_df = df[mask]
        else:
            filtered_df = df

        # 3. Display data summary metric card
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Total Records Parsed", len(df))
        with col2:
            st.metric("Filtered Records Displayed", len(filtered_df))

        # 4. Display the resulting interactive spreadsheet grid
        st.subheader("📋 Log Data Viewer")
        st.dataframe(filtered_df, use_container_width=True)
        
    except Exception as e:
        st.error(f"Error parsing file structure: {e}")
        st.warning("Please ensure the file is a properly formatted CSV, Excel sheet, or plaintext log file.")
