import pandas as pd
import os

from modules.bedops import bedops_coincidence, bedops_stopping

# Compare the coordinates
def coincidence_counter (df1, df2):
    coincidence = 0
    for _, row in df2.iterrows():
        if row["sstart"] in list(df1["sstart"].values) and row["send"] in list(df1["send"].values):
            coincidence += 1
    return coincidence

def stopping_main(data_df1, data_df2):
    data_df1_plus = data_df1[data_df1["sstrand"] == "plus"].copy()
    data_df1_minus = data_df1[data_df1["sstrand"] == "minus"].copy()
    data_df2_plus = data_df2[data_df2["sstrand"] == "plus"].copy()
    data_df2_minus = data_df2[data_df2["sstrand"] == "minus"].copy()

    coincidence_plus = coincidence_counter(data_df1_plus, data_df2_plus)
    coincidence_minus = coincidence_counter(data_df1_minus, data_df2_minus)

    total_coincidence = coincidence_plus + coincidence_minus
    perc_coincidence = total_coincidence / data_df2.shape[0] * 100
    print(f"\t\t\t- Coincidence with last corrected sequences:\n",
          f"\t\t\t\t- {total_coincidence}/{data_df2.shape[0]} - {perc_coincidence:.2f}%")

    if total_coincidence == data_df2.shape[0]:
        print(f"\t\t\t\t- TRUE")
        return True
    else:  # If the the coincidence is not the 100%
        print(f"\t\t\t\t- FALSE")
        return False
    
def stopping_bedops(new_data, old_data, folder_path, genome_fasta):
    new_data_plus = new_data[new_data["sstrand"] == "plus"].copy()
    new_data_minus = new_data[new_data["sstrand"] == "minus"].copy()
    old_data_plus = old_data[old_data["sstrand"] == "plus"].copy()
    old_data_minus = old_data[old_data["sstrand"] == "minus"].copy()

    # Prepare paths
    plus_path = os.path.join(folder_path, "plus")
    minus_path = os.path.join(folder_path, "minus")
    os.makedirs(plus_path, exist_ok=True)
    os.makedirs(minus_path, exist_ok=True)

    # -----------------------------------------------------------------------------
    ## Call BEDOPS on plus
    coincidence_plus, data_plus = bedops_stopping(new_data_plus, old_data_plus, plus_path, "plus", genome_fasta)

    # Call BEDOPS on minus. Special case, because BEDOPS reads the coordinates like the "+" strand.
    ## First modify the coordinates.
    new_data_minus[["sstart", "send"]] = new_data_minus[["send", "sstart"]]
    old_data_minus[["sstart", "send"]] = old_data_minus[["send", "sstart"]]

    # -----------------------------------------------------------------------------
## And now call BEDOPS on minus
    coincidence_minus, data_minus = bedops_coincidence(new_data_minus, old_data_minus, minus_path. "minus", genome_fasta)
    if not data_minus.empty:  # If the data frame is not empty
        data_minus[["sstart", "send"]] = data_minus[["send", "sstart"]]  # restore "data_minus" coordinates

    recapture_data = pd.concat([data_plus, data_minus], ignore_index=True)  # joins both Data Frames
    if not recapture_data.empty:  # If the data frame is not empty
        recapture_data = recapture_data.sort_values(by=["sstrand", "sseqid", "sstart"])  # Sort the data frame by the start coordinate
 
    # -----------------------------------------------------------------------------
    ## math part
    total_coincidence = coincidence_plus + coincidence_minus
    perc_coincidence = total_coincidence / old_data.shape[0] * 100
    print(f"\t\t\t- Coincidence (BEDOPS version):\n",
          f"\t\t\t\t- {total_coincidence}/{old_data.shape[0]} - {perc_coincidence:.2f}%")
    
    # -----------------------------------------------------------------------------
    ## Now recapture data

    
    if total_coincidence == old_data.shape[0]:
        print(f"\t\t\t\t- TRUE")
        return True, recapture_data
    else:  # If the the coincidence is not the 100%
        print(f"\t\t\t\t- FALSE")
        return False, recapture_data

