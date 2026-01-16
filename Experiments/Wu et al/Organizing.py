import os

# Define the possible values
N_values = ["N100", "N250", "N500"]
spec_values = ["spec1", "spec2", "spec3"]

# Loop through each configuration
for N in N_values:
    for spec in spec_values:
        folder_name = f"{N}_{spec}"
        
        # Check if the folder exists
        if os.path.isdir(folder_name):
            # Get all files in the folder
            files = os.listdir(folder_name)
            files = [file for file in files if file.endswith(".csv")]
            
            # Sort files to ensure correct ordering
            files.sort(key=lambda x: int(x.split('_')[2][7:]))
            
            # Rename files
            for index, file_name in enumerate(files):
                new_name = f"dataset_{index + 1}.csv"
                os.rename(os.path.join(folder_name, file_name), os.path.join(folder_name, new_name))

print("Files renamed in each folder.")

