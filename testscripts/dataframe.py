import pandas as pd
import re

# Read the content of the output file
with open("tempfile.txt", "r") as file:
    lines = file.readlines()

# Extract sequence
sequence = ""
for line in lines:
    line = line.strip()
    if line.startswith(('U', 'C', 'T', 'A')):  # Sequence line
        sequence = line.split()[0]  # Extract sequence
        break  # Exit the loop after finding the sequence

# Extract structure and free energy
data = []
for line in lines:
    line = line.strip()
    if line and not line.startswith(('U', 'C', 'T', 'A')):  # Skip empty lines and sequence line
        structure, free_energy = line.split()  # Split structure and free energy
        data.append((structure, float(free_energy)))

# Create a DataFrame without the 'Sequence' column
df = pd.DataFrame(data, columns=["Structure", "Minimum Free Energy"])

# Define a regular expression pattern to match structures with only one set of open and close brackets
pattern = r'\.*[(.]*\.*[).]*\.*'

# Filter out structures matching the pattern
df_filtered = df[df['Structure'].str.fullmatch(pattern)]

# Display the filtered DataFrame
print(df_filtered)    
