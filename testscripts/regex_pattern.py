# 04/07 Trying to work out which regex pattern will match the structure I want
# essentially want only one set of open and close brackets 
# next step: figure out how to match either ( or . (any number) 
import re

sequence = ".((((....(((((((..((..(((......)))))..)))))))......))))............."
sequence_test = "...((.((..)).)..."
pattern_test = r'\.*\(*[.()]*\)*\.*'

pattern = r'\.*[(.]*\.*[).]*\.*'
pattern2 = r'[^)]*[^(]*'

if re.fullmatch(pattern_test, sequence):
    print("The sequence matches the pattern.")
else:
    print("The sequence does not match the pattern.")

