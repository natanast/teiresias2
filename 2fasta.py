import pandas as pd
import sys

# define filename and columns
my_file = 'C:/Users/Nina/ptyxiaki/master-thesis/data/CLL_DB_data_aligned_62K.xlsx'
name_col = 'Sequence ID'
seq_col = 'aa sequence aligned' 
output_filename = 'C:/Users/Nina/ptyxiaki/master-thesis/data/62K.fa'

# Loading the xlsx file
sequence_list = pd.read_excel(my_file, index_col=None) 
# Converting column to a list
mylist_names = sequence_list[name_col].tolist()
mylist = sequence_list[seq_col].tolist()

# writing file opening
orig_stdout = sys.stdout
f = open(output_filename, 'w')
sys.stdout = f

# iterating through the variables
for line, i in zip(mylist, mylist_names):
    line = line.strip().split(',')
    header = '>' + '_' + i
    print(header)
    print(*line, sep='\n')

#writing and closing the file
sys.stdout = orig_stdout
f.close()