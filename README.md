# Fasta-Data---python-
Read fasta file using python and plot graphs
The code is divided into following major components
•	Reading input data 
•	Writing extracted data to respective module folders 
•	Plotting graphs 
We are using read_data() function to read the input file. This function makes use of starting characters of each line in the input file. If the line starts with ‘##’ we simple skip that line as it only states software name and version. If the line starts with ‘>>’ then this either represents start of a module or end of a module. We keep track of each module start and end using module_check variable, initialized to a 0. module_check = 0 represent module start/continue and module_check=1 represents module end. For each module, its respective column names are saved in data_headers dictionary, while the values are saved in data_values dictionary. The filter value for each module is saved in data_filter dictionary. In each of these dictionaries we use module name as the key to identify module under consideration. These are also filled based on the first check we discussed above i.e. check on starting characters of each line. If the line starts with ‘#’ (not with ‘##’) we save it in data_headers and if all of above initial line character check conditions fails, then we save that line in data_values. 
Now that we have separated data for each module, we save it in their respective folders. After this we combine the above created dictionaries into a single dictionary (data_extracted), where module_name is the key, and its value is a list, having first element as a 2d list of all columns and the next element as a 2d list of all values. 
After this we use get_frames() function to create pandas dataframes for each module, which are then used by plot_graphs() function to plot graphs.
The program starts by taking input from user which are passed to read_data() function. read_data() takes 2 input arguments
o	Input file name 
o	Output folder 
And returns a dictionary with extracted data. This extracted data is then passed to get_frames() function which returns a dictionary of dataframes identified by module_name as keys. 
This dictionary is then passed to plot_graphs() function which plot various graphs for different modules.   
