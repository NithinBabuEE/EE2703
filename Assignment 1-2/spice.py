
#																				SPICE PROGRAM - PART 1
#																				NITHIN BABU [EE18B021]	
# Importing the necessary modules.
from sys import argv, exit

# The program will throw in an error if there isn't exactly 2 arguments in the commandline.
if len(argv)!=2:
	print("Please provide the correct 2 arguments in the commandline.")
	exit()

# Assigning constants variables to .circuit and .end 
CIRCUIT = ".circuit"
END = ".end"

try:

# Opening the file mentioned in the commandline.
	with open(argv[1]) as f:
		lines = f.readlines()

# These are parameters to check the errors in the file format. 		
		start = -1; start_check = -1; end = -2; end_check = -1 

# The program will traverse through the file and take out only the required part.
		for line in lines:
			if CIRCUIT == line[:len(CIRCUIT)]:
				start = lines.index(line)
				start_check = 0;

			elif END == line[:len(END)]:
				end = lines.index(line)
				end_check = 0;

# The program will throw in an error if the circuit definition format is not proper.		
		if start >= end or start_check == -1 or end_check == -1:
			print("Invalid circuit definition.")
			exit()
			
# These lines of code will reverse and print the required result.
		for lines in reversed(lines[start+1:end]):
			c=reversed(lines.split('#')[0].split())	
			c = "  ".join(c)
			print(c)

# The program will throw in this error if the name of the netlist file is not proper 
# or if the netlist file is not found in the same directory as the program.

except FileNotFoundError:
	print("Invalid File.")
	exit()