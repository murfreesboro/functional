"""
This module is used to form the Fortran codes from the tmp maple file
"""
import re
import sys
import os

__author__  = "Fenglai Liu"
__version__ = "0.0"
__date__    = "Oct, 2010"



#########################################################
#
#  this function is used to copy the codes from tmp 
#  fortran into the result Fortran file. At the same
#  time, we get the result expression and replace the
#  "%%%%%%%" line with the result
# 
#########################################################
def form_result_file(output_file):

	tmp_fortran = open("tmp_tmp_fortran.F", "r")
	outfile = open(output_file, "a")

	result = "   "
	while True:
		line = tmp_fortran.readline()
		if not line: break
		if re.search(r"######", line) is not None:
			string = line.split("######")
			name = string[1]
			name.strip()

			# now we need to calculate the variable order
			order = 0
			if re.search(r"D1F", name) is not None:
				order = 1
			elif re.search(r"D2F", name) is not None:
				order = 2
			elif re.search(r"D3F", name) is not None:
				order = 3

			# now we calculate position variable
			# and print it out
			# this is only for derivatives
			if order > 0:
				begin = name.find("ID")
				end   = name.find("_POS")
				var   = name[begin:end]
				varPos= var + "_POS"
				varay = "D1VARS"
				if order == 2:
					varay = "D2VARS"
				elif order == 3:
					varay = "D3VARS"
				new_line = varPos + "=" + varay + "(" + var +")" 
				outfile.write("      ")
				outfile.write(new_line)
				outfile.write("\n")

			# now print out the code
			new_line = name + "=" + name + "+" + result
			outfile.write("     ")
			outfile.write(new_line)
			outfile.write("\n")
		elif re.search(r"^     #", line) is not None:
			new_line = line.replace("     #", "     &")
			outfile.write(new_line)
		else:
			outfile.write(line)
			if re.search("=", line) is not None: 
				result = get_expression(line)

	tmp_fortran.close()
	outfile.close()


#########################################################
# get the expression at the "=" line
#########################################################
def get_expression(code):

	string = code.split("=")
	result = string[0]
	result.strip()

	# check the final result, we require the final result should not be empty
	if re.search("\w+", result) is None:
		print "something wrong happended in getting the expression. Please check it!"
		print result
		sys.exit()

	return result





