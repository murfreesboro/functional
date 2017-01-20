"""
This is some python script used to generate the Fortran 
codes for the giving maple file so that to produce the 
functional and functional derivatives up to the third 
order.

How to use this program?
python main input_maple_file
or 
python main input_maple_file control_file
"""
import sys
import os
import infor
import generation
import head_print
import maple
import fortran

__author__  = "Fenglai Liu"
__version__ = "0.0"
__date__    = "Oct, 2010"

# firstly, check the argv list
has_maple_input = 0
has_control_infor = 0
if len(sys.argv) == 1:
	print "We need at least one argument, which is the input maple file!\n"
	sys.exit()
elif len(sys.argv) == 2:
	has_maple_input = 1
elif len(sys.argv) == 3:
	has_maple_input = 1
	has_control_infor = 1
else:
	print "Wrong argv list! We only support two arguments! Please check it!\n"
	sys.exit()

# obtain the maple input file name, and get the corresponding functional information
if has_maple_input:
	file_name = sys.argv[1]
	if os.path.exists(file_name):
		infor.set_functional_type(file_name)
		infor.set_non_equivalent_var_list(file_name)
		infor.set_full_var(file_name)
		infor.set_project(file_name)
		if os.path.exists(infor.output_file):
			os.remove(infor.output_file)
	else:
		print "We can not find the input maple file!\n"
		sys.exit()

# obtain the control information, all of the information also write into infor module
if has_control_infor:
	file_name = sys.argv[2]
	if os.path.exists(file_name):
		infor.set_func_deriv(file_name)
		infor.set_pointwise(file_name)
	else:
		print "We can not find the input control file!\n"
		sys.exit()


# generate the D1F etc.
generation.set_D1F()
if infor.functional_order >= 2:
	generation.set_D2F()
if infor.functional_order >= 3:
	generation.set_D3F()

# head printing
head_print.print_ref(infor.input_file, infor.output_file)
if not infor.in_pointwise:
	head_print.print_main_function(infor.output_file)

# begin to generate the maple codes
routines = ["close", "open_RA_RB", "open_RA", "open_RB"]
for routine in routines:

	# first, generate the maple codes
	if os.path.exists("tmp_tmp_fortran.F"):
		os.system("rm -rf tmp_tmp_fortran.F")
	maple.maple_code_gen(routine)

	# execute the maple file so that to produce the fortran codes
	os.system("maple -q -w0 tmp_tmp.mpl")

	# print the work function
	if infor.in_pointwise:
		head_print.print_work_function_pointwise(infor.output_file, routine)
	else:
		head_print.print_work_function(infor.output_file, routine)

	# copy the fortran codes from tmp file into the result file
	fortran.form_result_file(infor.output_file)

	# finally, clear the tmp files
	if os.path.exists("tmp_tmp.mpl"):
		os.system("rm -rf tmp_tmp.mpl")
	if os.path.exists("tmp_tmp_fortran.F"):
		os.system("rm -rf tmp_tmp_fortran.F")




