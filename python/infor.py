"""
This module is used to handle the input file so that to obtain the basic information
about the result functional, namely:

the name of the input file, output file, project name etc.;	

functional type, which is indicated by four possible numbers:	
1  LDA
2  GGA
3  META
4  LAP

the number of functional derivatives, used for debugging;

the variable list which contains the non-equivalent variables, this is characterized
by the keywords of "rho, grho, tau, laplacian" in the first line of input maple file;

whether to use the "full var" for the 1st order of functional derivatives, this is 
characterized by "full" in the first line;

The control information for processing the functionals;

lambda-dependent functionals
"""
import re
import sys
import os

__author__  = "Fenglai Liu"
__date__    = "Oct, 2010"

# functional related infor
functional_type = 0
non_equi_var_list = []
use_full_var = 0
use_lambda   = 0
lambda_deriv = 0
lambda_value = "  "
project      = "  "
input_file   = "  "
output_file  = "  "

# job related infor
in_pointwise = 0        # calculate the functional in pointwise way 
functional_order = 1    # in default we only calculate the first order derivatives


#########################################################
#
#  set project name
# 
#########################################################
def set_project (fname):

	global project
	global output_file
	global input_file

	dirname, filename = os.path.split(fname)
	name, suffix = os.path.splitext(filename)
	input_file   = filename

	# form the output file name etc. Pay attention to the lambda
	set_lambda(fname)
	if lambda_deriv:
		new_name     = name + "_lambda"
		project      = new_name 
		output_file  = new_name + ".f"
	else:
		project      = name 
		output_file  = name + ".f"


#########################################################
#
#  set functional type
# 
#########################################################
def set_functional_type (fname):

	# read in first line
	name = open(fname, "r")
	line = name.readline()
	global functional_type

	# get the number of variable
	if re.search("(?i)LDA", line) is not None:
		functional_type = 1
	elif re.search("(?i)GGA", line) is not None:
		functional_type = 2
	elif re.search("(?i)META", line) is not None:
		functional_type = 3
	elif re.search("(?i)LAP", line) is not None:
		functional_type = 4
	else:
		print "Illegal type of functional! Please check the first line!\n"
		sys.exit()

	name.close()

def get_type():

	if functional_type > 0 and functional_type < 5:
		return functional_type
	else:
		print "functional_type should be within 1-4!\n"
		sys.exit()
		

#########################################################
#
#  set functional derivatives order
#  some times we only need to calculate 1st, or up to
#  2ed order functional derivatives
# 
#########################################################
def set_func_deriv(fname):

	# read in first line
	name = open(fname, "r")
	global functional_order

	has_key = 0
	while True:
		line = name.readline()
		if not line: break 
		if re.search(r"(?i)\bfunctional_order\b", line) is not None: 
			has_key = 1
			break

	if has_key:
		line = name.readline()
		if re.search("1\s", line) is not None:
			functional_order = 1
		elif re.search("2\s", line) is not None:
			functional_order = 2
		elif re.search("3\s", line) is not None:
			functional_order = 3

	name.close()

#########################################################
#
# Whether to calculate the functional and functional
# derivatives in pointwise way?
# 
#########################################################
def set_pointwise(fname):

	# read in first line
	name = open(fname, "r")
	global in_pointwise

	has_key = 0
	while True:
		line = name.readline()
		if not line: break 
		if re.search(r"(?i)\bpointwise\b", line) is not None: 
			has_key = 1
			break

	if has_key:
		in_pointwise = 1

	name.close()

#########################################################
#
#  whether we use the lambda dependent functional?
#  whether we derive the functional derivatives with
#  respect to the lambda?
# 
#########################################################
def set_lambda(fname):

	# read in first line
	name = open(fname, "r")
	line = name.readline()
	global lambda_deriv
	global use_lambda
	global lambda_value

	# get the lambda, do lambda expression or lambda derivatives?
	if re.search(r"\slambda_deriv\s", line) is not None:
		lambda_deriv = 1
	elif re.search(r"\slambda\s", line) is not None:
		use_lambda = 1
	elif re.search(r"(?i)lambda\b", line) is not None and re.search(r"lambda_deriv", line) is not None:
		print "you can only set either lambda or lambda_deriv\n"
		print "they are corresponding to different functionals, so each time we only deal with one\n"
		sys.exit()

	# for each of the case, we all need the lambda value
	# it will be defined in the result Fortran file
	if use_lambda or lambda_deriv:
		has_key = 0
		while True:
			line = name.readline()
			if not line: break 
			if re.search(r"(?i)lambda\b", line) is not None and re.search(r"=", line) is not None:
				has_key = 1
				break
		if has_key:
			value = line.split("=")
			print value
			lambda_value = value[1]
			# finally check whether it's a fortran double precision type data
			if re.search(r"(?i)D", lambda_value) is None:
				print "# you should have D symbol in the lambda_value expression, since it will be directly used in Fortran codes!!\n"
				sys.exit()
		else:
			print "why we have lambda in the first line, but no lambda value can be found??\n"
			print "please check the maple file, and use:\n"
			print "# lambda = .... to set up the lambda value\n"
			sys.exit()

	name.close()

#########################################################
#
#  set non-equivalent varibles list
# 
#########################################################

def set_non_equivalent_var_list(fname):

	# read in first line
	name = open(fname, "r")
	line = name.readline()
	global non_equi_var_list

	if re.search("(?i)\srho\s", line) is not None:
		non_equi_var_list.append("RA")
		non_equi_var_list.append("RB")
	if re.search("(?i)\sgrho\s", line) is not None:
		non_equi_var_list.append("GAA")
		non_equi_var_list.append("GBB")
		non_equi_var_list.append("GAB")
	if re.search("(?i)\stau\s", line) is not None:
		non_equi_var_list.append("TA")
		non_equi_var_list.append("TB")
	if re.search("(?i)\slaplacian\s", line) is not None:
		non_equi_var_list.append("LA")
		non_equi_var_list.append("LB")

	name.close()

def get_equi_list():
	
	return non_equi_var_list


#########################################################
#
#  set full varible so that to judge whether we use
#  all of the variables in the 1st functional derivatives
#  calculation
# 
#########################################################
def set_full_var(fname):

	# read in first line
	name = open(fname, "r")
	line = name.readline()
	global use_full_var

	if re.search("(?i)\sfull\s", line) is not None:
		use_full_var = 1

	name.close()

def use_full_var_in_1st_derivatives():

	if use_full_var == 0 or use_full_var == 1:
		return use_full_var
	else:
		print "the number in the full_var is wrong. Should be set to 0 or 1\n"
		sys.exit()


#########################################################
#
#  get functional derivatives number according to the 
#  functional type, for debugging use
# 
#########################################################

def get_deriv_number (order):

	# get the basical var number
	if functional_type == 1:
		number = 2
	elif functional_type == 2:
		number = 5
	elif functional_type == 3:
		number = 7
	else:
		number = 9

	# result	
	if order == 1:
		return number
	elif order == 2:
		return ((number+1)*number)/2
	elif order == 3:
		return ((number+2)*(number+1)*number)/6
	else:
		print "Illegal type of order! It should be within 1-3!\n"
		sys.exit()




