"""
This module is used to generate the maple codes for the given functional and
functional derivatives
"""
import re
import sys
import os
import infor
import generation
import var

__author__  = "Fenglai Liu"
__version__ = "0.0"
__date__    = "Oct, 2010"




#########################################################
#
#  the main function for this module
# 
#########################################################
def maple_code_gen (ftype):


	# create tmp maple file, get the name of functional
	tmp_maple = open("tmp_tmp.mpl", "w")
	name = get_functional_name(ftype)

	# copy the original input maple into tmp file
	infile = open(infor.input_file, "r")
	has_name = 0
	while True:
		line = infile.readline()
		if not line: break
		tmp_maple.write(line)
		if re.search("(?i)%s"%name, line) is not None:
			has_name = 1
	infile.close()

	# check that whether we have the functional name in the input file
	if not has_name:
		print "why we can not find the functional name in the input maple file???\n"
		sys.exit()

	# print out the general subroutines in maple
	tmp_maple.write("# define procedure Fortran\n")
	tmp_maple.write("Fortran:=proc(a)\n")
	tmp_maple.write("   codegen[fortran](evalf(a),precision=double,filename=\"tmp_tmp_fortran.F\",optimized,mode=double):\n")
	tmp_maple.write("end proc:\n")
	tmp_maple.write("\n")
	tmp_maple.write("# define the text writing for the fortran file\n")
	tmp_maple.write("Text:=proc(str)\n")
	tmp_maple.write("   local fd:\n")
	tmp_maple.write("   fd:=fopen(\"tmp_tmp_fortran.F\",APPEND):\n")
	tmp_maple.write("   fprintf(fd,\"      %s\\n\",str):\n")
	tmp_maple.write("   fclose(fd):\n")
	tmp_maple.write("end proc:\n")
	tmp_maple.write("\n")
	tmp_maple.write("# assign the precision for the output \n")
	tmp_maple.write("Digits := 21:\n")

	# printing the head of the codes, and get the routine type
	routine_type = "close"
	if ftype == "open_RA_RB":
		tmp_maple.write("Text(\"IF (RA.GT.TOL .AND. RB.GT.TOL) THEN\");\n")
		routine_type = "open"
	elif ftype == "open_RA":
		tmp_maple.write("Text(\"ELSE IF (RA.GT.TOL .AND. RB.LE.TOL) THEN\");\n")
		routine_type = "open"
	elif ftype == "open_RB":
		tmp_maple.write("Text(\"ELSE IF (RA.LE.TOL .AND. RB.GT.TOL) THEN\");\n")
		routine_type = "open"
	else:
		if ftype != "close":
			print "The routine type is not correct???\n"
			sys.exit()
	
	# generate functional value
	tmp_maple.write("# functional value calculation\n")
	if infor.lambda_deriv:
		lambda_name = name + "_" + "lambda"
		line = lambda_name + ":= " + "diff (" + name + ", " + "lambda" + "):" 
		tmp_maple.write(line)
		tmp_maple.write("\n")
		tmp_maple.write("Fortran (%s):\n"%lambda_name)
	else:
		tmp_maple.write("Fortran (%s):\n"%name)

	# write the corresponding name	
	if infor.in_pointwise:
		tmp_maple.write("Text(\"! ###### F ###### \"):\n\n\n")
	else:
		tmp_maple.write("Text(\"! ###### F(i) ###### \"):\n\n\n")


	# get the first order derivatives
	tmp_maple.write("Text(\"IF (IDERIV .GE. 1) THEN\");\n")
	deriv_gen(1, routine_type, name, tmp_maple)
	tmp_maple.write ("Text(\"END IF\");\n")

	# get the second order derivatives
	if infor.functional_order >= 2:
		tmp_maple.write("Text(\"IF (IDERIV .GE. 2) THEN\");\n")
		deriv_gen(2, routine_type, name, tmp_maple)
		tmp_maple.write ("Text(\"END IF\");\n")

	# get the second order derivatives
	if infor.functional_order >= 3:
		tmp_maple.write("Text(\"IF (IDERIV .GE. 3) THEN\");\n")
		deriv_gen(3, routine_type, name, tmp_maple)
		tmp_maple.write ("Text(\"END IF\");\n")

	# for the close shell case, we have additional treatment
	if routine_type == "close":
		generate_rest_func_deriv_for_close(1, tmp_maple)
		if infor.functional_order >= 2:
			generate_rest_func_deriv_for_close(2, tmp_maple)
		if infor.functional_order >= 3:
			generate_rest_func_deriv_for_close(3, tmp_maple)

	# the ending of the subroutine 
	if ftype == "open_RB" or ftype == "close":
		if ftype == "open_RB":
			tmp_maple.write("Text(\"END IF\");\n")
		if not infor.in_pointwise:
			tmp_maple.write("Text(\"END DO\");\n")
		tmp_maple.write("Text(\"RETURN\");\n")
		tmp_maple.write("Text(\"END\");\n")
		tmp_maple.write("Text(\"             \");\n")
		tmp_maple.write("Text(\"             \");\n")
		tmp_maple.write("Text(\"             \");\n")
		tmp_maple.write("Text(\"             \");\n")
	tmp_maple.close()


#####################################################################
#
#  generate the additional functional derivatives, these functional
#  derivatives are not calcualted through maple. See the working
#  function in var.py
#
#####################################################################
def generate_rest_func_deriv_for_close(order, tmp_maple):

	# get the data
	if order == 1:
		DF_close = generation.get_D1F("close")
		DF_open  = generation.get_D1F("open")
		DF_type  = "D1F"
		ordern   = "1"
	elif order == 2:
		DF_close = generation.get_D2F("close")
		DF_open  = generation.get_D2F("open")
		DF_type  = "D2F"
		ordern   = "2"
	elif order == 3:
		DF_close = generation.get_D3F("close")
		DF_open  = generation.get_D3F("open")
		DF_type  = "D3F"
		ordern   = "3"
	else:
		print "order is not correct in the generate_rest_func_deriv_for_close!\n"
		sys.exit()
	DF_rest = var.get_rest_func_deriv_in_close(DF_close, DF_open)

	# do the real work by calling the work function
	if len(DF_rest) == 0:
		return
	else:
		tmp_maple.write("Text(\"IF (IDERIV .GE. %s) THEN\");\n"%ordern)
		for eachKey in DF_rest.keys():
			new_name = "ID" + "_" + eachKey + "_POS"
			old_name = "ID" + "_" + DF_rest[eachKey] + "_POS"

			# now we calculate position variable
			# and print it out
			# this is only for derivatives
			if order > 0:
				begin = new_name.find("ID")
				end   = new_name.find("_POS")
				v     = new_name[begin:end]
				varPos= v + "_POS"
				varay = "D1VARS"
				if order == 2:
					varay = "D2VARS"
				elif order == 3:
					varay = "D3VARS"
				new_line = varPos + "=" + varay + "(" + v +")" 
				tmp_maple.write("Text(\"%s\"):\n"%new_line)

			# now it's the functional assignment
			if infor.in_pointwise:
				new_functional = DF_type + "( " + new_name +")"
				old_functional = DF_type + "( " + old_name +")"
			else:
				new_functional = DF_type + "(i, " + new_name +")"
				old_functional = DF_type + "(i, " + old_name +")"
			line = new_functional + "= " + old_functional
			tmp_maple.write("Text(\"%s\"):\n"%line)
		tmp_maple.write ("Text(\"END IF\");\n")


#####################################################################
#
# judge the functional name in the original maple file according
# to the routine type
#
#####################################################################
def get_functional_name(ftype):

	if ftype == "close":
		return "functional_RA_RB"
	elif ftype == "open_RA_RB":
		return "functional_RA_RB"
	elif ftype == "open_RA":
		return "functional_RA"
	elif ftype == "open_RB":
		return "functional_RB"
	else:
		print "incorrect routine type provided in the get_functional_name!\n"
		sys.exit()



#####################################################################
#
# maple section to produce the functional derivatives
#
#####################################################################
def deriv_gen(order, rtype, fname, outfile):

	# get the variable list, here we only pay attention to the close shell case
	if rtype == "close":
		if order == 1:
			DF_open  = generation.get_D1F("open")
			DF_close = generation.get_D1F("close")
		elif order == 2:
			DF_open  = generation.get_D2F("open")
			DF_close = generation.get_D2F("close")
		elif order == 3:
			DF_open  = generation.get_D3F("open")
			DF_close = generation.get_D3F("close")
		else:
			print "error in the deriv_gen! Order is not correct!\n"
			print order
			sys.exit()
	else:
		if order == 1:
			DF_open  = generation.get_D1F("open")
		elif order == 2:
			DF_open  = generation.get_D2F("open")
		elif order == 3:
			DF_open  = generation.get_D3F("open")
		else:
			print "error in the deriv_gen! Order is not correct!\n"
			print order
			sys.exit()


	# the real working loop begins.....
	for fvar in DF_open:

		# let's examine something
   	# if only alpha exists, we do not evaluate derivatives containing 
   	# beta part; similarly, if only beta part exists, we do not evaluate
   	# derivatives containg alpha part
		if fname == "functional_RA":
			if re.search("B", fvar) is not None:
				continue
		elif fname == "functional_RB":
			if re.search("A", fvar) is not None:
				continue


		# get the functional name, new_functional is the derivatives from the
		# new functional
		new_functional = "ID"+"_"+fvar+"_"+"POS"
		if order > 1:
			string = fvar.split("_") 
			if order == 2:
				var = string[1]
				old_functional = "ID"+"_"+string[0]+"_"+"POS"		
			elif order == 3:
				var = string[2]
				old_functional = "ID"+"_"+string[0]+"_"+string[1]+"_"+"POS"
		elif order == 1:
			var = fvar
			if infor.lambda_deriv:
				old_functional = fname + "_" + "lambda"
			else:
				old_functional = fname


		# name_DF specified the "DnF" is the name
		if order == 1:
			name_DF = "D1F"
		elif order == 2:
			name_DF = "D2F"
		elif order == 3:
			name_DF = "D3F"
	
		# get the name
		if infor.in_pointwise:
			name_part2 = "( "
		else:
			name_part2 = "(i, "
		name_part3 = ")"
		name = name_DF + name_part2 + new_functional + name_part3

		# working codes
		outfile.write("# generate the functional values\n")
		line = new_functional + ":= " + "diff (" + old_functional + ", " + var + "):" 
		outfile.write(line)
		outfile.write("\n")

		# form the fortran codes
		if rtype == "close":
			if fvar in DF_close:
				outfile.write("Fortran(%s):\n"%new_functional)
				outfile.write("Text(\"! ###### %s ######  \"):\n"%name)
		else:
			outfile.write("Fortran(%s):\n"%new_functional)
			outfile.write("Text(\"! ###### %s ######  \"):\n"%name)



