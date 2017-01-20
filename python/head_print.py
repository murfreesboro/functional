"""
This module is used print out the fixing part of the Fortran codes, such
as title, comments etc.
"""

import re
import sys
import os
import infor

__author__  = "Fenglai Liu"
__date__    = "Oct, 2010"



#########################################################
#
#  print out reference section etc.
# 
#########################################################
def print_ref(input_file, output_file):

	output = open(output_file, 'a')

	# detect that whether we have reference in input file
	infile = open(input_file, 'r')
	has_ref = 0
	while True:
		line = infile.readline()
		if not line: break 
		if re.search(r"(?i)ref\b", line) is not None and re.search(r"^#", line) is not None:
			has_ref = 1
			break

	if has_ref:
		output.write("!--------------------------------------------------------------------------\n")
		output.write("!                 Reference section for the functional:\n")
		output.write("!--------------------------------------------------------------------------\n")

		while True:
			line = infile.readline()
			if re.search(r"(?i)end\b", line) is not None and re.search(r"^#", line) is not None: break
			line2 = line.replace("#", "!")
			output.write(line2)
		output.write("!\n!\n\n\n\n")

	output.close()
	infile.close()


#########################################################
#
#  print out main function
# 
#########################################################
def print_main_function(output_file):


	outfile = open(output_file, 'a')

	# set the project name
	project_name = infor.project
	if infor.functional_order == 2:
		project_name = project_name + "d2"
	elif infor.functional_order == 3:
		project_name = project_name + "d3"

	# the functional name and parameters
	title = "      subroutine functional_" + project_name
	outfile.write(title)
	outfile.write("\n")
	if infor.functional_order == 1:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
			outfile.write("     & F,D1F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
			outfile.write("     &LapA,LapB,F,D1F)\n")
	elif infor.functional_order == 2:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F,D2F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F,D2F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
			outfile.write("     & F,D1F,D2F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
			outfile.write("     &LapA,LapB,F,D1F,D2F)\n")
	else:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F,D2F,D3F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F,D2F,D3F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
			outfile.write("     & F,D1F,D2F,D3F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
			outfile.write("     &LapA,LapB,F,D1F,D2F,D3F)\n")
	outfile.write("      IMPLICIT NONE\n")


	# comment section
	outfile.write("!--------------------------------------------------------\n")
	outfile.write("! This routine is used to calculate the functional and\n")
	outfile.write("! functional derivatives values up to the third order.\n")
	outfile.write("! It's only the interface function, the real calculation\n")
	outfile.write("! work will be done by seperate routine.\n")
	outfile.write("! INPUT :\n")
	outfile.write("! INFOR : information related to the variables\n")
	outfile.write("! NG    : the number of grid points\n")
	outfile.write("! NDEN  : the number of densities\n")
	outfile.write("! TOL   : tolerance value for error\n")
	outfile.write("! rhoA  : the alpha electron density\n")
	outfile.write("! rhoB  : the beta  electron density\n")
	outfile.write("! DRhoA : the alpha rho\'\n")
	outfile.write("! DRhoB : the beta  rho\'\n")
	outfile.write("! TauA  : the alpha kinetic energy density\n")
	outfile.write("! TauB  : the beta  kinetic energy density\n")
	outfile.write("! LapA  : the alpha laplacian density\n")
	outfile.write("! LapB  : the beta  laplacian density\n")
	outfile.write("! OUTPUT:\n")
	outfile.write("! F     : functional values\n")
	outfile.write("! D1F   : the first  order functional derivatives\n")
	outfile.write("! D2F   : the second order functional derivatives\n")
	outfile.write("! D3F   : the third  order functional derivatives\n")
	outfile.write("!--------------------------------------------------------\n")

	# data specification
	outfile.write("      INTEGER INFOR(*)\n")
	outfile.write("      INTEGER NG,NDEN\n")
	outfile.write("      DOUBLE PRECISION TOL\n")
	outfile.write("      DOUBLE PRECISION rhoA(NG),rhoB(NG)\n")
	if infor.get_type() >= 2:
		outfile.write("      DOUBLE PRECISION DRhoA(NG,3),DRhoB(NG,3)\n")
	if infor.get_type() >= 3:
		outfile.write("      DOUBLE PRECISION TauA(NG),TauB(NG)\n")
	if infor.get_type() == 4:
		outfile.write("      DOUBLE PRECISION LapA(NG),LapB(NG)\n")

	if infor.functional_order == 1:
		outfile.write("      DOUBLE PRECISION F(NG), D1F(NG,*)\n")
	elif infor.functional_order == 2:
		outfile.write("      DOUBLE PRECISION F(NG), D1F(NG,*), D2F(NG,*)\n")
	else:
		outfile.write("      DOUBLE PRECISION F(NG), D1F(NG,*), D2F(NG,*),D3F(NG,*)\n")
	outfile.write("      \n")

	# real calculation begins 
	outfile.write("      IF (NDEN .EQ. 1) THEN \n")
	title = "      CALL functional_" + project_name + "_close"
	outfile.write(title)
	outfile.write("\n")
	if infor.functional_order == 1:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,F,D1F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,LapA,F,D1F)\n")
	elif infor.functional_order == 2:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F,D2F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F,D2F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,F,D1F,D2F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,LapA,F,D1F,D2F)\n")
	else:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F,D2F,D3F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F,D2F,D3F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,F,D1F,D2F,D3F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,LapA,F,D1F,D2F,D3F)\n")
	outfile.write("      ELSE  \n")
	title = "      CALL functional_" + project_name + "_open"
	outfile.write(title)
	outfile.write("\n")
	if infor.functional_order == 1:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,F,D1F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,LapA,LapB,\n")
			outfile.write("     & F,D1F)\n")
	elif infor.functional_order == 2:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F,D2F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F,D2F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,F,D1F,D2F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,LapA,LapB,\n")
			outfile.write("     & F,D1F,D2F)\n")
	else:
		if infor.get_type() == 1:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F,D2F,D3F)\n")
		elif infor.get_type() == 2:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F,D2F,D3F)\n")
		elif infor.get_type() == 3:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,F,D1F,D2F,D3F)\n")
		else:
			outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,LapA,LapB,\n")
			outfile.write("     & F,D1F,D2F,D3F)\n")

	outfile.write("      END IF  \n")
	outfile.write("      RETURN  \n")
	outfile.write("      END     \n\n\n\n")

	outfile.close()



#########################################################
#
#  print out the head for work function
# 
#########################################################
def print_work_function(output_file, routine_type):

	# firstly, check the routine type
	# we only print the function title for close shell
	# and RA_RB case
	if routine_type == "open_RA" or routine_type == "open_RB":
		return
	elif routine_type == "open_RA_RB":
		shell_type = "open"
	else:
		shell_type = "close"

	outfile = open(output_file, 'a')

	# set the project name
	project_name = infor.project
	if infor.functional_order == 2:
		project_name = project_name + "d2"
	elif infor.functional_order == 3:
		project_name = project_name + "d3"

	# the functional name and parameters
	if shell_type == "close": 
		title = "      subroutine functional_" + project_name + "_close"
		outfile.write(title)
		outfile.write("\n")
		if infor.functional_order == 1:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,F,D1F)\n")
			else:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,LapA,F,D1F)\n")
		elif infor.functional_order == 2:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F,D2F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F,D2F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,F,D1F,D2F)\n")
			else:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,LapA,F,D1F,D2F)\n")
		else:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,F,D1F,D2F,D3F)\n")
			else:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,DRhoA,TauA,LapA,F,D1F,D2F,D3F)\n")
		outfile.write("      IMPLICIT DOUBLE PRECISION (A-H,O-Z)\n")
		outfile.write("!--------------------------------------------------------\n")
		outfile.write("! This routine is used to calculate the functional and\n")
		outfile.write("! functional derivatives values up to the third order.\n")
		outfile.write("! We note that this routine is used in the close shell case,\n")
		outfile.write("! that means, we have RA=RB, GAA=GAB=GBB,TA=TB and LA=LB\n")
		outfile.write("! for each grid point.\n") 
		outfile.write("!--------------------------------------------------------\n")
		if infor.functional_order == 1:
			outfile.write("#include \"fderiv1.inc\"\n")
		elif infor.functional_order == 2:
			outfile.write("#include \"fderiv1.inc\"\n")
			outfile.write("#include \"fderiv2.inc\"\n")
		else:
			outfile.write("#include \"fderiv1.inc\"\n")
			outfile.write("#include \"fderiv2.inc\"\n")
			outfile.write("#include \"fderiv3.inc\"\n")
	elif shell_type == "open":
		title = "      subroutine functional_" + project_name + "_open"
		outfile.write(title)
		outfile.write("\n")
		if infor.functional_order == 1:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
				outfile.write("     & F,D1F)\n")
			else:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
				outfile.write("     & LapA,LapB,F,D1F)\n")
		elif infor.functional_order == 2:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F,D2F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F,D2F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
				outfile.write("     & F,D1F,D2F)\n")
			else:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
				outfile.write("     & LapA,LapB,F,D1F,D2F)\n")
		else:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
				outfile.write("     & F,D1F,D2F,D3F)\n")
			else:
				outfile.write("     &(INFOR,NG,NDEN,TOL,rhoA,rhoB,DRhoA,DRhoB,TauA,TauB,\n")
				outfile.write("     & LapA,LapB,F,D1F,D2F,D3F)\n")
		outfile.write("      IMPLICIT DOUBLE PRECISION (A-H,O-Z)\n")
		outfile.write("!--------------------------------------------------------\n")
		outfile.write("! This routine is used to calculate the functional and\n")
		outfile.write("! functional derivatives values up to the third order.\n")
		outfile.write("! We note that this routine is used in the open shell case.\n")
		outfile.write("!--------------------------------------------------------\n")
		if infor.functional_order == 1:
			outfile.write("#include \"fderiv1.inc\"\n")
		elif infor.functional_order == 2:
			outfile.write("#include \"fderiv1.inc\"\n")
			outfile.write("#include \"fderiv2.inc\"\n")
		else:
			outfile.write("#include \"fderiv1.inc\"\n")
			outfile.write("#include \"fderiv2.inc\"\n")
			outfile.write("#include \"fderiv3.inc\"\n")
	else:
		print "Illegal type of shell type in print_work_function!\n"
		sys.exit()

	# data specification
	outfile.write("      INTEGER INFOR(*)\n")
	outfile.write("      INTEGER IDERIV, NG\n")
	outfile.write("      DOUBLE PRECISION rhoA(NG)\n")
	outfile.write("      DOUBLE PRECISION RA, RB ! rho vaule at each point\n")
	if infor.get_type() >= 2:
		outfile.write("      DOUBLE PRECISION DRhoA(NG,3)\n")
		outfile.write("      DOUBLE PRECISION GAA, GAB, GBB ! the gamma value at each point\n")
	if infor.get_type() >= 3:
		outfile.write("      DOUBLE PRECISION TauA(NG)\n")
		outfile.write("      DOUBLE PRECISION TA, TB     ! the tau value at each point\n")
	if infor.get_type() == 4:
		outfile.write("      DOUBLE PRECISION LapA(NG)\n")
		outfile.write("      DOUBLE PRECISION LA, LB     ! the laplacian value at each point\n")
	if shell_type == "open":
		outfile.write("      DOUBLE PRECISION rhoB(NG)\n")
		if infor.get_type() == 2:
			outfile.write("      DOUBLE PRECISION DRhoB(NG,3)\n")
		elif infor.get_type() == 3:
			outfile.write("      DOUBLE PRECISION DRhoB(NG,3),TauB(NG)\n")
		else:
			outfile.write("      DOUBLE PRECISION DRhoB(NG,3),TauB(NG),LapB(NG)\n")
	if infor.functional_order == 1:
		outfile.write("      DOUBLE PRECISION F(NG),D1F(NG,*)\n")
	elif infor.functional_order == 2:
		outfile.write("      DOUBLE PRECISION F(NG),D1F(NG,*),D2F(NG,*)\n")
	else:
		outfile.write("      DOUBLE PRECISION F(NG),D1F(NG,*),D2F(NG,*),D3F(NG,*)\n")
	outfile.write("      DOUBLE PRECISION TOL ! tolerance\n")
	outfile.write("      DOUBLE PRECISION PI\n")
	if infor.use_lambda or infor.lambda_deriv:
		outfile.write("      DOUBLE PRECISION lambda\n")

	# finally it's the variable position information
	if infor.functional_order == 1:
		outfile.write("      INTEGER  D1VARS(N_FUNC_DERIV_1)\n")
	elif infor.functional_order == 2:
		outfile.write("      INTEGER  D1VARS(N_FUNC_DERIV_1)\n")
		outfile.write("      INTEGER  D2VARS(N_FUNC_DERIV_2)\n")
	else:
		outfile.write("      INTEGER  D1VARS(N_FUNC_DERIV_1)\n")
		outfile.write("      INTEGER  D2VARS(N_FUNC_DERIV_2)\n")
		outfile.write("      INTEGER  D3VARS(N_FUNC_DERIV_3)\n")
	outfile.write("      \n")

	# now initilize the variable array
	outfile.write("      ! firstly initilize variable position information\n")
	if infor.functional_order >= 1:
		outfile.write("      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)\n")
	if infor.functional_order >= 2:
		outfile.write("      CALL INIT_FUNC_DERIV_2(INFOR,D2VARS)\n")
	if infor.functional_order >= 3:
		outfile.write("      CALL INIT_FUNC_DERIV_3(INFOR,D3VARS)\n")
	outfile.write("      \n")

	# real codes section
	outfile.write("      ! deal with the constants\n")
	outfile.write("      PI    = 4.D0*DATAN(1.D0)\n")
	if infor.functional_order == 1:
		outfile.write("      IDERIV    = 1\n")
	elif infor.functional_order == 2:
		outfile.write("      IDERIV    = 2\n")
	elif infor.functional_order == 3:
		outfile.write("      IDERIV    = 3\n")
	if infor.use_lambda or infor.lambda_deriv:
		lambda_value = infor.lambda_value
		outfile.write("      lambda    = %s\n"%lambda_value)
	outfile.write("\n")

	# real calculation begins
	outfile.write("      ! the real calculation begins\n")
	outfile.write("      DO i = 1,NG\n")
	if shell_type == "close":
		outfile.write("      RA = rhoA(i)\n")
		outfile.write("      RB = RA \n")
		outfile.write("      IF (RA.LE.TOL) CYCLE\n")
		if infor.get_type() >= 2:
			outfile.write("      GAA = DRhoA(i,1)*DRhoA(i,1) + DRhoA(i,2)*DRhoA(i,2)\n")
			outfile.write("     &+ DRhoA(i,3)*DRhoA(i,3)\n")
			outfile.write("      GBB = GAA\n")
			outfile.write("      GAB = GAA\n")
		if infor.get_type() >= 3:
			outfile.write("      TA = TauA(i)\n")
			outfile.write("      TB = TA\n")
		if infor.get_type() >= 4:
			outfile.write("      LA = LapA(i)\n")
			outfile.write("      LB = LA\n")
	elif shell_type == "open":
		outfile.write("         RA = rhoA(i)\n")
		outfile.write("         RB = rhoB(i)\n")
		if infor.get_type() >= 2:
			outfile.write("         GAA = DRhoA(i,1)*DRhoA(i,1) + DRhoA(i,2)*DRhoA(i,2)\n")
			outfile.write("     &+ DRhoA(i,3)*DRhoA(i,3)\n")
			outfile.write("         GBB = DRhoB(i,1)*DRhoB(i,1) + DRhoB(i,2)*DRhoB(i,2)\n")
			outfile.write("     &+ DRhoB(i,3)*DRhoB(i,3)\n")
			outfile.write("         GAB = DRhoA(i,1)*DRhoB(i,1) + DRhoA(i,2)*DRhoB(i,2)\n")
			outfile.write("     &+ DRhoA(i,3)*DRhoB(i,3)\n")
		if infor.get_type() >= 3:
			outfile.write("         TA = TauA(i)\n")
			outfile.write("         TB = TauB(i)\n")
		if infor.get_type() >= 4:
			outfile.write("         LA = LapA(i)\n")
			outfile.write("         LB = LapB(i)\n")
		outfile.write("      IF (RB.LE.TOL) THEN\n")
		outfile.write("         RB = 0.0D0 \n")
		if infor.get_type() >= 2:
			outfile.write("         GBB = 0.0D0\n")
			outfile.write("         GAB = 0.0D0\n")
		if infor.get_type() >= 3:
			outfile.write("         TB = 0.0D0\n")
		if infor.get_type() >= 4:
			outfile.write("         LB = 0.0D0\n")
		outfile.write("      END IF \n")
		outfile.write("      IF (RA.LE.TOL) THEN\n")
		outfile.write("         RA = 0.0D0 \n")
		if infor.get_type() >= 2:
			outfile.write("         GAA = 0.0D0\n")
			outfile.write("         GAB = 0.0D0\n")
		if infor.get_type() >= 3:
			outfile.write("         TA = 0.0D0\n")
		if infor.get_type() >= 4:
			outfile.write("         LA = 0.0D0\n")
		outfile.write("      END IF \n")
		outfile.write("      IF (RA.LE.TOL .AND. RB.LE.TOL) THEN\n")
		outfile.write("         CYCLE \n")
		outfile.write("      END IF \n")
		outfile.write("             \n")
	else:
		print "Illegal type of shell type in print_work_function!\n"
		sys.exit()


	outfile.close()


#########################################################
#
#  print out the head for work function in pointwise
#  calculation
# 
#########################################################
def print_work_function_pointwise(output_file, routine_type):

	# firstly, check the routine type
	# we only print the function title for close shell
	# and RA_RB case
	if routine_type == "open_RA" or routine_type == "open_RB":
		return
	elif routine_type == "open_RA_RB":
		shell_type = "open"
	else:
		shell_type = "close"

	outfile = open(output_file, 'a')

	# set the project name
	project_name = infor.project
	if infor.functional_order == 2:
		project_name = project_name + "d2"
	elif infor.functional_order == 3:
		project_name = project_name + "d3"

	# the functional name and parameters
	if shell_type == "close": 
		title = "      subroutine functional_" + project_name + "_pw" + "_close"
		outfile.write(title)
		outfile.write("\n")
		if infor.functional_order == 1:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,F,D1F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,F,D1F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,F,D1F)\n")
			else:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,LapA,F,D1F)\n")
		elif infor.functional_order == 2:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,F,D1F,D2F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,F,D1F,D2F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,F,D1F,D2F)\n")
			else:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,LapA,F,D1F,D2F)\n")
		else:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,F,D1F,D2F,D3F)\n")
			else:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,DRhoAX,DRhoAY,DRhoAZ,TauA,LapA,F,D1F,D2F,D3F)\n")
		outfile.write("      IMPLICIT DOUBLE PRECISION (A-H,O-Z)\n")
		outfile.write("!--------------------------------------------------------\n")
		outfile.write("! This routine is used to calculate the functional and\n")
		outfile.write("! functional derivatives pointwisely.\n")
		outfile.write("! We note that this routine is used in the close shell case,\n")
		outfile.write("! that means, we have RA=RB, GAA=GAB=GBB,TA=TB and LA=LB.\n")
	elif shell_type == "open":
		title = "      subroutine functional_" + project_name + "_pw" + "_open"
		outfile.write(title)
		outfile.write("\n")
		if infor.functional_order == 1:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,F,D1F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,F,D1F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,F,D1F)\n")
			else:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,LapA,LapB,F,D1F)\n")
		elif infor.functional_order == 2:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,F,D1F,D2F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,F,D1F,D2F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,F,D1F,D2F)\n")
			else:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,LapA,LapB,F,D1F,D2F)\n")
		else:
			if infor.get_type() == 1:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 2:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,F,D1F,D2F,D3F)\n")
			elif infor.get_type() == 3:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,F,D1F,D2F,D3F)\n")
			else:
				outfile.write("     &(INFOR,NDEN,TOL,rhoA,rhoB,DRhoAX,DRhoAY,DRhoAZ,\n")
				outfile.write("     & DRhoBX,DRhoBY,DRhoBZ,TauA,TauB,LapA,LapB,F,D1F,D2F,D3F)\n")
		outfile.write("      IMPLICIT DOUBLE PRECISION (A-H,O-Z)\n")
		outfile.write("!--------------------------------------------------------\n")
		outfile.write("! This routine is used to calculate the functional and\n")
		outfile.write("! functional derivatives pointwisely.\n")
		outfile.write("! We note that this routine is used in the open shell case.\n")
	else:
		print "Illegal type of shell type in print_work_function!\n"
		sys.exit()

	# print out the include files
	if infor.functional_order == 1:
		outfile.write("#include \"fderiv1.inc\"\n")
	elif infor.functional_order == 2:
		outfile.write("#include \"fderiv1.inc\"\n")
		outfile.write("#include \"fderiv2.inc\"\n")
	else:
		outfile.write("#include \"fderiv1.inc\"\n")
		outfile.write("#include \"fderiv2.inc\"\n")
		outfile.write("#include \"fderiv3.inc\"\n")

	# comment section
	outfile.write("! INTPUT :\n")
	outfile.write("! INFOR  : variable information\n")
	outfile.write("! NDEN   : number of densities\n")
	outfile.write("! rhoA   : the alpha electron density\n")
	outfile.write("! rhoB   : the beta  electron density\n")
	outfile.write("! DRhoAX : the gradient of rho alpha on x direction\n")
	outfile.write("! DRhoAY : the gradient of rho alpha on y direction\n")
	outfile.write("! DRhoAZ : the gradient of rho alpha on z direction\n")
	outfile.write("! DRhoBX : the gradient of rho beta  on x direction\n")
	outfile.write("! DRhoBY : the gradient of rho beta  on y direction\n")
	outfile.write("! DRhoBZ : the gradient of rho beta  on z direction\n")
	outfile.write("! TauA   : the alpha kinetic energy density\n")
	outfile.write("! TauB   : the beta  kinetic energy density\n")
	outfile.write("! LapA   : the alpha laplacian density\n")
	outfile.write("! LapB   : the beta  laplacian density\n")
	outfile.write("! OUTPUT :\n")
	outfile.write("! F      : functional values\n")
	outfile.write("! D1F    : the first  order functional derivatives\n")
	outfile.write("! D2F    : the second order functional derivatives\n")
	outfile.write("! D3F    : the third  order functional derivatives\n")
	outfile.write("!--------------------------------------------------------\n")


	# data specification
	outfile.write("      INTEGER INFOR(*)\n")
	outfile.write("      INTEGER IDERIV\n")
	outfile.write("      DOUBLE PRECISION rhoA\n")
	outfile.write("      DOUBLE PRECISION RA, RB \n")
	if infor.get_type() >= 2:
		outfile.write("      DOUBLE PRECISION DRhoAX, DRhoAY, DRhoAZ\n")
		outfile.write("      DOUBLE PRECISION GAA, GAB, GBB \n")
	if infor.get_type() >= 3:
		outfile.write("      DOUBLE PRECISION TauA\n")
		outfile.write("      DOUBLE PRECISION TA, TB\n")
	if infor.get_type() == 4:
		outfile.write("      DOUBLE PRECISION LapA\n")
		outfile.write("      DOUBLE PRECISION LA, LB\n")
	if shell_type == "open":
		outfile.write("      DOUBLE PRECISION rhoB\n")
		if infor.get_type() >= 2:
			outfile.write("      DOUBLE PRECISION DRhoBX, DRhoBY, DRhoBZ\n")
		if infor.get_type() >= 3:
			outfile.write("      DOUBLE PRECISION TauB\n")
		if infor.get_type() == 4:
			outfile.write("      DOUBLE PRECISION LapB\n")
	if infor.functional_order == 1:
		outfile.write("      DOUBLE PRECISION F,D1F(*)\n")
	elif infor.functional_order == 2:
		outfile.write("      DOUBLE PRECISION F,D1F(*),D2F(*)\n")
	else:
		outfile.write("      DOUBLE PRECISION F,D1F(*),D2F(*),D3F(*)\n")
	outfile.write("      DOUBLE PRECISION TOL ! tolerance\n")
	outfile.write("      DOUBLE PRECISION PI\n")
	if infor.use_lambda or infor.lambda_deriv:
		outfile.write("      DOUBLE PRECISION lambda\n")

	# finally it's the variable position information
	if infor.functional_order == 1:
		outfile.write("      INTEGER  D1VARS(N_FUNC_DERIV_1)\n")
	elif infor.functional_order == 2:
		outfile.write("      INTEGER  D1VARS(N_FUNC_DERIV_1)\n")
		outfile.write("      INTEGER  D2VARS(N_FUNC_DERIV_2)\n")
	else:
		outfile.write("      INTEGER  D1VARS(N_FUNC_DERIV_1)\n")
		outfile.write("      INTEGER  D2VARS(N_FUNC_DERIV_2)\n")
		outfile.write("      INTEGER  D3VARS(N_FUNC_DERIV_3)\n")
	outfile.write("      \n")

	# now initilize the variable array
	outfile.write("      ! firstly initilize variable position information\n")
	if infor.functional_order >= 1:
		outfile.write("      CALL INIT_FUNC_DERIV_1(INFOR,D1VARS)\n")
	if infor.functional_order >= 2:
		outfile.write("      CALL INIT_FUNC_DERIV_2(INFOR,D2VARS)\n")
	if infor.functional_order >= 3:
		outfile.write("      CALL INIT_FUNC_DERIV_3(INFOR,D3VARS)\n")
	outfile.write("      \n")

	# real codes section
	outfile.write("      ! firstly deal with the constants\n")
	outfile.write("      PI    = 4.D0*DATAN(1.D0)\n")
	if infor.functional_order == 1:
		outfile.write("      IDERIV    = 1\n")
	elif infor.functional_order == 2:
		outfile.write("      IDERIV    = 2\n")
	elif infor.functional_order == 3:
		outfile.write("      IDERIV    = 3\n")
	if infor.use_lambda or infor.lambda_deriv:
		lambda_value = infor.lambda_value
		outfile.write("      lambda    = %s\n"%lambda_value)
	outfile.write("\n")

	# real calculation begins
	outfile.write("      ! the real calculation begins\n")
	if shell_type == "close":
		outfile.write("      RA = rhoA\n")
		outfile.write("      RB = RA \n")
		outfile.write("      IF (RA.LE.TOL) RETURN\n")
		if infor.get_type() >= 2:
			outfile.write("      GAA = DRhoAX*DRhoAX + DRhoAY*DRhoAY\n")
			outfile.write("     &+ DRhoAZ*DRhoAZ\n")
			outfile.write("      GBB = GAA\n")
			outfile.write("      GAB = GAA\n")
		if infor.get_type() >= 3:
			outfile.write("      TA = TauA\n")
			outfile.write("      TB = TA\n")
		if infor.get_type() >= 4:
			outfile.write("      LA = LapA\n")
			outfile.write("      LB = LA\n")
	elif shell_type == "open":
		outfile.write("         RA = rhoA\n")
		outfile.write("         RB = rhoB\n")
		if infor.get_type() >= 2:
			outfile.write("         GAA = DRhoAX*DRhoAX + DRhoAY*DRhoAY\n")
			outfile.write("     &+ DRhoAZ*DRhoAZ\n")
			outfile.write("         GBB = DRhoBX*DRhoBX + DRhoBY*DRhoBY\n")
			outfile.write("     &+ DRhoBZ*DRhoBZ\n")
			outfile.write("         GAB = DRhoAX*DRhoBX + DRhoAY*DRhoBY\n")
			outfile.write("     &+ DRhoAZ*DRhoBZ\n")
		if infor.get_type() >= 3:
			outfile.write("         TA = TauA\n")
			outfile.write("         TB = TauB\n")
		if infor.get_type() >= 4:
			outfile.write("         LA = LapA\n")
			outfile.write("         LB = LapB\n")
		outfile.write("      IF (RB.LE.TOL) THEN\n")
		outfile.write("         RB = 0.0D0 \n")
		if infor.get_type() >= 2:
			outfile.write("         GBB = 0.0D0\n")
			outfile.write("         GAB = 0.0D0\n")
		if infor.get_type() >= 3:
			outfile.write("         TB = 0.0D0\n")
		if infor.get_type() >= 4:
			outfile.write("         LB = 0.0D0\n")
		outfile.write("      END IF \n")
		outfile.write("      IF (RA.LE.TOL) THEN\n")
		outfile.write("         RA = 0.0D0 \n")
		if infor.get_type() >= 2:
			outfile.write("         GAA = 0.0D0\n")
			outfile.write("         GAB = 0.0D0\n")
		if infor.get_type() >= 3:
			outfile.write("         TA = 0.0D0\n")
		if infor.get_type() >= 4:
			outfile.write("         LA = 0.0D0\n")
		outfile.write("      END IF \n")
		outfile.write("      IF (RA.LE.TOL .AND. RB.LE.TOL) THEN\n")
		outfile.write("         RETURN \n")
		outfile.write("      END IF \n")
		outfile.write("             \n")
	else:
		print "Illegal type of shell type in print_work_function!\n"
		sys.exit()


	outfile.close()
