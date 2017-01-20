"""
This module is used to build up the variable list for functional derivatives used
in the future maple file. The variable list is like GAA_GBB_TA in the 
D3F(i, POS_GAA_GBB_TA) etc. Here we also take care of the open shell and 
close shell difference so that to speed up the resulting codes.
"""

import sys
import var
import infor

__author__  = "Fenglai Liu"
__version__ = "0.0"
__date__    = "Oct, 2010"



D1F_close = []
D1F_open  = []
D2F_close = []
D2F_open  = []
D3F_close = []
D3F_open  = []


###############################################################################
#
#  setting the 1st order functional derivatives
#
###############################################################################

def set_D1F():


	global D1F_close
	global D1F_open

	if infor.functional_type == 1:
		D1F_open= ["RA", "RB"]
	elif infor.functional_type == 2:
		D1F_open= ["RA", "RB", "GAA", "GBB", "GAB"]
	elif infor.functional_type == 3:
		D1F_open= ["RA", "RB", "GAA", "GBB", "GAB", "TA", "TB"]
	elif infor.functional_type == 4:
		D1F_open= ["RA", "RB", "GAA", "GBB", "GAB", "TA", "TB", "LA", "LB"]
	else:
		print "functional type error in set_d1F of generation module!\n"
		sys.exit()


	# here we have to make a statement, that we believe all of the first
	# order functional derivatives should be symmetric between alpha and 
	# beta density(more information please check the additional doc). 
	# However, if this is not the case, please use "full var" in the maple 
	# input section so that we will do infor.use_full_var
	if infor.use_full_var:
		D1F_close = D1F_open
	else:
		if infor.functional_type == 1:
			D1F_close = ["RA"]
		elif infor.functional_type == 2:
			D1F_close = ["RA", "GAA", "GAB"]
		elif infor.functional_type == 3:
			D1F_close = ["RA", "GAA", "GAB", "TA"]
		else:
			D1F_close = ["RA", "GAA", "GAB", "TA", "LA"]


def get_D1F (ftype):

	if ftype == "close":
		return D1F_close
	elif ftype == "open":
		return D1F_open
	else:
		print "Illegal functional type required in get_D1F of generation module!\n"
		sys.exit()


###############################################################################
#
#  setting the 2ed order functional derivatives
#
###############################################################################

def set_D2F():

	rho  = ("RA","RB")
	grho = ("GAA","GBB","GAB")
	tau  = ("TA","TB")
	lap  = ("LA","LB")

	# generate both derivative variables for close and open shell cases
	i = 0
	while i < 2:

		if i == 0:
			state = "open"
		else:
			state = "close"

		if infor.functional_type >= 1:
			D2F_var_generate (rho,rho, state)
		if infor.functional_type >= 2:
			D2F_var_generate (rho,grho, state)
			D2F_var_generate (grho,grho, state)
		if infor.functional_type >= 3:
			D2F_var_generate (rho,tau, state)
			D2F_var_generate (grho,tau, state)
			D2F_var_generate (tau,tau, state)
		if infor.functional_type == 4:
			D2F_var_generate (rho,lap, state)
			D2F_var_generate (grho,lap, state)
			D2F_var_generate (tau,lap, state)
			D2F_var_generate (lap,lap, state)

		i = i + 1

	# debug codes
	#print "the final length for D2F_open is", len(D2F_open)
	#print "final length for D2F_open should be", infor.get_deriv_number(2)

def D2F_var_generate(var_list1, var_list2, state):

	global D2F_close
	global D2F_open

	# firstly generate the D2F_open, that is for the open shell case
	D2F = []
	for i in var_list1:
		for j in var_list2:
			varible = i+"_"+j
			if not var.is_var_inside(varible,D2F):
				D2F.append(varible)

	# if it's close state, then we have to delete the equivalent part
	# according to the non-equi-var-list, that is done is the arrange
	# var for close 
	if state == "open":
		D2F_open  += D2F
	else:
		var.arrange_var_for_close (D2F)
		D2F_close += D2F


def get_D2F (ftype):

	if ftype == "close":
		return D2F_close
	elif ftype == "open":
		return D2F_open
	else:
		print "Illegal functional type required in get_D2F of generation module!\n"
		sys.exit()



###############################################################################
#
#  setting the 3rd order functional derivatives
#
###############################################################################

def set_D3F():

	rho  = ("RA","RB")
	grho = ("GAA","GBB","GAB")
	tau  = ("TA","TB")
	lap  = ("LA","LB")

	# generate both derivative variables for close and open shell cases
	i = 0
	while i < 2:

		if i == 0:
			state = "open"
		else:
			state = "close"

		if infor.functional_type >= 1:
			D3F_var_generate (rho, rho, rho, state)
		if infor.functional_type >= 2:
			D3F_var_generate (rho,  rho,  grho, state)
			D3F_var_generate (rho,  grho, grho, state)
			D3F_var_generate (grho, grho, grho, state)
		if infor.functional_type >= 3:
			D3F_var_generate (rho,  rho,  tau, state)
			D3F_var_generate (rho,  grho, tau, state)
			D3F_var_generate (grho, grho, tau, state)
			D3F_var_generate (rho,  tau,  tau, state)
			D3F_var_generate (grho, tau,  tau, state)
			D3F_var_generate (tau,  tau,  tau, state)
		if infor.functional_type == 4:
			D3F_var_generate (rho,  rho,  lap, state)
			D3F_var_generate (rho,  grho, lap, state)
			D3F_var_generate (grho, grho, lap, state)
			D3F_var_generate (rho,  tau,  lap, state)
			D3F_var_generate (grho, tau,  lap, state)
			D3F_var_generate (tau,  tau,  lap, state)
			D3F_var_generate (rho,  lap,  lap, state)
			D3F_var_generate (grho, lap,  lap, state)
			D3F_var_generate (tau,  lap,  lap, state)
			D3F_var_generate (lap,  lap,  lap, state)

		i = i + 1

	
	# debug codes
	#print "the final length for D3F_open is", len(D3F_open)
	#print "final length for D3F_open should be", infor.get_deriv_number(3)


def D3F_var_generate(var_list1, var_list2, var_list3, state):

	global D3F_close
	global D3F_open

	# firstly generate the D2F_open, that is for the open shell case
	D3F = []
	for i in var_list1:
		for j in var_list2:
			for k in var_list3:
				varible = i+"_"+j+"_"+k
				if not var.is_var_inside(varible,D3F):
					D3F.append(varible)

	# if it's close state, then we have to delete the equivalent part
	# according to the non-equi-var-list, that is done is the arrange
	# var for close 
	if state == "open":
		D3F_open  += D3F
	else:
		var.arrange_var_for_close (D3F)
		D3F_close += D3F


def get_D3F (ftype):

	if ftype == "close":
		return D3F_close
	elif ftype == "open":
		return D3F_open
	else:
		print "Illegal functional type required in get_D3F of generation module!\n"
		sys.exit()


