"""
This module is to provide a grounp of functions which are used in dealing with 
variables handling
"""

import sys
import infor

__author__  = "Fenglai Liu"
__version__ = "0.0"
__date__    = "Oct, 2010"



def is_var_inside (var, DF):
	"""
	This subroutine is used to check that whether we have 
	the same variable inside the list. 
	"""

	tmp_var_list = var.split("_")
	if len(tmp_var_list) == 2:
		tmp_var = tmp_var_list[1] + "_" + tmp_var_list[0]
		if tmp_var in DF:
			return 1
		else:
			return 0
	elif len(tmp_var_list) == 3:
		tmp1 = tmp_var_list[0] + "_" + tmp_var_list[2] + "_" + tmp_var_list[1]
		tmp2 = tmp_var_list[1] + "_" + tmp_var_list[0] + "_" + tmp_var_list[2]
		tmp3 = tmp_var_list[2] + "_" + tmp_var_list[0] + "_" + tmp_var_list[1]
		tmp4 = tmp_var_list[1] + "_" + tmp_var_list[2] + "_" + tmp_var_list[0]
		tmp5 = tmp_var_list[2] + "_" + tmp_var_list[1] + "_" + tmp_var_list[0]
		if tmp1 in DF or tmp2 in DF or tmp3 in DF or tmp4 in DF or tmp5 in DF:
			return 1
		else:
			return 0
	else:
		print "check the length of DF!\n"
		sys.exit()


def arrange_var_for_close (DF):
	"""
	This function is used to re-arrange the var list for the close shell case
	"""

	# first of all, reverse the whole array
	DF.reverse()
	DF_original = DF[:]

	# exchange the alpha and beta names, then detect whether we have it
	for i in DF_original:
		if not var_in_non_equi_list (i):
			j = exchange_alpha_beta(i)
			k = re_arrange_var(j)
			if k in DF and k != i:
				DF.remove(i)
	
	# finally, reverse it again
	DF.reverse()


def get_rest_func_deriv_in_close (DF_close, DF_open):
	"""
	This function is used to get the rest of the functional derivatives
	which is not calculated directly in maple. For example,
	GAA_GAA is symmetrical with GBB_GBB
	So we do not calculate GBB_GBB, instead we have 
	D2F(i,GBB_GBB) = D2F(i,GAA_GAA)
	Hence here we just get the list of functional derivatives for the 
	examples of GBB_GBB
	"""

	DF = {}
	# let's make up the dict
	keys = []
	for i in DF_open:
		if i in DF_close:
			continue
		else:
			keys.append(i)
	# actually here means the open shell and close shell list are equal to each other
	if len(keys) == 0:
		return DF
	else:
		DF = DF.fromkeys(keys)

	# whether it's the first order functional derivatives?
	# we can know it by testing the first element
	var_list = DF_close[0].split("_")
	if len(var_list) == 1:
		for key in DF.keys():
			i = key
			if i is not None:
				j = exchange_alpha_beta(i)
				k = j[0]
				if k not in DF_close:
					print "something wrong in matching DF_close!", i, k
					sys.exit()
				else:
					DF[i] = k
	else:
		for key in DF.keys():
			i = key
			if i is not None:
				j = exchange_alpha_beta(i)
				k = re_arrange_var(j)
				if k not in DF_close:
					print "something wrong in matching DF_close!", i, k 
					sys.exit()
				else:
					DF[i] = k

	return DF


def var_in_non_equi_list(i):
	"""
	test that for the given varible list whether we have var in the 
	non-quivalent list (see infor.py). If we have it, then we will
	pass the rest of the steps in the close shell procedure.
	"""
	var_list = i.split("_")
	for j in var_list:
		if j in infor.non_equi_var_list:
			return 1
	
	return 0


def exchange_alpha_beta(i):
	"""
	This function is used to change the one pure spin state into the 
	other, then we return a new list
	"""

	var_list = i.split("_")
	result   = [] 
	for j in var_list:
		if j == "RA": 
			result.append("RB")
		elif j == "RB": 
			result.append("RA")
		elif j == "GAA": 
			result.append("GBB")
		elif j == "GBB": 
			result.append("GAA")
		elif j == "TA":
			result.append("TB")
		elif j == "TB":
			result.append("TA")
		elif j == "LA":
			result.append("LB")
		elif j == "LB":
			result.append("LA")
		elif j == "GAB": 
			result.append("GAB")
		else:
			print "something wrong in the exchange_alpha_beta, no correct var got\n"
			sys.exit()

	return result


def re_arrange_var(i):
	"""
	This function is used to re-arrange the varible list so that 
	to adapt it for following rules:
	1  alpha is prior to beta;
	2  density > gradient density > tau > lap
	"""

	# first step, find the variable in the list according to the order,
	# and re-arrange it 
	result_list = []
	for l in range(len(i)):
		if "RA" in i and i.index("RA") >= 0:
			result_list.append("RA")
			i.remove("RA")
		elif "RB" in i and i.index("RB") >= 0:
			result_list.append("RB")
			i.remove("RB")
		elif "GAA" in i and i.index("GAA") >= 0:
			result_list.append("GAA")
			i.remove("GAA")
		elif "GAB" in i and i.index("GAB") >= 0:
			result_list.append("GAB")
			i.remove("GAB")
		elif "GBB" in i and i.index("GBB") >= 0:
			result_list.append("GBB")
			i.remove("GBB")
		elif "TA" in i and i.index("TA") >= 0:
			result_list.append("TA")
			i.remove("TA")
		elif "TB" in i and i.index("TB") >= 0:
			result_list.append("TB")
			i.remove("TB")
		elif "LA" in i and i.index("LA") >= 0:
			result_list.append("LA")
			i.remove("LA")
		elif "LB" in i and i.index("LB") >= 0:
			result_list.append("LB")
			i.remove("LB")
		else:
			print "something wrong in the re_arrange_var, no correct var got\n"
			sys.exit()


	# second step, organize it into a standard string
	if len(result_list) == 2:
		result = result_list[0] + "_" + result_list[1]
	elif len(result_list) == 3:
		result = result_list[0] + "_" + result_list[1] + "_" + result_list[2]
	else:
		print "something wrong in the re_arrange_var, the length should be 2 or 3\n"
		sys.exit()

#	print result
	return result

