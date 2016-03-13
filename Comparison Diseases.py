from __future__ import division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t

CATERGORIS = ('No', 'Py2.7', 'Py3')
ALPHA = 0.05

def _filter_data(frame, col, names):
	"""Return multiple dataframes based on the column value"""
	return map(lambda x: frame.loc[frame[col]==x], names)

def get_data():
	fever = pd.read_csv('fever.csv')
	wbc = pd.read_csv('WBC.csv')
	feverdata = _filter_data(fever, 'Diagnosis', CATERGORIS)
	wbcdata = _filter_data(wbc, 'Diagnosis', CATERGORIS)
	return feverdata, wbcdata

def _numpy_data(frames, name):
	"""Return multiple numpy array based on the number of frames"""
	return map(lambda frame: np.array(frame[name]), frames)

def t_test(x1, x2, s1, s2, n1, n2):
	return abs(x1 - x2) / (s1 ** 2 / n1 + s2 ** 2 / n2) ** 0.5

def _stats(arr1, arr2):
	"""Return the stats which t_test function need"""
	return arr1.mean(), arr2.mean(), arr1.std(), arr2.std(), len(arr1), len(arr2)

def cohen_d(x1, x2, s1, s2, n1, n2):
	return abs(x1 - x2) / (((n1 - 1) * (s1 ** 2) + (n2 - 1) * (s2 ** 2)) / (n1 + n2 - 2)) ** 0.5

def comparison_test(arr1, arr2, name1, name2, tail=1, bonf=1):
	stats = _stats(arr1, arr2)
	p_value = 1 - t.cdf(t_test(*stats), min(len(arr1), len(arr2)) - 1)
	if tail == 2: #In case for two tail test
		p_value *= 2
	alpha = ALPHA / bonf
	print "The p-value between {0} and {1} is {2:.5f} and it {3} the null hypothesis".format(
			name1, name2, p_value,
			"reject" if p_value < alpha else "does not reject")
	print "Actual Alpha is", format(alpha, '.4f')
	if p_value < alpha:
		print "The effect size is {0:.5f}".format(cohen_d(*stats))

def confi_interval(x, s, n, alpha):
	margin = t.ppf(1 - (alpha / 2), n - 1) * (s / n ** 0.5)
	print "Lower bound is {0:.2f} and higher bound is {1:.2f}".format(x - margin, x + margin)
	return x - margin, x, x, x, x + margin #the return structure is for visualization purpose

def _stats_confi(arr):
	"""Return parameters that can be used for function confi_interval"""
	return arr.mean(), arr.std(), len(arr), ALPHA

def subnormal_rate(arr, groupname):
	print "Below are WBC Group who has", groupname
	print "The overall subnormal rate(<4500) is ",\
			format(len(arr[arr < 4500]) / len(arr), '.5f')
	print "The abnormally low(<1700) is ",\
			format(len(arr[arr < 1700]) / len(arr), '.5f')
	print "The high risk of secondary infection(<500) rate is ",\
			format(len(arr[arr < 500]) / len(arr), '.5f')
	print 

def conditional(no, py2, py3, fever):
	num_py2 = len(py2[py2 > fever])
	num_py3 = len(py3[py3 > fever])
	num_no = len(no[no > fever])
	return (num_py2 + num_py3) / (num_py2 + num_py3 + num_no)

def analysis(data):
	f, w = data
	fno, fpy2, fpy3 = _numpy_data(f, 'High Temperature')
	wno, wpy2, wpy3 = _numpy_data(w, 'WBC')
	#run t-test for the difference of mean. Note that below are sparated tests
	print "*****************\nBelow are Body Temperature under Fever comparion test\n*****************\n"
	comparison_test(fno, fpy2, 'No disease', 'Py2', bonf=3)
	comparison_test(fno, fpy3, 'No disease', 'Py3', bonf=3)
	comparison_test(fpy2, fpy3, 'Py2', 'Py3', 2, 3)
	print "\n*****************\nBelow are WBC comparion test\n*****************\n"
	comparison_test(wno, wpy2, 'No disease', 'Py2', bonf=3)
	comparison_test(wno, wpy3, 'No disease', 'Py3', bonf=3)
	comparison_test(wpy2, wpy3, 'Py2', 'Py3', 2, 3)
	print "\n*****************\nBelow are WBC abnormal rates\n*****************\n"
	subnormal_rate(wno, "No disease")
	subnormal_rate(wpy2, "Py2 disease")
	subnormal_rate(wpy3, "Py3 disease")

	print "The probabiliy of having either diseases given fever > 102 F is",\
			format(conditional(fno, fpy2, fpy3, 102), '.5f'), '\n'
	
	fig, (ax1, ax2) = plt.subplots(2)
	graphdataf = [fno, fpy2, fpy3,\
				'The Confident Intervals For Diseases VS. Fever',\
				'The Body Temperature (Fahrenheit)']
	graphdataw = [wno, wpy2, wpy3,\
				'The Confident Intervals For Diseases VS. White Blood Cells',\
				'The Number of White Blood Cells']
	#Calculate the confident interval data
	for ax, data in {ax1: graphdataf, ax2: graphdataw}.iteritems():
		data1 = confi_interval(*_stats_confi(data[0]))
		data2 = confi_interval(*_stats_confi(data[1]))
		data3 = confi_interval(*_stats_confi(data[2]))	
		myax = ax.boxplot([data1, data2, data3])
		ax.set_title(data[3])
		ax.set_xlabel('The name of the diseases')
		ax.set_ylabel(data[4])
		ax.set_xticklabels(CATERGORIS)
		plt.setp(myax['whiskers'], color='k', linestyle='-')
	plt.subplots_adjust(hspace=0.3, top=0.95, bottom=0.07)
	plt.show()

analysis(get_data())
