import re

def _read_cz_info(cus_file):
	file_ = open(cus_file, 'r')
	Lats =  []
	Longts = []
	Demands = []
	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)
	nb_cus = int(list_[0])
	while True:
		line_ = file_.readline()
		list_ = re.split(" |\t|\n", line_)
		if list_[0] == '':
			break
		else:
			Lats.append(float(list_[0]))
			Longts.append(float(list_[1]))
			Demands.append(float(list_[2]))
	return Lats, Longts, Demands


def _read_fac_info(fac_file):
	file_ = open(fac_file, 'r')
	Lats =  []
	Longts = []
	line_ = file_.readline()
	list_ = re.split(" |\t|\n", line_)
	nb_fac = int(list_[0])
	itr = 1
	while True:
		line_ = file_.readline()
		list_ = re.split(" |\t|\n", line_)
		if list_[0] == '':
			break
		else:
			Lats.append(float(list_[0]))
			Longts.append(float(list_[1]))
		itr += 1
		if itr > nb_fac:
			break
	return Lats, Longts

from math import sin, cos, sqrt, atan2, radians
def _calc_geo_distance(lat1, lon1, lat2, lon2):
	# approximate radius of earth in mile
	R = 3956.0
	lat1 = radians(lat1)
	lon1 = radians(lon1)
	lat2 = radians(lat2)
	lon2 = radians(lon2)
	dlon = lon2 - lon1
	dlat = lat2 - lat1
	a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
	c = 2 * atan2(sqrt(a), sqrt(1 - a))
	distance = R * c
	return distance

def _calc_dist_matx(czLats, czLongts, czDemands, facLats, facLongts):
	dist_matx = []
	nbCz = len(czLats)
	nbFac = len(facLats) 
	for c in range(nbCz):
		tmp = []
		for f in range(nbFac):
			gdist = _calc_geo_distance(czLats[c], czLongts[c], facLats[f], facLongts[f])
			tmp.append(gdist * czDemands[c])
		dist_matx.append(tmp)

	return dist_matx