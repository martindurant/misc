import string
from string import atof

"""
Float <-> sexagessimal string conversion
"""

def todec(inlist,sep = ":"):
	"turn string of sexagesimal strigns into decimal float"
	words=inlist.split(sep)
	deg=atof(words[0])
	min=atof(words[1])/60.0
	sec=atof(words[2])/3600.0
	if deg<0:
            return deg-min-sec
	return deg+min+sec

def tosex(inlist,sep=":"):
	"turn decimal into sexagesial string"
	deg=int(number)
	min=int((number-deg)*60)
	sec=int(number-deg-min/60.0)*60
	outlist=sep.join([deg,min,sec])
	return outlist

