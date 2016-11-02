from cnsmd cimport getmd,getns


def wrapper_getmd(char * utc_start_str, char * ra, char * dec):
	cdef float md
	md = getmd(utc_start_str, ra, dec)
	return md

def wrapper_getns(char * utc_start_str, char * ra, char * dec):
	cdef float ns
	ns = getmd(utc_start_str, ra, dec)
	return ns
