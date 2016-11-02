cdef extern from "/home/wfarah/realtime_model/src/get_nsmd.c":
	float getmd(char * utc_start_str, char * ra, char * dec)
	float getns(char * utc_start_str, char * ra, char * dec)
