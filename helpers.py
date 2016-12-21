import os
import sys
import atexit
from signal import SIGTERM,SIGKILL
import numpy as np
import time
import logging
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
sys.path.append("/home/observer/Python/dev/sigpyproc/lib/python")
from sigpyproc.Readers import FilReader
import scipy
import time

class Features:
	def __init__(self,beam,sample,sn,dm,box,F1,F2,F3,event,left,right,
			event_div,mean_off,std_off,sig_0,sig_1,sig_2,ks_d,ks_p,
			sw_w,sw_p,Mod_ind,Mod_indT,time):
		self.beam = beam
		self.sample = sample
		self.sn = sn
		self.dm = dm
		self.box = box
		self.F1 = F1
		self.F2 = F2
		self.F3 = F3
		self.event = event
		self.left = left
		self.right = right
		self.event_div = event_div
		self.mean_off = mean_off
		self.std_off = std_off
		self.sig_0 = sig_0
		self.sig_1 = sig_1
		self.sig_2 = sig_2
		self.ks_d = ks_d
		self.ks_p = ks_p
		self.sw_w = sw_w
		self.sw_p = sw_p
		self.Mod_ind = Mod_ind
		self.Mod_indT = Mod_indT
		self.time = time
		self.utc = ""
		if (F1 > 70) or (F2 > 70) or (F3 > 70):
			self.isphonecall = True
		else:
			self.isphonecall = False
	def str_fmt(self,pred_class):
		#pred_class is either RFI or PULSES, as determined by the classifier
		#pred_class is a str
		out_str = str([self.beam,self.sample,self.sn,self.dm,self.box,self.F1,
			self.F2,self.F3,self.event,self.left,self.right,self.event_div,
			self.mean_off,self.std_off,self.sig_0,self.sig_1,self.sig_2,
			self.ks_d,self.ks_p,self.sw_w,self.sw_p,self.Mod_ind,
			self.Mod_indT,self.time,self.utc]).strip("[]").replace(", "," ")+\
					"\n"
		return str(prec_class)+" "+out_str


"beam sample sn dm box F1 F2 F3 event left right event_div "+\
		            "mean_off std_off sig_0 sig_1 sig_2 ks_d ks_p sw_w sw_p "+\
					            "Mod_ind Mod_indT time utc\n"

fil = FilReader("/home/wfarah/highres_1644/2016-11-10-04:27:01/FB/BEAM_177/2016-11-10-04:27:01.fil")
fbottom = fil.header.ftop
foff = fil.header.foff

f1_min = 840
f1_max = 845
f2_min = 835
f2_max = 840
f3_min = 830
f3_max = 835
F1_ch = []
F2_ch = []
F3_ch = []
F_rest_ch = []
for i in range(320):
	f = fbottom + i*foff
	if f<f1_max and f>f1_min:
		F1_ch.append(i)
	elif f<f2_max and f>f2_min:
		F2_ch.append(i)
	elif f<f3_max and f>f3_min:
		F3_ch.append(i)
	else:
		F_rest_ch.append(i)


###############################################################################
#
# Turn the calling process into a daemon
#
def daemonize(pidfile, logfile):
# standard input will always be directed to /dev/null
	stdin = "/dev/null"
	stdout = logfile
	stderr = logfile

	try:
		pid = os.fork()
		if pid > 0:
			# exit first parent
			sys.exit("Fork #1")
	except OSError, e:
		sys.stderr.write("fork #1 failed: %d (%s)\n" % (e.errno, e.strerror))
		sys.exit(1)

# decouple from parent environment
	os.chdir("/")
	os.setsid()
	os.umask(0)

# do second fork
	try:
		pid = os.fork()
		if pid > 0:
# exit from second parent
			sys.exit("Fork #2")
	except OSError, e:
		sys.stderr.write("fork #2 failed: %d (%s)\n" % (e.errno, e.strerror))
		sys.exit(1)

# redirect standard file descriptors
	sys.stdout.flush()
	sys.stderr.flush()
	si = file(stdin, 'r')
	so = file(stdout, 'a+')
	se = file(stderr, 'a+', 0)
	os.dup2(si.fileno(), sys.stdin.fileno())
	os.dup2(so.fileno(), sys.stdout.fileno())
	os.dup2(se.fileno(), sys.stderr.fileno())
	atexit.register(delpid,pidfile)
# write pidfile, enable a function to cleanup pid file upon crash
	pid = str(os.getpid())
	logging.debug("Writing pid file (pid %s)",pid)
	file(pidfile,'w+').write("%s\n" % pid)




def parse_cfg(cfg_file,tags=None):
	"""Function that returns config file with given tags as dictionar

	Args:
		cfg_file (str): full directory to config file
		tags (list): list of tags to search the cgf_file

	Returns:
		config_dict (dict): dictionary with keys given in tags, and values
							extracted from cfg_file. If one tag doesn't exist,
							value corresponded will be None, else value is of
							type str, or list if multiple values exist for 
							same key.
	"""
	if tags == None:
		tags = []
		with open(cfg_file) as o:
			for line in o:
				if line[0] in ["\n","#"]: continue
				tags.append(line.split()[0])
	config_dict = {}
	with open(cfg_file) as o:
		for line in o:
			if line[0] in ["\n","#"]: continue
			for tag in tags:
				if tag in line:
					i = line.split()
					assert tag == i[0]
					config_dict[tag] = []
					for ii in i[1:]:
						if "#" in ii: break
						config_dict[tag].append(ii)
					if len(config_dict[tag]) == 1:
						config_dict[tag] = config_dict[tag][0]
					tags.remove(tag)
	for tag in tags:
		logging.warning("Couldn't parse <"+tag+"> from "+cfg_file)
		config_dict[tag] = None
	return config_dict


def control_monitor(control_dir,script_name):
	""" Function that writes pid into control folder, and constantly lists
		directory for .quit file, and kills the script
		
	Args:
		control_dir (str): control directory specified in mopsr.cfg
		script_name (str): name of script as it should show in ctrl direc
							(without extension)
		pid (int): process id

	"""
	while True:
		time.sleep(0.5)
		lst = os.listdir(control_dir)
		if script_name+".quit" in lst:
			logging.critical("Read .quit file, cleaning and exiting")
			stop_daemon(control_dir+"/"+script_name+".pid")


def client_control_monitor(control_dir,script_name,bf_numb):
	""" Same as control_monitor, but edited to check for the quit file
	    without the bp number as suffix """
	script_name_suffix = script_name+"_"+str(bf_numb)
	while True:
		time.sleep(0.5)
		lst = os.listdir(control_dir)
		if script_name+".quit" in lst or\
				script_name_suffix+".quit" in lst:
			logging.critical("Read .quit file, cleaning and exiting")
			stop_daemon(control_dir+"/"+script_name_suffix+".pid")


def delpid(pidfile):
	logging.critical("Deleting pidfile")
	os.remove(pidfile)

def stop_daemon(pidfile):
	pf = file(pidfile,'r')
	pid = int(pf.read().rstrip())
	pf.close()
	logging.critical("Trying to kill main thread %s ",pid)
	os.kill(pid,SIGTERM)
	time.sleep(4)
	logging.critical("Main thread %s is still alive, sending SIGKILL",pid)
	os.kill(pid,SIGKILL)

def sigHandler(signo,frame):
	logging.critical("%s Recieved a SIGTERM, cleaning...",os.getpid())
	sys.exit(1)

def create_xml_elem(msg_type,dump_dict=None):
	if msg_type == "ok":
		ok_tag = Element('frb_detector_message')
		child = SubElement(ok_tag, 'reply')
		child.text = "ok"
		return tostring(ok_tag,encoding ='ISO-8859-1').replace("\n","")
	elif msg_type == "fail":
		fail_tag = Element('frb_detector_message')
		child = SubElement(fail_tag, 'reply')
		child.text = "fail"
		return tostring(fail_tag,encoding='utf8').replace("\n","")
	elif msg_type == "dump":
		dump_tag = Element('frb_detector_message')
		cmd = SubElement(dump_tag,'cmd')
		cmd.text = 'dump'
		start_utc = SubElement(dump_tag,'start_utc')
		start_utc.text = dump_dict['start_utc']
		end_utc = SubElement(dump_tag,'end_utc')
		end_utc.text = dump_dict['end_utc']
		dm = SubElement(dump_tag,'dm')
		dm.text = dump_dict['dm']
		beam_number = SubElement(dump_tag,'beam_number')
		beam_number.text = dump_dict['beam_number']
		utc = SubElement(dump_tag,'utc')
		utc.text = dump_dict['utc']
		snr = SubElement(dump_tag,'snr')
		snr.text = dump_dict['snr']
		probability = SubElement(dump_tag,'probability')
		probability.text = dump_dict['probability']
		return tostring(dump_tag,encoding='ISO-8859-1').replace("\n","")

def mod_index(event,t_crunch=False):
	#event should be median subtracted, crunches in time when t_crunch is True
	"""Computes the modulation index of a 2D event window. 
	Crunches time if t_crunch=Trues"""
	if t_crunch:
		event=event.sum(axis=1)
	return np.sqrt((event**2).mean()-(event.mean())**2)/event.mean()


def get_features(beam,t_sample,sn,H_dm,H_w,file_directory):
	timer = time.time()
	fil = FilReader(file_directory)
	t_smear = np.ceil(((31.25*8.3*H_dm)/(0.840)**3)/(fil.header.tsamp*1000000))
	w = 2**H_w
	block = fil.readBlock(int(t_sample-3.5*w),int(t_smear+6.5*w))
	disp_block = block.dedisperse(H_dm)
	t = disp_block.sum(axis=0)
	event = disp_block[:,3*w:4*w]
	event_flatten = event.flatten()
	off_event = np.column_stack((disp_block[:,:3*w],disp_block[:,4*w:]))
	off_event_flatten = off_event.flatten()
	mean_offevent = off_event_flatten.mean()
	std_offevent = off_event_flatten.std()
	sig_0 = np.where(event_flatten>mean_offevent)[0].shape[0]/float(event.size)
	sig_1 = np.where(event_flatten>mean_offevent + std_offevent)[0].shape[0]/float(event.size)
	sig_2 = np.where(event_flatten>mean_offevent + 2*std_offevent)[0].shape[0]/float(event.size)
	event_left = disp_block[:,w:2*w]
	event_right = disp_block[:,5*w:6*w]
	event_s = event.mean(axis=1)
	event_s_mean = event_s.mean()
	event_s_std = event_s.std()
	ks = scipy.stats.kstest(event_s,'norm',[event_s_mean,event_s_std])
	ks_d = ks[0]
	ks_pvalue = ks[1]
	sw_w,sw_pvalue = scipy.stats.shapiro(event_s)
	event_left_s = event_left.mean(axis=1)
	event_right_s = event_right.mean(axis=1)
	med = np.median(block)
	mod_ind = mod_index(event-med)
	mod_indT = mod_index(event-med,True)
	m = (event_s/event_left_s + event_s/event_right_s).mean()/2
	event_med_sub = event_s - med
	all_ch = event_med_sub.sum()
	F1 = 100*event_med_sub[F1_ch].sum() / all_ch
	F2 = 100*event_med_sub[F2_ch].sum() / all_ch
	F3 = 100*event_med_sub[F3_ch].sum() / all_ch
	timer = time.time() - timer
	ftrs = Features(beam=int(beam),sample=int(t_sample),sn=sn,dm=H_dm,box=H_w,
			F1=float(F1),F2=float(F2),F3=float(F3),
			event=float(np.sum(event_s-med)),
			left=float(np.sum(event_left_s-med)),
			right=float(np.sum(event_right_s-med)),event_div=float(m),
			mean_off=float(mean_offevent),std_off=float(std_offevent),
			sig_0=float(sig_0),sig_1=float(sig_1),sig_2=float(sig_2),ks_d=ks_d,
			ks_p=ks_pvalue,sw_w=sw_w,sw_p=sw_pvalue,Mod_ind=float(mod_ind),
			Mod_indT=float(mod_indT),time=timer)
	return ftrs



def get_feature_names():
	return "beam sample sn dm box F1 F2 F3 event left right event_div "+\
			"mean_off std_off sig_0 sig_1 sig_2 ks_d ks_p sw_w sw_p "+\
			"Mod_ind Mod_indT time utc\n"
