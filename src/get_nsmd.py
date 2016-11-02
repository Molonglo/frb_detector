import os
from subprocess import Popen
import socket
import time

def spawn_it_once(port):
	process = Popen(args=["./get_nsmd",str(port)],shell = False)
	time.sleep(0.05)
	return process

def get_nsmd(utc,ra,dec,port):
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	addr = ("localhost",port)
	s.connect(addr)
	s.sendall(utc+","+ra+","+dec)
	ns,md = [float(i) for i in (s.recv(32).rstrip("\x00").split(","))]
	s.close()
	return ns,md

if __name__=="__main__":
	port = 30022 #Port number
	process = spawn_it_once(port)
	
	utc = "2016-09-14-11:43:54.3"
	ra = "03:45:43.1"
	dec = "-34:54:12.43"

	ns,md = get_nsmd(utc,ra,dec,port)
	print ns,md

	process.kill()
