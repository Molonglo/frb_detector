/***************************************************************************
 *  
 *    Copyright (C) 2013 by Andrew Jameson
 *    Licensed under the Academic Free License version 2.1
 * 
 ****************************************************************************/

#include "/home/dada/psrdada/mopsr/src/mopsr_delays.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>



float getmd(char * utc_start_str, char * ra, char * dec)
{
  int arg = 0;

  float fractional = 0;

  char negative_dec = 0;

  // convert UTC_START to a unix UTC
  time_t utc_start = str2utctime (utc_start_str);

  mopsr_source_t source;

  mopsr_delays_hhmmss_to_rad(ra,&(source.raj));
  mopsr_delays_ddmmss_to_rad (dec, &(source.decj));
  
  struct timeval timestamp;
  timestamp.tv_sec = utc_start;
  timestamp.tv_usec = (uint64_t) (fractional * 1e6);
  
  struct tm * utc = gmtime (&utc_start);
  cal_app_pos_iau (source.raj, source.decj, utc,
                   &(source.ra_curr), &(source.dec_curr));

  double jer_delay = calc_jer_delay (source.ra_curr, source.dec_curr, timestamp);
  double md_angle = asin(jer_delay) * (180.0 / M_PI);

  return md_angle;
}



float getns(char * utc_start_str, char * ra, char * dec) 
{
  int arg = 0;

  int device = 0;

  float fractional = 0;

  char negative_dec = 0;

  // convert UTC_START to a unix UTC
  time_t utc_start = str2utctime (utc_start_str);

  mopsr_source_t source;

  mopsr_delays_hhmmss_to_rad(ra,&(source.raj));
  mopsr_delays_ddmmss_to_rad (dec, &(source.decj));

  struct timeval timestamp;
  timestamp.tv_sec = utc_start;
  timestamp.tv_usec = (uint64_t) (fractional * 100000);
  
  struct tm * utc = gmtime (&utc_start);
  cal_app_pos_iau (source.raj, source.decj, utc,
                   &(source.ra_curr), &(source.dec_curr));


  double jer_delay = calc_jer_delay (source.ra_curr, source.dec_curr, timestamp);
  double md_angle = asin(jer_delay);
  double ha_source = calc_ha_source (source.ra_curr, source.dec_curr, timestamp);

  double ns_angle = ns_tilt(ha_source, source.dec_curr, md_angle);


  ns_angle *= (180.0 / M_PI);

  return ns_angle;
}

void error(const char *msg)
{
	perror(msg);
	exit(1);
}

int main(int argc, char ** argv){
  //char * utc = "2016-09-05-11:43:54";
  //char * ra = "11:23:54.3";
  //char * dec = "-43:23:43";
  char * port = argv[1];
  char * utc;
  char * ra;
  char * dec;

  utc = malloc(sizeof(char)*16);
  ra = malloc(sizeof(char)*16);
  dec = malloc(sizeof(char)*16);

  double md,ns;
  int sockfd, newsockfd, portno;
  socklen_t clilen;
  char buffer[64];
  char out[32];
  struct sockaddr_in serv_addr, cli_addr;
  int n;
  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd <0)
	  error("ERROR opening socket");
  if (setsockopt(sockfd, SOL_SOCKET, SO_REUSEADDR, &(int){ 1 }, sizeof(int)) < 0)
	  error("setsockopt(SO_REUSEADDR) failed");
  bzero((char *) &serv_addr, sizeof(serv_addr));

  portno = atoi(port);
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  serv_addr.sin_port = htons(portno);
  if (bind(sockfd, (struct sockaddr *) &serv_addr,
			  sizeof(serv_addr)) < 0)
	  error("ERROR on binding");
  listen(sockfd,5);
  clilen = sizeof(cli_addr);
  while (1){
    newsockfd = accept(sockfd, 
  			(struct sockaddr *) &cli_addr,
		  &clilen);
  	if (newsockfd < 0) 
  		error("ERROR on accept");
  	bzero(buffer,64);
	bzero(out,sizeof(out));
  	n = read(newsockfd,buffer,64);
  	if (n < 0) error("ERROR reading from socket");



  	utc = strtok(buffer,",");
  	ra = strtok(NULL,",");
  	dec = strtok(NULL,",");

  	md = getmd(&utc[0],&ra[0],&dec[0]);
  	ns = getns(&utc[0],&ra[0],&dec[0]);

	snprintf(out,sizeof(out),"%.10f,%.10f",ns,md);
	n = write(newsockfd,out,sizeof(out));
	close(newsockfd);


  	bzero(utc, 16);
  	bzero(ra, 16);
  	bzero(dec, 16);

  }
  close(sockfd);
  free(utc);
  free(ra);
  free(dec);
  return 0;
}
