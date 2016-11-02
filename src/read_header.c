/* read_header.c - general (modified) handling routines for SIGPROC headers */
#include "read_header.h"


/* read a string from the input which looks like nchars-char[1-nchars] */
static inline void get_string(FILE *inputfile, int *nbytes, char string[]) /* includefile */
{
  int nchar;
  strcpy(string,"ERROR");
  fread(&nchar, sizeof(int), 1, inputfile);
  *nbytes=sizeof(int);
  if (feof(inputfile)) exit(0);
  if (nchar>80 || nchar<1) return;
  fread(string, nchar, 1, inputfile);
  string[nchar]='\0';
  *nbytes+=nchar;
}


/* attempt to read in the general header info from a pulsar data file */
struct Header get_header(char *fil_direc) /* includefile */
{
  struct Header header;
  char string[80], message[80];
  header.expecting_rawdatafile=0;
  header.expecting_source_name=0; 
  header.expecting_frequency_table=0;
  int nbytes;
  FILE * inputfile = fopen(fil_direc,"rb");
  if (inputfile==NULL){
	  printf("Error opening file \n");
	  exit(-1);
  }
  /* try to read in the first line of the header */
  get_string(inputfile,&nbytes,string);
  if (strcmp(string,"HEADER_START")!=0) {
	/* the data file is not in standard format, rewind and return */
	rewind(inputfile);
    exit(0);
  }
  /* store total number of bytes read so far */
  header.totalbytes=nbytes;
  /* loop over and read remaining header lines until HEADER_END reached */
  while (1) {
    get_string(inputfile,&nbytes,string);
    if (strcmp(string,"HEADER_END")==0) break;
    header.totalbytes+=nbytes;
    if (strcmp(string,"rawdatafile")==0) {
      header.expecting_rawdatafile=1;
    } else if (strcmp(string,"source_name")==0) {
      header.expecting_source_name=1;
    } else if (strcmp(string,"FREQUENCY_START")==0) {
      header.expecting_frequency_table=1;
      header.channel_index=0;
    } else if (strcmp(string,"FREQUENCY_END")==0) {
      header.expecting_frequency_table=0;
    } else if (strcmp(string,"az_start")==0) {
      fread(&header.az_start,sizeof(header.az_start),1,inputfile);
      header.totalbytes+=sizeof(header.az_start);
    } else if (strcmp(string,"za_start")==0) {
      fread(&header.za_start,sizeof(header.za_start),1,inputfile);
      header.totalbytes+=sizeof(header.za_start);
    } else if (strcmp(string,"src_raj")==0) {
      fread(&header.src_raj,sizeof(header.src_raj),1,inputfile);
      header.totalbytes+=sizeof(header.src_raj);
    } else if (strcmp(string,"src_dej")==0) {
      fread(&header.src_dej,sizeof(header.src_dej),1,inputfile);
      header.totalbytes+=sizeof(header.src_dej);
    } else if (strcmp(string,"tstart")==0) {
      fread(&header.tstart,sizeof(header.tstart),1,inputfile);
      header.totalbytes+=sizeof(header.tstart);
    } else if (strcmp(string,"tsamp")==0) {
      fread(&header.tsamp,sizeof(header.tsamp),1,inputfile);
      header.totalbytes+=sizeof(header.tsamp);
    } else if (strcmp(string,"period")==0) {
      fread(&header.period,sizeof(header.period),1,inputfile);
      header.totalbytes+=sizeof(header.period);
    } else if (strcmp(string,"fch1")==0) {
      fread(&header.fch1,sizeof(header.fch1),1,inputfile);
      header.totalbytes+=sizeof(header.fch1);
    } else if (strcmp(string,"fchannel")==0) {
      fread(&header.frequency_table[header.channel_index++],sizeof(double),1,inputfile);
      header.totalbytes+=sizeof(double);
      header.fch1=header.foff=0.0; /* set to 0.0 to signify that a table is in use */
    } else if (strcmp(string,"foff")==0) {
      fread(&header.foff,sizeof(header.foff),1,inputfile);
      header.totalbytes+=sizeof(header.foff);
    } else if (strcmp(string,"nchans")==0) {
      fread(&header.nchans,sizeof(header.nchans),1,inputfile);
      header.totalbytes+=sizeof(header.nchans);
    } else if (strcmp(string,"telescope_id")==0) {
      fread(&header.telescope_id,sizeof(header.telescope_id),1,inputfile);
      header.totalbytes+=sizeof(header.telescope_id);
    } else if (strcmp(string,"machine_id")==0) {
      fread(&header.machine_id,sizeof(header.machine_id),1,inputfile);
      header.totalbytes+=sizeof(header.machine_id);
    } else if (strcmp(string,"data_type")==0) {
      fread(&header.data_type,sizeof(header.data_type),1,inputfile);
      header.totalbytes+=sizeof(header.data_type);
    } else if (strcmp(string,"ibeam")==0) {
      fread(&header.ibeam,sizeof(header.ibeam),1,inputfile);
      header.totalbytes+=sizeof(header.ibeam);
    } else if (strcmp(string,"nbeams")==0) {
      fread(&header.nbeams,sizeof(header.nbeams),1,inputfile);
      header.totalbytes+=sizeof(header.nbeams);
    } else if (strcmp(string,"nbits")==0) {
      fread(&header.nbits,sizeof(header.nbits),1,inputfile);
      header.totalbytes+=sizeof(header.nbits);
    } else if (strcmp(string,"barycentric")==0) {
      fread(&header.barycentric,sizeof(header.barycentric),1,inputfile);
      header.totalbytes+=sizeof(header.barycentric);
    } else if (strcmp(string,"pulsarcentric")==0) {
      fread(&header.pulsarcentric,sizeof(header.pulsarcentric),1,inputfile);
      header.totalbytes+=sizeof(header.pulsarcentric);
    } else if (strcmp(string,"nbins")==0) {
      fread(&header.nbins,sizeof(header.nbins),1,inputfile);
      header.totalbytes+=sizeof(header.nbins);
    } else if (strcmp(string,"nsamples")==0) {
      /* read this one only for backwards compatibility */
      fread(&header.itmp,sizeof(header.itmp),1,inputfile);
      header.totalbytes+=sizeof(header.itmp);
    } else if (strcmp(string,"nifs")==0) {
      fread(&header.nifs,sizeof(header.nifs),1,inputfile);
      header.totalbytes+=sizeof(header.nifs);
    } else if (strcmp(string,"npuls")==0) {
      fread(&header.npuls,sizeof(header.npuls),1,inputfile);
      header.totalbytes+=sizeof(header.npuls);
    } else if (strcmp(string,"refdm")==0) {
      fread(&header.refdm,sizeof(header.refdm),1,inputfile);
      header.totalbytes+=sizeof(header.refdm);
    } else if (header.expecting_rawdatafile) {
      strcpy(header.rawdatafile,string);
      header.expecting_rawdatafile=0;
    } else if (header.expecting_source_name) {
      strcpy(header.source_name,string);
      header.expecting_source_name=0;
    } else {
      sprintf(message,"read_header - unknown parameter: %s\n",string);
      fprintf(stderr,"ERROR: %s\n",message);
      exit(1);
    }
    if (header.totalbytes != ftell(inputfile)){
        fprintf(stderr,"ERROR: Header bytes does not equal file position\n");
        fprintf(stderr,"String was: '%s'\n",string);
        fprintf(stderr,"       header: %d file: %d\n",header.totalbytes,ftell(inputfile));
        exit(1);
    }
  } 
  /* add on last header string */
  header.totalbytes+=nbytes;
  fclose(inputfile);
  /* return total number of bytes read */
  //return totalbytes;
  return header;
}

/*
int main(int argc, char ** argv){
    struct Header current_header;
    FILE * fptr = fopen(argv[1],"rb");
    if (fptr==NULL){
        fprintf(stderr,"Error opening file%s\n",argv[1]);
        exit(-1);
    }

    read_header(fptr,&current_header);
    printf("%i\n",current_header.totalbytes);
    printf("%f\n",current_header.tsamp*1000);
    return 0;
}*/
