#include "functions.h"

#define MAX(a,b)((a>b)?a:b)
#define version 1.01

struct BandPower {
    int F1;
    int F2;
    int F3;
    int F_rst;
    int F_all;
    float R1;
    float R2;
    float R3;
};

struct CandFeatures {
    float sn;
    float sn_0;
    int width;
    int nstart;
    int nend;
    float F1;
    float F2;
    float F3;
    int n_rms_mask;
    float sn_rms;
    float mod_ind;
    float mod_indT;
	int isphonecall;
};

const int F1[7] = {7,8,9,10,11,12,-99};
const int F2[6] = {14,15,16,17,18,-99};
const int F3[6] = {20,21,22,23,24,-99};
const int F_rst[25] = {0,1,2,3,4,5,6,13,19,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,-99};

const int widths_1[7] = {1,2,3,4,-99,-99,-99};
const int widths_2[7] = {2,3,4,5,6,7,8};
const int widths_3[7] = {4,6,8,10,12,14,16};
const int widths_4[7] = {8,12,16,20,24,28,32};
const int widths_5[7] = {16,24,32,40,48,56,64};



/* Function to be implemented to update the variables in the
 * header files */
/*
void get_header(){
}
*/

/* Function that returns the median of a array of long */
static inline float median(long x[], int n) {
    int i,j;
    long * temp_arr;
    temp_arr = malloc(n*sizeof(long));
    for (i=0;i<n;i++) temp_arr[i]=x[i];
    float temp;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++){
        for(j=i+1; j<n; j++) {
            if(temp_arr[j] < temp_arr[i]) {
                // swap elements
                temp = temp_arr[i];
                temp_arr[i] = temp_arr[j];
                temp_arr[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        temp = (temp_arr[n/2] + temp_arr[n/2 - 1]) / 2.0;
        free(temp_arr);
        return(temp);
    } else {
        // else return the element in the middle
        temp = x[n/2];
        free(temp_arr);
        return temp;
    }
}

/* Function that returns the median of an array of floats */
static inline float median_float(float x[], int n){
    int i,j;
    float * temp_arr;
    temp_arr = malloc(n*sizeof(float));
    for (i=0;i<n;i++) temp_arr[i]=x[i];
    float temp;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++){
        for(j=i+1; j<n; j++) {
            if(temp_arr[j] < temp_arr[i]) {
                // swap elements
                temp = temp_arr[i];
                temp_arr[i] = temp_arr[j];
                temp_arr[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        temp = (temp_arr[n/2] + temp_arr[n/2 - 1]) / 2.0;
        free(temp_arr);
        return(temp);
    } else {
        // else return the element in the middle
        temp = x[n/2];
        free(temp_arr);
        return temp;
    }
}

static inline float median_char(unsigned char x[], int n){
    int i,j;
    unsigned char * temp_arr;
    temp_arr = malloc(n*sizeof(char));
    for (i=0;i<n;i++) temp_arr[i]=x[i];
    float temp;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++){
        for(j=i+1; j<n; j++) {
            if(temp_arr[j] < temp_arr[i]) {
                // swap elements
                temp = temp_arr[i];
                temp_arr[i] = temp_arr[j];
                temp_arr[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        temp = (temp_arr[n/2] + temp_arr[n/2 - 1]) / 2.0;
        free(temp_arr);
        return(temp);
    } else {
        // else return the element in the middle
        temp = x[n/2];
        free(temp_arr);
        return temp;
    }
}



/*  Returns the median and the MAD of 1D array */
static inline void MAD_1D(float time_series[], int n, float *md, float *mad){
    int i;
    float t_median = median_float(time_series,n);
    float * temp;
    temp = malloc(n * sizeof(float));
    for (i=0;i<n;i++){
        temp[i] = fabs(t_median - time_series[i]);
        //printf("%f ",temp[i]);
    }
    *md = t_median;
    *mad = median_float(temp,n);
    free(temp);
}

static inline float get_mad(float array[], int nsamp, float med){
    int i;
    float * temp;
    float mad;
    temp = malloc(nsamp * sizeof(float));
    //printf("med: %f \n",med);
    for (i=0;i<nsamp;i++){
        temp[i] = fabs(array[i] - med);
        //printf("%f ",temp[i]);
    }
    mad = median_float(temp,nsamp);
    free(temp);
    return mad;
}



static inline void time_to_extract(int H_w, float H_dm, int *time_extract,
        int *disp_delay, float sampling_time){
    *disp_delay = ceil( ((31.25*0.0000083*H_dm)/pow(0.840,3))/sampling_time );
    int w_samp = pow(2,H_w);
    *time_extract = MAX(min_after_event+backstep,(*disp_delay+w_samp)*2+backstep);
}


/* Removes baseline of time_series and normalizes it, time_series_norm is casted to float 
void norm_baseline(long * time_series, int size, float med, float mad, float * time_series_norm){
    int i;
    for (i=0;i<size;i++) time_series_norm[i] = (float)((time_series[i] - med)/(1.4826*mad));
}*/

static inline void norm_baseline(float * time_series, int nsamp, int med, float mad){
    int i;
    for (i=0;i<nsamp;i++){
        time_series[i] = ((time_series[i] - med)/(1.4826*mad));
    }
}

/* Convolves 1d array with boxcar of width = width, from start till start + bin search */
static inline void convolve(float * time_series, int ndat, int width,
       int nbin_start, int nbin_search, float * ans){
    int i,j;
    for (i=0;i<ndat;i++) ans[i]=0.0;
    for (i=nbin_start;i<nbin_start+nbin_search;i++){
        float ss=0.0;
        for (j=0;j<width;j++) ss+=time_series[i+j];
        ans[i]=ss;
    }
}



/* Function that returns list of dispersion delays */
static inline void get_dm_delay(struct Header *header, float dm, 
        int *dm_delay){
  float nu1,nu2,nu1_2,nu2_2,shift;
  int i;
  for (i=0; i<header->nchans; i++){
      nu1 = header->fch1;
      nu2 = header->fch1 + i*header->foff; // df usually -ve
      nu1_2 = 1.0e6/(nu1*nu1);
      nu2_2 = 1.0e6/(nu2*nu2);
      shift = 4.148808e-3 * dm * (nu2_2-nu1_2);  // secs
      dm_delay[i] = round(shift / header->tsamp);
  }
  return ;
}

/* Function that reads filterbank files, skips nsamp_skip bytes and reads nsamp_read number 
 * samples. data is the pointer to the 1D data, and is column major.
 */
static inline void read_block(char * file_direc, unsigned long nsamp_skip,
        int nsamp_read, unsigned char * data, struct Header *header){

    unsigned long long bytes_skip = header->totalbytes + 
        nsamp_skip*header->nchans;
    FILE * fptr = fopen(file_direc,"rb");
    if (fptr==NULL){
        printf("Error opening file \n");
        exit(-1);
    }
    fseek(fptr,bytes_skip,SEEK_CUR);
    fread((unsigned char *) data,header->nbits/8,header->nchans*nsamp_read,fptr);
    fclose(fptr);
    return;
}



/* Performs a transformation of index from row-major to column-major */ 
static inline int index_transform(int old_index, int rows, int cols){
    int real_index, old_col, old_row;
    old_row = old_index % cols;
    old_col = old_index / cols;
    real_index = old_row * rows + old_col;
    return real_index;
}
/*
 Computes the time_series 
static inline void get_time_series(unsigned char * data, float dm, long * time_series, int nsamp_read){
    int i,j;
    unsigned long index,real_index,sum;
    int dm_delay[nchans];

    get_dm_delay(f0,df,dm,dm_delay);
    

    for (i=0;i<nsamp_read;i++){
        sum = 0;
        for (j=0;j<nchans;j++){
            index = j*nsamp_read + (i+dm_delay[j])%nsamp_read;
            //if (index >= nsamp_read * nchans){
            //    printf("ERROR, READING OUT OF BOUNDS");
            //    exit(-1);
            //}
            real_index = index_transform(index,nchans,nsamp_read);
            sum += data[real_index];
        }
        time_series[i] = sum;
    }
   return;
}
*/
/*
void get_time_series2(unsigned char * data, float dm, long * time_series, int nsamp_read){
    int i,j;
    unsigned long index,sum;
    int dm_delay[nchans];
    int total = nsamp_read * nchans;
    get_dm_delay(f0,df,dm,dm_delay);

    for (i=0;i<nsamp_read;i++){
        sum = 0;
        for (j=0;j<nchans;j++){
            sum += data[(i*nchans + j + nchans*dm_delay[j])%total];
        }
        time_series[i] = sum;
    }
    return;
}
*/


static inline void get_time_series2(unsigned char * data, int * dm_delay, 
        float * time_series, int nsamp_read, int nchans){
    int i,j;
    int total = nsamp_read * nchans;
    long sum;

    for (i=0;i<nsamp_read;i++){
        sum = 0;
        for (j=0;j<nchans;j++){
            sum += data[(i*nchans + j + nchans*dm_delay[j])%total];
        }
        time_series[i] = sum;
    }
    return;
}


static inline float block_median(unsigned char * block, int nsamp_read, 
        int nchans){
    float * n_medians = malloc(nsamp_read*sizeof(int));
    int i;
    float med_est;
    for (i=0;i<nsamp_read;i++){
        n_medians[i] = median_char(&block[i*nchans],nchans);
    }
    med_est = median_float(n_medians,nsamp_read);
    free(n_medians);
    return med_est;
}

/* Function that returns a boolean if integer is in array */
static inline bool inarray(int j, const int* F_n){
    int i;
    while (F_n[i]!=-99){
        if (j==F_n[i]) return true;
        i+=1;
    }
    return false;
}

/* Returns an int telling in which frequency band the int j is*/
static inline int get_freq_band(int j){
    if (inarray(j,F1)) return 1;
    if (inarray(j,F2)) return 2;
    if (inarray(j,F3)) return 3;
    return 0;
}

/* Initializing the struct power in bands to 0 */
static inline void init_struct(struct BandPower *freq_chan){
    freq_chan->F1 = 0;
    freq_chan->F2 = 0;
    freq_chan->F3 = 0;
    freq_chan->F_rst = 0;
    freq_chan->F_all = 0;
    return ;
}

/* Initializing the struct CandFeatures */
static inline void init_features(struct CandFeatures *features){
	float sn = -1.;
	float sn_0 = -1.;
    int width = -1;
    int nstart = -1;
    int nend = -1;
    float F1 = -1.;
    float F2 = -1.;
    float F3 = -1.;
    int n_rms_mask = -1;
    float sn_rms = -1.;
    float mod_ind = -1.;
    float mod_indT = -1.;
    int isphonecall = -1;
}

static inline int find_max(float * array, int n){
    int c, index = 0;
    float max;

    max = array[0];

    for (c=1; c<n; c++){
        if (array[c] > max){
            index = c;
            max = array[c];
        }
    }
    return index;
}

/* Function that computes the power in each frequency band, stores information
 * in the struct BandPower */
static inline void get_power(struct BandPower *freq_chan, int *dm_delay,
        int nchans, unsigned char *block, int start_bin, 
        int end_bin, int median){
    int temp;
    int i,j;
    int which_freq_band;

    //printf("%i %i\n\n\n",start_bin,end_bin);
    for (j=0;j<nchans;j++){
        temp = 0;
        for (i=start_bin;i<end_bin;i++){
            temp += (int)block[i*nchans+j+dm_delay[j]*nchans] - median;
        }
        which_freq_band = get_freq_band(j);
        switch(which_freq_band){
            case 1:
                freq_chan->F1 += temp;
                break;
            case 2:
                freq_chan->F2 += temp;
                break;
            case 3:
                freq_chan->F3 += temp;
                break;
            case 0:
                freq_chan->F_rst += temp;
                break;
        }
        freq_chan->F_all += temp;
    }
    //printf("%i %i %i %i %i\n",freq_chan->F1,freq_chan->F2,freq_chan->F3,freq_chan->F_rst,freq_chan->F_all);
    freq_chan->R1 = (float)freq_chan->F1 / freq_chan->F_all;
    freq_chan->R2 = (float)freq_chan->F2 / freq_chan->F_all;
    freq_chan->R3 = (float)freq_chan->F3 / freq_chan->F_all;
    //printf("\n%f %f %f\n",freq_chan->R1,freq_chan->R2,freq_chan->R3);
}

/* Function that returns whether the candidate pass the phone call test */
static inline bool is_phone_call(struct BandPower *phone_bands){
    if (phone_bands->R1 > 0.7 || phone_bands->R2 > 0.7 || phone_bands->R3 > 0.7
            || phone_bands->R1 + phone_bands->R2 > 0.9 || phone_bands->R2 + 
            phone_bands->R3 > 0.9 || phone_bands->R1 + phone_bands->R3 > 0.9){
        return true;
    }
    else return false;
}

static inline float get_std(unsigned char * array, int size){
    float mean = 0.0, sum_deviation = 0.0;
    int i;
    for (i=0; i<size; i++){
        mean += array[i];
    }
    mean = mean/size;
    for (i=0; i<size; i++) sum_deviation += (array[i]-mean)*(array[i]-mean);
    return sqrt(sum_deviation/size);
}

static inline void compute_rms(unsigned char * block, int nsamp, int * dm_delay, 
        int nstart, int nend, float * rms_freq, int nchans,
        int * event){
    int nsamp_noevent = nsamp-(nend-nstart);
    unsigned char * freq_channel;
    freq_channel = malloc(nsamp_noevent*sizeof(int));
    int i,ii,j,jj;
    int low_lim, high_lim;

    ii = 0;
    for (i=0;i<nchans;i++){
        low_lim = nstart+dm_delay[i];
        high_lim = nend+dm_delay[i];
        jj=0;
        for (j=0;j<nsamp;j++){
            if (j>=low_lim && j<high_lim){
                //Warning 'event' is row-major
                event[ii] = (int)block[i+j*nchans];
                ii+=1;
                continue;
            }
            freq_channel[jj] = block[i+j*nchans];
            jj+=1;
        }
        rms_freq[i] = get_std(freq_channel,nsamp_noevent);
    }
    free(freq_channel);
    return;
}

static inline float get_mean(int * array, int size){
    int i;
    float mean = 0;
    for (i=0;i<size;i++) mean += array[i];
    mean = mean / size;
    return mean;
}

static inline void get_mod_ind(int * event, int nchans, int width, 
        float * mod_ind, float * mod_indT){
    float mean_event = 0;
    float mean_event_squared = 0;
    int total_samp = width * nchans;
    int * t_crunched;
    t_crunched = malloc(nchans*sizeof(int));
    int * event_squared;
    event_squared = malloc(total_samp*sizeof(int));

    int i,j;

    for (i=0;i<nchans;i++) t_crunched[i]=0;

    for (i=0;i<total_samp;i++){
        event_squared[i] = event[i]*event[i];
        mean_event += event[i];
    }
    j=0;
    for (i=0;i<total_samp;i++){
        mean_event_squared += event_squared[i];
        t_crunched[j] += event[i];
        if ((i+1)%width==0) j++;
    }
    mean_event = mean_event/total_samp;
    mean_event_squared = mean_event_squared/total_samp;

    *mod_ind = sqrt(mean_event_squared-mean_event*mean_event)/mean_event;

    int * t_crunched_squared;
    t_crunched_squared = malloc(nchans*sizeof(int));

    float t_crunched_mean = 0;
    float t_crunched_mean_squared = 0;

    for (i=0;i<nchans;i++){
        t_crunched_squared[i] = t_crunched[i] * t_crunched[i];
        t_crunched_mean += t_crunched[i];
        //printf("%i ",t_crunched[i]);
    }
    for (i=0;i<nchans;i++) t_crunched_mean_squared += t_crunched_squared[i];

    t_crunched_mean = t_crunched_mean / nchans;
    t_crunched_mean_squared = t_crunched_mean_squared / nchans;

    *mod_indT = sqrt(t_crunched_mean_squared - t_crunched_mean*
            t_crunched_mean)/t_crunched_mean;

    free(t_crunched_squared);
    free(event_squared);
    free(t_crunched);
    return;
}


static inline bool is_in_array(int x, int *array, int size){
    int i;
    for (i=0;i<size;i++){
        if (x==array[i]) return true;
    }
    return false;
}

static inline float get_snr_rms(unsigned char * block, float * time_series, int nsamp_read,
        int * dm_delay, int * contam_chans, int n_contam_chans, int nstart, 
        int width, int nchans){
    int total = nsamp_read * nchans;
    int i,j;
    float median_ts,mad,sn_rms;

    for (i=0;i<nchans;i++){
        if (is_in_array(i,contam_chans,n_contam_chans)){
            for (j=0;j<nsamp_read;j++){
                time_series[j] -= block[(j*nchans + i + nchans*dm_delay[i])%total];
            }
        }
    }

    MAD_1D(time_series,nsamp_read,&median_ts,&mad);
    norm_baseline(time_series, nsamp_read, median_ts, mad);
    sn_rms = 0;
    for (i=nstart;i<nstart+width;i++) sn_rms += time_series[i];
    sn_rms = sn_rms / sqrt(width);
    return sn_rms;
}


static inline void HUNTER(unsigned char * block, float * time_series_raw, int nsamp, 
        int * dm_delay, int H_w, int disp_delay, int nchans,
        struct CandFeatures * features){

    int i, max_index_conv;
    int n_w_trials = 7;
    const int * width_trials;
    float sn_new;
    //float sn[7] = {0.,0.,0.,0.,0.,0.,0.};
    //int nstart[7] = {0,0,0,0,0,0,0};
    int * dm_0;

    float * time_series;
    time_series = malloc(nsamp*sizeof(float));
    float * convolved;
    convolved = malloc(nsamp * sizeof(float));

    float median_ts,mad;
    switch(H_w){
        case 1:;
            width_trials = &widths_1[0];
            break;
        case 2:;
            width_trials = &widths_2[0];
            break;
        case 3:;
            width_trials = &widths_3[0];
            break;
        case 4:;
            width_trials = &widths_4[0];
            break;
        case 5:;
            width_trials = &widths_5[0];
            break;
    }

    for (i=0;i<nsamp;i++) time_series[i] = time_series_raw[i];
    MAD_1D(time_series,nsamp,&median_ts,&mad);

    //printf("\n\n%f\n%f\n\n",mad,median_ts);
    norm_baseline(time_series, nsamp, median_ts, mad);

    features->sn = 0.;

    for (i=0;i<n_w_trials;i++){
        if (width_trials[i]==-99) continue;
        convolve(time_series, nsamp, width_trials[i], 
        backstep - 50, 100, convolved); //the 50 and 100 is to search +- 100
                                        //samples around HEIMDAL start sample
        max_index_conv = find_max(&convolved[backstep - 50], 100);
        sn_new = convolved[backstep - 50 + max_index_conv] / sqrt(width_trials[i]);
        if (sn_new > features->sn){
            features->sn = sn_new;
            features->nstart = max_index_conv + 50;
            features->nend = features->nstart + width_trials[i];
            features->width = width_trials[i];
        }
    }

    //Getting sn_0:
    dm_0 = malloc(nchans*sizeof(int));
    for (i=0;i<nchans;i++) dm_0[i]=0;
    for (i=0;i<nsamp;i++){
        time_series[i] = 0.;
        convolved[i] = 0.;
    }
    get_time_series2(block, dm_0, time_series, nsamp, nchans);
    //printf("nsamples: %i \n",nsamp);
    MAD_1D(time_series,nsamp,&median_ts,&mad);
    norm_baseline(time_series, nsamp, median_ts, mad);
    if (disp_delay == 0) printf("Warning, dispersion delay = 0\n");
    convolve(time_series, nsamp, features->width,
            backstep, disp_delay, convolved);
    max_index_conv = find_max(&convolved[backstep],disp_delay);
    features->sn_0 = convolved[backstep + max_index_conv] / sqrt(features->width);

    free(time_series);
    free(convolved);
    free(dm_0);
    return;
} 


static inline struct CandFeatures process_candidate(int sample, float H_dm, 
        int H_w, char * fil_direc){
    struct CandFeatures features;
	init_features(&features);
    struct Header header;
    int i,ii,total,nsamp_read, disp_delay, median_est;
    unsigned char * block;
    //int dm_delay[nchans]; 
    int event_size;
    struct BandPower freq_chan;
    init_struct(&freq_chan);
    header = get_header(fil_direc);

    int * dm_delay;
    dm_delay = malloc(header.nchans*sizeof(int));

    time_to_extract(H_w, H_dm, &nsamp_read, &disp_delay, header.tsamp);
    get_dm_delay(&header,H_dm,dm_delay);

    total = nsamp_read*header.nchans;
    block = malloc(header.nbits*total/8);

    read_block(fil_direc,sample-backstep,nsamp_read,block,&header);
    median_est = round(block_median(block,nsamp_read,header.nchans));

    int start_bin = backstep - pow(2,H_w)/2;
    int end_bin = backstep + pow(2,H_w)/2;
    get_power(&freq_chan, dm_delay, header.nchans, block, start_bin, 
            end_bin,median_est);

    if (is_phone_call(&freq_chan)){
		features.F1 = freq_chan.R1;
		features.F2 = freq_chan.R2;
		features.F3 = freq_chan.R3;
        free(block);
        free(dm_delay);
		features.isphonecall = 1;
        return features;
    }

	features.isphonecall = 0;

    float rms_freq_median, rms_freq_mad;
    int n_contam_chans = 0;
    float * rms_freq;
    int * contam_chans;
    rms_freq = malloc(header.nchans*sizeof(float));
    float * time_series_raw;
    time_series_raw = malloc(nsamp_read*sizeof(float));
    int * event;

    get_time_series2(block, dm_delay, time_series_raw, 
            nsamp_read, header.nchans);

    HUNTER(block, time_series_raw, nsamp_read, dm_delay, 
            H_w, disp_delay, header.nchans, &features);


    event_size = header.nchans*features.width;
    event = malloc(event_size*sizeof(int));
    
    compute_rms(block, nsamp_read, dm_delay, features.nstart, features.nend, 
            rms_freq, header.nchans, event);

    rms_freq_median = median_float(rms_freq, header.nchans);
    rms_freq_mad = get_mad(rms_freq, header.nchans, rms_freq_median);
    
    for (i=0;i<header.nchans;i++){
        rms_freq[i] = (rms_freq[i] - rms_freq_median)/rms_freq_mad;
        if (rms_freq[i] > 3) n_contam_chans += 1;
    }

    contam_chans = malloc(n_contam_chans * sizeof(int));
    ii=0;
    for (i=0;i<header.nchans;i++){
        if (rms_freq[i]>3){
            contam_chans[ii] = i;
            ii+=1;
        }
    }

    features.n_rms_mask = n_contam_chans;
    if (n_contam_chans == 0){
        features.sn_rms = features.sn;
    }
    else{
        features.sn_rms = get_snr_rms(block, time_series_raw, nsamp_read, 
                dm_delay, contam_chans, n_contam_chans, 
                features.nstart, features.width, header.nchans);
    }
    //printf("\nPrinting event:\n");
    if (features.sn_rms==-99) features.sn_rms = features.sn;
    for (i=0;i<event_size;i++){
        event[i] -= median_est;
    }
    float mod_ind = 0.;
    float mod_indT = 0.;
    get_mod_ind(event,header.nchans,features.width,&mod_ind,&mod_indT);
    features.mod_ind = mod_ind;
    features.mod_indT = mod_indT;
    
    struct BandPower event_band_power;
    init_struct(&event_band_power);

    get_power(&event_band_power, dm_delay, header.nchans, block, 
            features.nstart, features.nend,median_est);

    features.F1 = event_band_power.R1;
    features.F2 = event_band_power.R2;
    features.F3 = event_band_power.R3;

    free(contam_chans);
    free(rms_freq);
    free(time_series_raw);
    free(block);
    free(event);
    free(dm_delay);
    return features;
}


void get_version(){
	printf("Version: %.2f",version);
}


/* Link to python module */
struct CandFeatures get_features(int sample, float H_dm, int H_w, 
        char * file_direc){
    //printf("Sample: %i\n",sample);
    //printf("DM: %f\n",H_dm);
    //printf("Width: %i\n",H_w);
    //printf("directory: %s",file_direc);
    struct CandFeatures features;
    features = process_candidate(sample, H_dm, H_w, file_direc);
    return features;
}

/*
void test_header(char * file_direc){
    struct Header current_header; 
    
    current_header = get_header(file_direc);
    printf("Total header bytes: %i\n",current_header.totalbytes);
    return;
}
*/

int main(){
    //char *file_direc = "/home/wfarah/disp_test/2016-08-05-07:48:39.fil";
	char * file_direc = "/data/mopsr/archives/2016-09-04-16:12:06/BEAM_034/2016-09-04-16:12:06.fil";
    struct CandFeatures features;
    //test_header(file_direc);
    //features = get_features(51072,58.1129,3,file_direc);
    features = get_features(1924753, 69.042,5,file_direc);
	printf("\n\n\n\nsn: %f\nsn_0: %f\nwidth: %i\nstart_bin: %i\nend bin: %i\n",features.sn,features.sn_0,features.width,features.nstart,features.nend);
    printf("sn_rms: %f\n",features.sn_rms);
    printf("n_rms_mask: %i\n",features.n_rms_mask);
    printf("mod_ind: %f\n",features.mod_ind);
    printf("mod_indT: %f\n",features.mod_indT);
    printf("F1: %f, F2: %f, F3: %f\n",features.F1,features.F2,features.F3);
    return 1;
}
