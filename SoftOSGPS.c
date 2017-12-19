
#include "SoftOSGPS.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


double samp_rate=16.3676e6; /*  sampling rate */
double carrier_IF=4.1304e6; /*  carrier Intermediate Frequency */

/***********************************************************************
  Software GPS RECEIVER (SoftOSGPS)



***********************************************************************/


#include  "serport.h"           /* NMEA */
#include  "nmea.h"              /* NMEA */
#include  <time.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  "file.h"

#ifdef __linux__
#include <termios.h>
#include <unistd.h>
/*#include <sys/io.h>*/   /* for iopl() */
#include <errno.h>
extern void restore_term (void);
#endif

#define MAIN
#include  "consts.h"
#include  "structs.h"
#include  "interfac.h"
#include  "gpsfuncs.h"
#include  "globals.h"
#include "rinex.h"
#undef  MAIN
FILE *output, *debug, *in, *out, *kalm, *data_bits;
FILE *rinex_obs, *rinex_nav, *rinex_par;
char output_file[40],almanac_file[40],ephemeris_file[40],receiver_file[40];
char location_file[40],ion_utc_file[40],data_bits_file[40],kalman_file[40];
char rinex_obs_file[40],rinex_nav_file[40],rinex_par_file[40],debug_file[40];

/* Internal functions */
static void read_rcvr_par (void);
static void display       (void);
static void chan_allocate (void);
static void cold_allocate (void);
static void read_rinex_par (void);
static void read_filenames (void);

/* External functions */
extern void nav_fix (void);
extern void Interrupt_Install (void);
extern void Interrupt_Remove  (void);

#if (defined VCPP)
#define clear_screen()              _clearscreen (_GCLEARSCREEN);
#define goto_screen_location(x, y)  _settextposition(x, y)
#define check_for_keyboard_press()  (kbhit() ? getch() : '\0')

#elif ((defined BCPP) || (defined __TURBOC__))

#include <conio.h> /* For getch()/kbhit() */
#define clear_screen()              clrscr ()
#define goto_screen_location(x, y)  gotoxy(x, y)
#define check_for_keyboard_press()  (kbhit() ? getch() : '\0')

#elif (defined __linux__)

/* Basic terminal escape sequence (note: ASCII 0x1b == ESC) */
#define clear_screen()              printf ("\n%c[2J", 0x1b);
#define goto_screen_location(x, y)  printf ("%c[%d;%dH", 0x1b, x, y)
#define check_for_keyboard_press()  getchar()

#endif


void generate_prn_codes(void);
void create_wave(float [], int, double, double, double, int);
void gensig(float [], int, int, double, double, double, int);
void Sim_GP2021(char IF[],long nsamp,long *tic);
void gpsisr (void);

char prn_code[33][1023];

int inpwd(unsigned short int add)
{
  return (REG_read[add]);
}

void outpwd(unsigned short int add, unsigned short int data)
{
  REG_write[add]=data;
}

float i_prompt[N_channels],q_prompt[N_channels],i_track[N_channels],q_track[N_channels];
float code_phase[N_channels],carrier_phase[N_channels],old_sin_carrier_phase[N_channels];
int   dump[N_channels],ms_counter[N_channels],bit_counter[N_channels],last_ms_set[N_channels];
long  carrier_cycle[N_channels];
long tic_ref;
int  dump_ref;

struct stat buffer;
int         status;
time_t utctime;

long  tic_counter;

/******************************************************************************
FUNCTION main()
RETURNS  None.

PARAMETERS None.

PURPOSE
	This is the main program to control the GPS receiver

WRITTEN BY
	Clifford Kelley

******************************************************************************/
int
main ( void )
{
	FILE *ifdata;
	long tic;
    char IF[16368];  /* reserve enough space for 1 ms at sampling rate  */
	int i,ch,nsamp;
	float time;
	int wrd;
	char wrdh, wrdl;

        read_filenames();

	out=fopen("test.dat","w+");
	
	time=0.0;
	nsamp = samp_rate*interr_int/1.0e6;
	tic_ref = samp_rate*0.1-1;
    tic = tic_ref;

	generate_prn_codes();

	ifdata=fopen("gnss_RPV_60s.bin","rb");
	/*    ifdata=fopen("prn6_test.log","rb");  */
	if (ifdata == NULL) {
	  fprintf (stderr, "Error: Unable to open IF file: %s\n",
		   strerror (errno));
	  fprintf (stderr, "Exiting...\n");
	  exit (EXIT_FAILURE);
	}

    status = stat("gnss_RPV_60s.bin", &buffer);

	tic_counter=0;

    read_rcvr_par ();

    rec_pos_xyz.x = 0.0;
    rec_pos_xyz.y = 0.0;
    rec_pos_xyz.z = 0.0;

    if (out_kalman == 1)
     kalm = fopen (kalman_file, "w+");
    if (out_rinex == 1)
    {
      rinex_obs = fopen (rinex_obs_file, "w+");
      rinex_nav = fopen (rinex_nav_file, "w+");
      read_rinex_par();
      write_rinex_obs_head = 1;
      write_rinex_nav_head = 1;
    }
    if (out_pos == 1 || out_vel == 1 || out_time == 1)
     output = fopen (output_file, "w+");
    if (out_pos == 1)
     fprintf (output, "time (seconds), latitude (degrees), "
 	     "longitude (degrees), hae (meters), ");
    if (out_vel == 1) 
     fprintf (output, "velocity north, velocity east, velocity up, ");
    if (out_time == 1)
     fprintf (output, "clock offset, ");
    if (out_pos || out_vel || out_time)
     fprintf (output, "hdop, vdop, tdop\n");

    if (out_debug == 1)
      debug = fopen ("debug.log", "w+");
    if (out_data == 1)
      data_bits = fopen ("navdata.log", "w+");
    read_initial_data ();
    current_loc = receiver_loc ();
    rec_pos_llh.lon = current_loc.lon;
    rec_pos_llh.lat = current_loc.lat;
    rec_pos_llh.hae = current_loc.hae;
    nav_tic = nav_up * 10;

    thetime=buffer.st_mtime - 843 + 7201; /* for this data set  */

	/*	utctime = thetime - dtls + 28800; */

	/*    clear_screen ();  */

    if (status != cold_start)
      chan_allocate ();
    else if (status == cold_start)
      cold_allocate ();
    m_time[1] = clock_tow;

  /* Initialize IODE and IODC to the invalid value of -1. */
    for (i = 0; i < 32; i++) {
      gps_eph[i].iode = -1;
      gps_eph[i].iodc = -1;
	}
    read_ephemeris ();
	{
	  /* int err; */
      /*      open_com (0, Com0Baud, 0, 1, 8, &err); */       /* NMEA */
	} 

    for (ch=0;ch<N_channels;ch++)
	{ 
		i_prompt[N_channels]=q_prompt[N_channels]=i_track[N_channels]=q_track[N_channels]=0.0;
        code_phase[N_channels]=carrier_phase[N_channels]=0.0;
		dump[ch]=dump_ref;  /* number of samples per dump */
		carrier_cycle[ch]=0;
		old_sin_carrier_phase[ch]=0.0;
		last_ms_set[ch]=999;
	}


    while (!feof(ifdata))
    { 
	   for (i=0;i<nsamp;i++)   /*  Simulate an interrupt every nsamp-les  */
	   { 
		  IF[i]=fgetc(ifdata);
	   } 
	   time=time+nsamp/samp_rate;
/* 	   printf("time = %f\n",time); */
       Sim_GP2021(IF, nsamp,&tic);
	   gpsisr();   /* Do the interrupt routine stuff (tracking loops etc.)   */
	   display();
       if (data_frame_ready == 1) 
	   {
	      for (ch = 0; ch < N_channels; ch++) 
		  {
	          if ((ichan[ch].state == track) && ichan[ch].tow_sync) 
			  {
	    /* decode the navigation message for this channel */
	             navmess (ichan[ch].prn, ch);
			  }
		  }
		  if (out_data)
		  {
	          for (ch = 0; ch < N_channels; ch++) 
			  {
	             wrdl = ichan[ch].prn;
				 fputc(wrdl,data_bits);
			  }
			  for (i=0;i<1500;i++)
			  {
				  wrd= data_message[i];
				  wrdh= (wrd & 0xff00) >> 8;
				  wrdl= wrd & 0xff;
				  fputc(wrdh,data_bits);
				  fputc(wrdl,data_bits);
			  }
			  fprintf(data_bits,"\n");
		  }
	      data_frame_ready = 0;   
       }

       if (sec_flag == 1)
       {
	 /*          SendNMEA ();  */
	      nav_fix();  
		  /* Do the main routine stuff (display, position fix etc.) */
          thetime++;
          clock_tow = (clock_tow+1) % 604800L;
          time_on++;
          sec_flag = 0;
          for (ch = 0; ch < N_channels; ch++)
          {
              if (ichan[ch].state == track)
              {
                  if (ichan[ch].CNo < 33)
                  {
                      /* calculate carrier clock and doppler correction */
                        long carrier_corr = (-xyz[ichan[ch].prn].doppler -
                                               clock_offset * 1575.42) /
                        42.57475e-3;
                      /* calculate code clock and doppler correction */
                        long code_corr =
                        clock_offset * 24. + xyz[ichan[ch].prn].doppler / 65.5;

		                setup_channel (ch, ichan[ch].prn, code_corr, -carrier_corr);
				  } 
			  } 
		  } 
	   } 


   }
	  /* Update the Almanac Data file */
   write_almanac ();

  /* Update the Ephemeris Data file */
   write_ephemeris ();

  /* Update the ionospheric model and UTC parameters */
   write_ion_utc ();

  /* Update the curloc file for the next run */
   if (status == navigating)
    {
      out = fopen ("curloc.dat", "w+");
      fprintf (out, "latitude  %f\n", rec_pos_llh.lat * r_to_d);
      fprintf (out, "longitude %f\n", rec_pos_llh.lon * r_to_d);
      fprintf (out, "hae       %f\n", rec_pos_llh.hae);
      fclose (out);
    }
   return 0;
}

void generate_prn_codes(void)
{
	int i,j,G1,G2,prn,chip;
	int G2_i[33] = {0x000, 0x3f6, 0x3ec, 0x3d8, 0x3b0, 0x04b, 0x096, 0x2cb, 0x196,
		                   0x32c, 0x3ba, 0x374, 0x1d0, 0x3a0, 0x340, 0x280, 0x100,
					       0x113, 0x226, 0x04c, 0x098, 0x130, 0x260, 0x267, 0x338,
					       0x270, 0x0e0, 0x1c0, 0x380, 0x22b, 0x056, 0x0ac, 0x158};


	for (prn=1;prn<33;prn++)
	{
      G1 = 0x1FF;
      G2 = G2_i[prn];
	  prn_code[prn][0]=1;
	  for (chip=1;chip<1023;chip++)
	  {
		  prn_code[prn][chip]=(G1^G2) & 0x1;  /* exor the right hand most bit  */
		  i=((G1<<2)^(G1<<9)) & 0x200;
		  G1=(G1>>1) | i;
		  j=((G2<<1)^(G2<<2)^(G2<<5)^(G2<<7)^(G2<<8)^(G2<<9)) & 0x200;
		  G2=(G2>>1) | j;
	  }
	}
}


void Sim_GP2021(char IF[],long nsamp, long *tic)
{
  int ch,i,reg,prompt_chip,track_chip,tic_count;
  int Accum_status_A;
  double carrier_f,code_f;
  float  carrier_ph,track_ph,prompt_ph;
  float sin_carrier_ph,cos_carrier_ph;
  double t;

  tic_count=-1;
  for (i=0;i<nsamp;i++)
    {
      if (*tic==0) /* save the carrier and code info for measurements */
	{ 
          *tic=tic_ref;
          tic_count = i;
	} 
      else *tic=*tic-1;
    }


  Accum_status_A = 0;

  for (ch=0;ch<N_channels;ch++)
    {
      char *this_prn_code = prn_code[ichan[ch].prn];

      /* implement epoch set */
      if (REG_write[ch*8+7] != -1) 
	{
	  REG_read[ch*8+7] = REG_write[ch*8+7];
	  ms_counter[ch]=REG_write[ch*8+7] & 0xff;
	  bit_counter[ch]=REG_write[ch*8+7]>>8;
	}
      REG_write[ch*8+7]=-1;

      if (ichan[ch].prn>0)
	{
	  reg=ch<<3;
	  carrier_f = (REG_write[reg+3]*65536+REG_write[reg+4])*0.04257475; /* read the DCO frequencies */
	  code_f    = (REG_write[reg+5]*65536+REG_write[reg+6])*0.04257475;
	  for (i=0;i<nsamp;i++)
	    {

	      t=i/samp_rate;
	      carrier_ph = 2.0*M_PI*(carrier_f*t+carrier_phase[ch]);
	      sin_carrier_ph=sin(carrier_ph);
	      cos_carrier_ph=cos(carrier_ph);

	      if ((old_sin_carrier_phase[ch]<=0 && sin_carrier_ph>0) | (old_sin_carrier_phase[ch]<0 && 
									sin_carrier_ph>=0)) carrier_cycle[ch]++;
	      old_sin_carrier_phase[ch]=sin_carrier_ph;

	      prompt_ph   = t*code_f+code_phase[ch];
	      prompt_chip = (int)prompt_ph%1023;

	      i_prompt[ch] += IF[i]*sin_carrier_ph*(this_prn_code[prompt_chip]*2.0-1.0);
	      q_prompt[ch] += IF[i]*cos_carrier_ph*(this_prn_code[prompt_chip]*2.0-1.0);

	      track_ph    = prompt_ph-0.5;
	      track_chip = (int)track_ph%1023;

	      i_track[ch] += IF[i]*sin_carrier_ph*(this_prn_code[track_chip]*2.0-1.0);
	      q_track[ch] += IF[i]*cos_carrier_ph*(this_prn_code[track_chip]*2.0-1.0);

	      reg=(ch<<2)+0x84;
	      if  (prompt_ph>=(1023.0+((float)REG_write[reg])/2.0)) /*  dump  */
		{
		  if (ch==2) fprintf(out," %d  %d  %d  %d  %15.5f  %15.5f\n",REG_read[reg],REG_read[reg+1],
				     REG_read[reg+2],REG_read[reg+3],carrier_f,code_f);

		  code_phase[ch]=-i/samp_rate*code_f+(prompt_ph-(1023.0+((float)REG_write[reg])/2.0));
            
		  REG_read[reg  ] = i_track[ch]/1.4;      /* dump correlators into the appropriate 2021 registers */
		  REG_read[reg+1] = q_track[ch]/1.4;      /* and scale them to match the 2021 levels */
		  REG_read[reg+2] = i_prompt[ch]/1.4;
		  REG_read[reg+3] = q_prompt[ch]/1.4;

		  REG_write[reg] = 0;                     /* reset slew to 0  */
		  i_track[ch]  = q_track[ch]  = 0.0;      /* reset the correlators */
		  i_prompt[ch] = q_prompt[ch] = 0.0;

		  Accum_status_A = Accum_status_A | (1<<ch);     /* set the bit if a dump occurs */
            
		  dump[ch]=dump_ref;

		  /*  increment ms and bit counters  */

		  ms_counter[ch] ++;                                                 /* ms counter  */
		  if (ms_counter[ch]==20) bit_counter[ch] = (++bit_counter[ch])%50;  /* bit counter */
		  ms_counter[ch] = ms_counter[ch]%20;
		  REG_read[ch*8+7] = ms_counter[ch]+(bit_counter[ch]<<8);
		}
	      else dump[ch]--;

	      if (ch==0) tic_counter++;
	      if (i==tic_count) /* at TIC save the carrier and code info for measurements */
		{ 
		  reg=ch<<3;
		  REG_read[reg+4] = REG_read[reg+7];
		  REG_read[reg+3] = fmod(carrier_ph/2.0/M_PI,1.0)*1024;
		  REG_read[reg+1] = (int)(prompt_ph*2.0)%2046;                  /* number of half chips  */
		  REG_read[reg+5] = (prompt_ph*2.0-REG_read[reg+1])*1024;       /* half chip phase */

		  REG_read[reg+2] = carrier_cycle[ch] & 0xffff; /* carrier cycle low  */
		  REG_read[reg+6] = carrier_cycle[ch] >> 16;    /* carrier cycle high */

		  carrier_cycle[ch]=0;
		  if (ch==0) tic_counter=0;
		} 
	    }


	  /* set the carrier and code phases */
	  carrier_phase[ch] = fmod(nsamp*carrier_f/samp_rate+carrier_phase[ch],1.0); /* mod 2pi   */

	  code_phase[ch] = 	fmod(nsamp/samp_rate*code_f+code_phase[ch],2046.0);      /* mod 2046  */
	}
    }
  REG_read[0x82]=Accum_status_A;
  if (tic_count > -1) REG_read[0x83]=0x2000; /* if a tic occurs set the accum status B bit  */
  else REG_read[0x83]=0x0;
}


/*
int
ch_carrier_DCO_phase (int ch)
{
  return (from_gps ((ch << 3) + 0x3));
}

long
ch_carrier_cycle (int ch)
{
  long result;
  result = from_gps ((ch << 3) + 6);
  result = result << 16;
  result = result + from_gps ((ch << 3) + 2);
  return (result);
}

int
ch_code_phase (int ch)
{
  return (from_gps ((ch << 3) + 0x1));
}

unsigned int
ch_epoch (int ch)
{
  return (from_gps ((ch << 3) + 4));
}

int
ch_code_DCO_phase (int ch)
{
  return (from_gps ((ch << 3) + 5));
}

#ifdef UNUSED
void
data_tst (int data)
{
  to_gps (0xf2, data);
}

void
ch_code_incr_hi (int ch, int data)
{
  to_gps ((ch << 3) + 0x5, data);
}

void
ch_code_incr_lo (int ch, int data)
{
  to_gps ((ch << 3) + 0x6, data);
}

void
carr_incr_hi (int ch, int data)
{
  to_gps ((ch << 3) + 0x3, data);
}

void
carr_incr_lo (char ch, int data)
{
  to_gps ((ch << 3) + 0x4, data);
}


int
meas_status (void)
{
  return (from_gps (0x81));
}


void
ch_off (int ch)
{
  ch_status = ch_status & !(0x1 << (ch+1));
  reset_cntl (ch_status);
}


#endif 



correlator correl(float IF[],int nsamp, double samp_rate,double carrier_f,double carrier_phase,double code_f,float code_phase,int prn)
{
   int i,chip;
   float I,Q;
   double t;
   correlator output;

   I=Q=0.0;
   for (i=0;i<nsamp;i++)
   {
      t=i/samp_rate;

	  chip=(int(t*code_f+code_phase))%1023;
      I+=IF[i]*(sin(2.0*3.1415926*(carrier_f*t+carrier_phase))*(prn_code[prn][chip]*2.0-1.0));
      Q+=IF[i]*(cos(2.0*3.1415926*(carrier_f*t+carrier_phase))*(prn_code[prn][chip]*2.0-1.0));
   }
   output.I=I;
   output.Q=Q;
   return(output);
}
*/

/******************************************************************************
FUNCTION cold_allocate()
RETURNS  None.

PARAMETERS None.

PURPOSE To allocate the PRNs to channels for a cold start, start by
	  searching for PRN 1 through 12 and cycling through all PRN
	  numbers skipping channels that are tracking

WRITTEN BY
	Clifford Kelley

******************************************************************************/
void
cold_allocate ()
{
  int ch, i, alloc;
  /* d_freq=4698; */           /* search 200 Hz intervals */
  search_max_f = 50;            /* widen the search for a cold start */
  satfind (0);
  /* almanac_valid=1; */
  for (i = 1; i <= 32; i++)
    {
      if (gps_alm[i].inc > 0.0 && gps_alm[i].week != gps_week % 1024)
        almanac_valid = 0;
    }
  if (al0 == 0.0 && b0 == 0.0)
    almanac_valid = 0;
  for (ch = 0; ch < N_channels; ch++)   /* if no satellite is being tracked */
    /* turn the channel off */
    {
      if (ichan[ch].CNo < 30)  /* if C/No is too low turn the channel off */
	channel_off (ch);
    }
  for (i = 0; i <= chmax; i++)
    {
      alloc = 0;
      for (ch = 0; ch < N_channels; ch++)
        {
          if (ichan[ch].prn == cold_prn)
            {
              alloc = 1;        /* satellite is already allocated a channel */
	      cold_prn = cold_prn % 32 + 1;
              break;
            }
        }
      if (alloc == 0)           /* if not allocated find an empty channel */
        {
          for (ch = 0; ch < N_channels; ch++)
            {
              if (ichan[ch].state == off)
                {
                  long carrier_corr = -clock_offset * 1575.42 / 42.57475e-3;
		  setup_channel (ch, cold_prn, 0, -carrier_corr); /* opposite sign of GP2021 */
                  cold_prn = cold_prn % 32 + 1;
                  break;
                }
            }
        }
    }
}

/******************************************************************************
FUNCTION read_rcvr_par(void)
RETURNS  None.

PARAMETERS None.

PURPOSE   To read in from the rcvr_par file the receiver parameters that
			control acquisition, tracking etc.

WRITTEN BY
	Clifford Kelley

******************************************************************************/
void
read_rcvr_par (void)
{
  char intext[40];
  if ((in = fopen ("rcvr_par.dat", "r")) == NULL)
    {
      printf ("Cannot open rcvr_par.dat file.\n");
      exit (0);
    }
  else
    {
      fscanf (in, "%s %s", intext, tzstr);
      fscanf (in, "%s %f", intext, &mask_angle);
      mask_angle = mask_angle / r_to_d;
      fscanf (in, "%s %f", intext, &clock_offset);
      fscanf (in, "%s %u", intext, &interr_int);
      fscanf (in, "%s %d", intext, &cold_prn);
      fscanf (in, "%s %d", intext, &ICP_CTL);
      fscanf (in, "%s %f", intext, &nav_up);
      fscanf (in, "%s %d", intext, &out_pos);
      fscanf (in, "%s %d", intext, &out_vel);
      fscanf (in, "%s %d", intext, &out_time);
      fscanf (in, "%s %d", intext, &out_kalman);
      fscanf (in, "%s %d", intext, &out_debug);
      fscanf (in, "%s %d", intext, &out_data);
      fscanf (in, "%s %d", intext, &m_tropo);
      fscanf (in, "%s %d", intext, &m_iono);
      fscanf (in, "%s %d", intext, &align_t);
      fscanf (in, "%s %u", intext, &Com0Baud);  /* NMEA */
      fscanf (in, "%s %u", intext, &Com1Baud);  /* NMEA */
      fscanf (in, "%s %u", intext, &GPGGA);     /* NMEA */
      fscanf (in, "%s %u", intext, &GPGSV);     /* NMEA */
      fscanf (in, "%s %u", intext, &GPGSA);     /* NMEA */
      fscanf (in, "%s %u", intext, &GPVTG);     /* NMEA */
      fscanf (in, "%s %u", intext, &GPRMC);     /* NMEA */
      fscanf (in, "%s %u", intext, &GPZDA);     /* NMEA */
	  fscanf (in, "%s %d", intext, &out_rinex); /* RINEX data logging */
    }
  fclose (in);
}

/******************************************************************************
FUNCTION read_rinex_par(void)
RETURNS  None.

PARAMETERS None.

PURPOSE   Reads in the RINEX file parameters from the rinexpar.dat file

WRITTEN BY
	Jonathan J. Makela

******************************************************************************/
void
read_rinex_par (void)
{
  char intext[40];
  char *p;

  if ((in = fopen (rinex_par_file, "r")) == NULL)
    {
      printf ("Cannot open %s file.\n",rinex_par_file);
      exit (0);
    }
  else
    {
      fscanf (in, "%s %s\n", intext, system_type);
      fscanf (in, "%s %s\n", intext, program_name);
      fscanf (in, "%s %s\n", intext, agency_name);
      fscanf (in, "%s %s\n", intext, marker_name);
      fscanf (in, "%s\t", intext);
      fgets(observer_name, 20, in);
      /* Remove newline character if it is there */
      p = strchr(observer_name, '\n');
      if(p != NULL)
	*p = '\0';
      else
	fscanf (in, "%[ a-z A-Z 0-9]\n", intext);
      fscanf (in, "%s %s\n", intext, receiver_number);
      fscanf (in, "%s %s\n", intext, receiver_type);
      fscanf (in, "%s %s\n", intext, receiver_version);
      fscanf (in, "%s %s\n", intext, antenna_number);
      fscanf (in, "%s %s\n", intext, antenna_type);
      fscanf (in, "%s %lf\n", intext, &loc_x);
      fscanf (in, "%s %lf\n", intext, &loc_y);
      fscanf (in, "%s %lf\n", intext, &loc_z);
      fscanf (in, "%s %lf\n", intext, &delx);
      fscanf (in, "%s %lf\n", intext, &dely);
      fscanf (in, "%s %lf\n", intext, &delz);
      fscanf (in, "%s %d\n", intext, &lamda_factor_L1);
      fscanf (in, "%s %d\n", intext, &lamda_factor_L2);
      fscanf (in, "%s %d\n", intext, &n_obs);
      fscanf (in, "%s %s\n", intext, obs1);
      fscanf (in, "%s %s\n", intext, obs2);
      fscanf (in, "%s %s\n", intext, time_system);
    }
  fclose (in);
}


/******************************************************************************
FUNCTION display()
RETURNS  None.

PARAMETERS None.

PURPOSE
	This function displays the current status of the receiver on the
	computer screen.  It is called when there is nothing else to do

WRITTEN BY
	Clifford Kelley

******************************************************************************/
void
display (void)
{
  int ch;
  time_t utctime = thetime - dtls;
  /*  goto_screen_location(1, 1);  */

  printf ("                   OpenSource GPS Software Version 2.00\n");
  printf ("%s", ctime (&utctime));
  printf ("TOW  %6ld\n", clock_tow);
  printf ("meas time %f  error %f  delta %f\n", m_time[1], m_error,
          delta_m_time);
  cur_lat.deg = rec_pos_llh.lat * r_to_d;
  cur_lat.min = (rec_pos_llh.lat * r_to_d - cur_lat.deg) * 60;
  cur_lat.sec =
    (rec_pos_llh.lat * r_to_d - cur_lat.deg - cur_lat.min / 60.) * 3600.;
  cur_long.deg = rec_pos_llh.lon * r_to_d;
  cur_long.min = (rec_pos_llh.lon * r_to_d - cur_long.deg) * 60;
  cur_long.sec =
    (rec_pos_llh.lon * r_to_d - cur_long.deg - cur_long.min / 60.) * 3600.;
  printf ("   latitude    longitude          HAE      clock error (ppm)\n");
  printf ("  %4d:%2d:%5.2f  %4d:%2d:%5.2f  %10.2f  %f\n",
          cur_lat.deg, abs (cur_lat.min), fabs (cur_lat.sec), cur_long.deg,
          abs (cur_long.min), fabs (cur_long.sec), rec_pos_llh.hae,
          clock_offset);
  printf (" Speed     Heading      TIC_dt\n");
  printf (" %f   %f   %f\n", speed, heading * r_to_d, TIC_dt);
  printf ("   \n");
  printf
    ("tracking %2d status %1d almanac valid %1d gps week %4d alm_page %2d\n",
     n_track, status, almanac_valid, gps_week % 1024, alm_page);
  if (display_page == 0)
    {
      printf
        (" ch prn state n_freq az el doppler t_count n_frame sfid ura page missed CNo\n");
      for (ch = 0; ch < N_channels; ch++)
        {
          printf
            (" %2d %2d  %2d  %3d   %4.0f  %3.0f   %6.0f   %4d  %4d  %2d  %3d  %3d%5d     %2d\n",
             ch, ichan[ch].prn, ichan[ch].state, ichan[ch].n_freq,
             xyz[ichan[ch].prn].azimuth * 57.3,
             xyz[ichan[ch].prn].elevation * 57.3, xyz[ichan[ch].prn].doppler,
             ichan[ch].frame_bit % 1500, ichan[ch].n_frame, ichan[ch].sfid,
             gps_eph[ichan[ch].prn].ura, schan[ch].page5, ichan[ch].missed,
             ichan[ch].CNo);
        }
      printf (" GDOP=%6.3f  HDOP=%6.3f  VDOP=%6.3f  TDOP=%6.3f\n", gdop, hdop,
	      vdop, tdop);
    }
  else if (display_page == 1)
    {
      printf (" ch prn state TLM      TOW  Health  Valid  TOW_sync offset\n");
      for (ch = 0; ch < N_channels; ch++)
        {
          printf (" %2d %2d  %2d  %6ld   %6ld   %2d     %2d     %2d   %4d\n",
                  ch, ichan[ch].prn, ichan[ch].state, ichan[ch].TLM,
                  ichan[ch].TOW, gps_eph[ichan[ch].prn].health,
                  gps_eph[ichan[ch].prn].valid, ichan[ch].tow_sync, 0);
        }
    }
  else if (display_page == 2)
    {
      printf (" ch prn state n_freq az  el        tropo        iono\n");
      for (ch = 0; ch < N_channels; ch++)
        {
          printf (" %2d %2d  %2d  %3d   %4.0f  %3.0f   %10.4f   %10.4f\n",
                  ch, ichan[ch].prn, ichan[ch].state, ichan[ch].n_freq,
                  xyz[ichan[ch].prn].azimuth * 57.3,
                  xyz[ichan[ch].prn].elevation * 57.3, 
		  satellite[ichan[ch].prn].Tropo * c,
                  satellite[ichan[ch].prn].Iono * c);
        }
    }
  else if (display_page == 3)
    {
      printf (" ch prn state      Pseudorange     delta Pseudorange\n");
      for (ch = 0; ch < N_channels; ch++)
        {
          printf (" %2d %2d  %2d  %20.10f   %15.10f\n",
                  ch, ichan[ch].prn, ichan[ch].state, 
		  satellite[ichan[ch].prn].Pr,
                  satellite[ichan[ch].prn].dPr);
        }
    }
  else if (display_page == 4)   /* can be used for debugging purposes */
    {
      printf (" ch prn state sfid page SF1  SF2  SF3  SF4  SF5\n");
      for (ch = 0; ch < N_channels; ch++)
        {
          printf (" %2d %2d   %2d   %2d  %2d  %3x  %3x  %3x  %3x  %3x\n",
                  ch, ichan[ch].prn, ichan[ch].state, ichan[ch].sfid,
                  schan[ch].page5, schan[ch].word_error[0],
                  schan[ch].word_error[1], schan[ch].word_error[2],
                  schan[ch].word_error[3], schan[ch].word_error[4]);
        }
    }

}
/******************************************************************************
FUNCTION chan_allocate()
RETURNS  None.

PARAMETERS None.

PURPOSE
	This function allocates the channels with PRN numbers

WRITTEN BY
	Clifford Kelley

******************************************************************************/

void
chan_allocate ()
{
  int ch, prnn, alloc;
  int i;
  /*       almanac_valid=1; */
  for (prnn = 1; prnn <= 32; prnn++)
    {
      xyz[prnn] = satfind (prnn);
      if (gps_alm[prnn].inc > 0.0 && gps_alm[prnn].week != gps_week % 1024)
        {
          almanac_valid = 0;
        }
    }
  if (al0 == 0.0 && b0 == 0.0)
    almanac_valid = 0;
  for (ch = 0; ch < N_channels; ch++)
    /* if the sat has dropped below mask angle turn the channel off */
    {
      if (ichan[ch].CNo < 30)
        {
          memset (&satellite[ichan[ch].prn], 0, sizeof(struct satellite));
	       channel_off (ch);
          for (i = 0; i < 5; i++)
            schan[ch].word_error[i] = 0;
        }
    }
  for (prnn = 1; prnn <= 32; prnn++)
    {
      if (xyz[prnn].elevation > mask_angle && gps_alm[prnn].health == 0 &&
          gps_alm[prnn].ety != 0.00)
        {
          alloc = 0;
          for (ch = 0; ch < N_channels; ch++)
            {
              if (ichan[ch].prn == prnn)
                {
                  alloc = 1;    /* satellite already allocated a channel */
                  break;
                }
            }
          if (alloc == 0)       /* if not allocated find an empty channel */
            {
              for (ch = 0; ch < N_channels; ch++)
                {
                  if (ichan[ch].state == off)
                    {
                      /* calculate carrier clock and doppler correction */
                      long carrier_corr = (-xyz[prnn].doppler -
					   clock_offset * 1575.42) /
                        42.57475e-3;
                      /* calculate code clock and doppler correction */
                      long code_corr =
                        clock_offset * 24. + xyz[prnn].doppler / 65.5;

		      setup_channel (ch, prnn, code_corr, -carrier_corr);  /* opposite sign of GP2021  */
                      break;
                    }
                }
            }
        }
    }
}

int
pcifind (void)
{
  return 0;
}

void read_filenames(void)
{
  char txt1[40],txt2[40];
  FILE *filenames;
  filenames = fopen("SW_files.def","r+");
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(almanac_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(ephemeris_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(receiver_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(location_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(ion_utc_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(output_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(kalman_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(data_bits_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(rinex_par_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(rinex_obs_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(rinex_nav_file,"%s",txt2);
  fscanf(filenames,"%s %s",txt1,txt2);
  sprintf(debug_file,"%s",txt2);

  fclose(filenames);
}
