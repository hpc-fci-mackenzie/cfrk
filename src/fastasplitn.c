/* FASTASPLITN.C
 * Version:   1.0.7  12-MAY-2015
 * Author:    David Mathog, Biology Division, Caltech
 * email:     mathog@caltech.edu
 * Copyright: 2015 David Mathog and California Institute of Technology (Caltech)
 *
 * Description:
 *    Splits a fasta file into a N files, with cycling entries among the N
 *       output  files. This does a better job of randomizing the split content
 *          than did FASTASPLIT.C, which doesn't work well when there is a systematic
 *             variation in the input file.  For instance, in the NCBI "nt" database the
 *                entries get  longer and longer from the beginning to the end, so that
 *                   regular FASTASPLIT  produced a final file 5 times bigger than the first one
 *                      (all with the same  number of sequences though!)
 *
 *                         If the optional P value is specified, it emits output to
 *                            stdout with contents identical to the Pth fragment of a full
 *                               N file set. 
 *
 *                                  If in addition the optional C value is specified the number of sequential
 *                                     entries emitted to each stream may be changed from the default of 1.
 *
 *                                        Very little error checking.
 *                                           
 *                                              No input line may exceed 1M  characters.
 *                                                 
 *                                                    Compiles cleanly with:
 *                                                       
 *                                                          gcc -Wall -std=c99 -pedantic -o fastasplitn fastasplitn.c
 *
 *                                                          Arguments are:
 *
 *                                                           1:  infile (- or STDIN mean read from stdin)
 *                                                            2:  N number of files to produce. 
 *                                                             3:  P phase. If specified and in range 1-N only a single
 *                                                                    file  is produced as if it was the Pth of the N files. 
 *                                                                           If P=0 then N files are produced.  Anything else is an error.
 *                                                                            4:  C cycle.  If specified must be >=1.  If C=3, N=2 then
 *                                                                                  1,2,3,7,8,9, etc.    -> file1
 *                                                                                        4,5,6,10,11,12, etc. -> file2
 *                                                                                              C=1 is the default.
 *                                                                                                    
 *                                                                                                     
 *                                                                                                     License terms:
 *                                                                                                         You may run this program on any platform. You may
 *                                                                                                             redistribute the source code of this program subject to
 *                                                                                                                 the condition that you do not first modify it in any way.
 *                                                                                                                     You may  distribute binary versions of this program so long
 *                                                                                                                         as they were compiled from unmodified source code.  There
 *                                                                                                                             is no charge for using this software.  You may not charge
 *                                                                                                                                 others for the use of this software.
 *
 *                                                                                                                                 Bug reports:  Please report bugs via email.
 *
 *                                                                                                                                 1.0.7  12-MAY-2015.  Fixed bug on default stdin (no -in)
 *                                                                                                                                 1.0.6  06-MAY-2015.  Added -2, -q2f, -i, -h, -bs, --help, -in, -t options.
 *                                                                                                                                        Changed input format so that it can use -p -n -c parameters.
 *                                                                                                                                               Old style parameters are supported but deprecated.  The use of any of the new
 *                                                                                                                                                      style ones means all must be new style.
 *                                                                                                                                                      1.0.5  13-JUN-2013.  Increased MAXSTRING size.
 *                                                                                                                                                      1.0.4  Slightly clarified command line usage, after 
 *                                                                                                                                                             a suggestion by Bernd Brandt.
 *                                                                                                                                                             1.0.3  Added code to implement C value.  Useful for splitting files
 *                                                                                                                                                                    into sequential sequences if one knows ahead of time how many
 *                                                                                                                                                                           sequences are present in the input file.  Ie, if there are 20
 *                                                                                                                                                                                  and N=4 C=5 P=0 then 4 files will be produced with the first
 *                                                                                                                                                                                         holding seqs 1->5, the second 6->10, and so forth.
 *                                                                                                                                                                                         1.0.2  Added code to read symbol "SPLITFRAGTEMPLATE" which replaces
 *                                                                                                                                                                                                default value template for output filename creation.  Cleaned
 *                                                                                                                                                                                                       up printf so that all messages go to stderr.
 *                                                                                                                                                                                                       */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>   /* for toupper */
#include <limits.h>  /* for INT_MAX */

/* definitions */
#define EXVERSTRING "1.0.7 12-MAY-2015"
#define COPYSTRING  "2015 David Mathog and California Institute of Technology"
#define BUGSTRING   "mathog@caltech.edu"
#define LICSTRING   "GNU General Public License 2"

#define MYMAXSTRING 2000000
#define MAXOUTFILES 1000  /* WAY more than should ever be needed.  ulimit may restrict to far less than this. */

/* Prototypes */
void emit_help(void);
void emit_hexamples(void);
void insane(char *string);
int  lcl_strcasecmp(const char *s1, const char *s2);
void process_command_line_args(int argc,char **argv);
void setnonnegnumeric(int *val,int *numarg,int argc,char **argv,char * label);
void process_command_line_args(int argc,char **argv);

/* globals */
int   gbl_streams; /* n */
int   gbl_phase;   /* p */
int   gbl_cyc_len;   /* c */
int   gbl_f2, gbl_q4, gbl_q2f, gbl_bs;
const char static_template[]="frag%3.3d";
const char *gbl_template = &static_template[0];
const char *gbl_infile   = NULL; /* stdin */

int  lcl_strcasecmp(const char *s1, const char *s2){
int c1;
int c2;
  for(; ;s1++,s2++){
    c1=toupper(*s1);
    c2=toupper(*s2);
    if(c1 < c2)return -1;
    if(c1 > c2)return  1;
    if(c1 == 0)return  0;  /*c2 also is 0 in this case */
  }
}

void emit_help(void){
(void) fprintf(stderr,"fastasplitn [options]             or \n");
(void) fprintf(stderr,"fastasplitn infile N [P [C]]     (deprecated form, fasta only)\n");
(void) fprintf(stderr,"   Distribute entries from a fasta or fastq file into N output streams. \n\n");
(void) fprintf(stderr,"Command line options:\n");
(void) fprintf(stderr,"  -in INFILE  Input file name.  \"-\" or \"stdin\" = read from stdin [Default: from stdin]\n");
(void) fprintf(stderr,"  -n N  Number of data streams to produce [No default: N>=1, <=%d",MAXOUTFILES);
(void) fprintf(stderr,"  -p P  Phase to emit.  If 0 [Default] emit N streams to N files\n");
(void) fprintf(stderr,"          If 1->N emit only contents of Pth stream to stdout.\n");
(void) fprintf(stderr,"  -c C  Transfer a block of C consecutive entries from input to current output stream.\n");
(void) fprintf(stderr,"          Then shift to the next output stream. C>=1, [Default: 1]\n");
(void) fprintf(stderr,"  -t TEMPLATE  Output stream file name template. [Default: frag%%3.3d]\n");
(void) fprintf(stderr,"           The contents of environmental symbol SPLITFRAGTEMPLATE overrides this setting\n");
(void) fprintf(stderr,"  -f2   Input fasta file has only one line of sequence per entry [Default: sequence may be multiline]\n");
(void) fprintf(stderr,"  -q4   Input is fastq file, each entry has exactly 4 lines. [Default: input is fasta]\n");
(void) fprintf(stderr,"  -q2f  Convert fastq to fasta. [Default: pass entries as is]\n");
(void) fprintf(stderr,"  -bs BS   Set the buffer size, all input lines must be shorter than this.\n");
(void) fprintf(stderr,"        BS cannot be more than %d.  [Default: %d]\n",INT_MAX, MYMAXSTRING);
(void) fprintf(stderr,"  -h    Print this help message (also -help --h --help -? --?)\n");
(void) fprintf(stderr,"  -hexamples  Print examples and notes\n");
(void) fprintf(stderr,"  -i    Emit version, copyright, license and contact information\n\n");
}

void emit_hexamples(void){
(void) fprintf(stderr,"Examples: in all examples the input has 20 sequences and is read from stdin\n");
(void) fprintf(stderr,"  fastasplitn -n 4 -p 0 -c 5 -t 'out%%1.1d'\n");
(void) fprintf(stderr,"     Entries  1->5 to file out1, 6->10 to out2, 11->15 to out3, 16->20 to out4\n");
(void) fprintf(stderr,"  fastasplitn -n 4 -p 2 -c 3:\n");
(void) fprintf(stderr,"     Defines 4 streams but only the 2nd one is used, blocks of 3 entries\n");
(void) fprintf(stderr,"     at a time are processed and these are sent to stdout.\n");
(void) fprintf(stderr,"     Entries 4->6 and 16->17 are written to stdout\n");
(void) fprintf(stderr,"  fastasplitn -q4 -q2f -n 4 -p 2 -c 3:\n");
(void) fprintf(stderr,"     Fastq entries 4->6 and 16->17 are read as fastq, converted to fasta,\n");
(void) fprintf(stderr,"     and written to stdout\n\n");
(void) fprintf(stderr,"Notes:\n");
(void) fprintf(stderr,"  Parallel processing of the output streams may be implemented\n");
(void) fprintf(stderr,"    on linux or unix systems as shown below. Take care with the FIFO names\n");
(void) fprintf(stderr,"    to avoid collisions with other processes, and be sure that all\n");
(void) fprintf(stderr,"    processes have exited before removing the FIFO files.\n\n");
(void) fprintf(stderr,"  mkfifo /tmp/stepA1 /tmp/stepA2 /tmp/stepA3 /tmp/stepB1 /tmp/stepB2 /tmp/stepB3\n");
(void) fprintf(stderr,"  nohup a_merge_program /tmp/stepB1  /tmp/stepB2 /tmp/stepB3 >final_results.txt 2>/dev/null &\n");
(void) fprintf(stderr,"  nohup a_processing_program </tmp/stepA1 >/tmp/stepB1 2>/dev/null &\n");
(void) fprintf(stderr,"  nohup a_processing_program </tmp/stepA2 >/tmp/stepB2 2>/dev/null &\n");
(void) fprintf(stderr,"  nohup a_processing_program </tmp/stepA3 >/tmp/stepB3 2>/dev/null &\n");
(void) fprintf(stderr,"  fastasplitn -in input.fa -t '/tmp/stepA%%1.1d' -n 3\n\n");
}

void insane(char *string){
 (void) fprintf(stderr,"%s\n",string);
 exit(EXIT_FAILURE);
}

void setnonnegnumeric(int *val,int *numarg,int argc,char **argv,char * label){
      (*numarg)++;
      if( ( *numarg >= argc ) || (argv[*numarg] == NULL)){
        (void) fprintf( stderr, "%s: missing argument\n",label);
        exit(EXIT_FAILURE);
      }
      if(sscanf(argv[*numarg],"%d",val) != 1){
        (void) fprintf(stderr,"Bad integer argument/parameter [%s %s] \n",label,argv[*numarg]);
        exit(EXIT_FAILURE);
      }
      if(*val < 0){
        (void) fprintf(stderr,"Illegal negative integer argument/parameter [%s %s] \n",label,argv[*numarg]);
        exit(EXIT_FAILURE);
      }
}

void process_command_line_args(int argc,char **argv){

   int numarg=0;
   char *env_template=NULL;
 
   gbl_phase    = 0;
   gbl_cyc_len    = 1;
   gbl_streams  = 0;
   gbl_f2       = 0;
   gbl_q4       = 0;
   gbl_q2f      = 0;
   gbl_bs       = MYMAXSTRING;
   /* gbl_infile   set above */
   /* gbl_template set above */

   while( ++numarg < argc){
      if( (lcl_strcasecmp(argv[numarg], "-h")==0)     ||
          (lcl_strcasecmp(argv[numarg], "--h")==0)    ||
          (lcl_strcasecmp(argv[numarg], "-?")==0)     ||
          (lcl_strcasecmp(argv[numarg], "--?")==0)    ||
          (lcl_strcasecmp(argv[numarg], "-help")==0)  ||
          (lcl_strcasecmp(argv[numarg], "--help")==0) ){
        emit_help();
        exit(EXIT_SUCCESS);
      }
      else if(lcl_strcasecmp(argv[numarg], "-hexamples")==0){
        emit_hexamples();
        exit(EXIT_SUCCESS);
      }
      else if(lcl_strcasecmp(argv[numarg], "-in")==0){
         gbl_infile=argv[++numarg];
         continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-t")==0){
         gbl_template=argv[++numarg];
         continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-i")==0){
        (void)fprintf(stderr,"Version:   %s\n",EXVERSTRING);
        (void)fprintf(stderr,"bugs to:   %s\n",BUGSTRING);
        (void)fprintf(stderr,"Copyright: %s\n",COPYSTRING);
        (void)fprintf(stderr,"License:   %s\n",LICSTRING);
        exit(EXIT_SUCCESS);
      }
      else if(lcl_strcasecmp(argv[numarg], "-bs")==0){
        setnonnegnumeric(&gbl_bs,&numarg,argc,argv,"-bs");
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-n")==0){
        setnonnegnumeric(&gbl_streams,&numarg,argc,argv,"-n");
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-p")==0){
        setnonnegnumeric(&gbl_phase,&numarg,argc,argv,"-p");
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-c")==0){
        setnonnegnumeric(&gbl_cyc_len,&numarg,argc,argv,"-c");
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-f2")==0){
        gbl_f2=1;
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-q4")==0){
        gbl_q4=1;
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-q2f")==0){
        gbl_q2f=1;
        continue;
      }
      else if(lcl_strcasecmp(argv[numarg], "-")==0 || lcl_strcasecmp(argv[numarg], "stdin")==0){
         /* deprecated format */
         if( (argc < 3)                                   ||
             (argc > 6)                                   ||
             (argc >=3 && (sscanf(argv[2],"%d",&gbl_streams)==EOF)) ||
             gbl_streams < 1                                        ||
             (argc >=4 && 
               (     (sscanf(argv[3],"%d",&gbl_phase)==EOF)  ||
                     gbl_phase < 0                           ||
                     gbl_phase > gbl_streams
               )
             )                                            ||
             (argc >=5 && 
               (     (sscanf(argv[4],"%d",&gbl_cyc_len)==EOF)  ||
                     gbl_cyc_len < 1
               )
             )
         ){
            emit_help();
            exit(EXIT_FAILURE);
        }
        gbl_infile=argv[1];
        break;
      }
      else {
        (void) fprintf(stderr,"Unknown command line argument: %s\n",argv[numarg]);
        emit_help();
        exit(EXIT_FAILURE);
        continue;
      }
   } 

   /* template can alos be set by an ENV symbol (deprecated) */

   env_template=getenv("SPLITFRAGTEMPLATE");
   if(env_template){ gbl_template = env_template; }
   else {
      env_template=getenv("splitfragtemplate");
      if(env_template){ gbl_template = env_template; }
   }

   /* sanity checking */

   if(!gbl_streams){
      if(argc==1){
         emit_help();
         exit(EXIT_FAILURE);
      }
      else {
         insane("fastasplitn: fatal error: -n parameter either 0 or not set");
      }
   }
   if(gbl_streams > MAXOUTFILES){
      (void) fprintf(stderr,"fastasplitn: fatal error: -n parameter too large, greater than %d\n",MAXOUTFILES);
      exit(EXIT_FAILURE);
   }
   if((gbl_phase > gbl_streams))insane("fastasplitn: fatal error: -p parameter may not be larger than -n parameter");
   if((gbl_cyc_len > gbl_streams))insane("fastasplitn: fatal error: -c parameter may not be larger than -n parameter");
   if(gbl_f2 && gbl_q4)insane("fastasplitn: fatal error: -f2 and -q4 cannot be used together");
   if(gbl_q2f && !gbl_q4)insane("fastasplitn: fatal error: -q2f requires -q4");
   if(!strchr(gbl_template,'%'))insane("fastasplitn: fatal error: output file template contains no %% characters");

   /* emit run parameters to stderr */
   (void) fprintf(stderr,"fastasplitn: status: Reading from:    %s\n",(gbl_infile ? gbl_infile : "stdin"));
   (void) fprintf(stderr,"fastasplitn: status: Output template: %s\n",gbl_template);
   (void) fprintf(stderr,"fastasplitn: status: Streams:         %d\n",gbl_streams);
   (void) fprintf(stderr,"fastasplitn: status: Phases:          %d\n",gbl_phase);
   (void) fprintf(stderr,"fastasplitn: status: Cycle length:    %d\n",gbl_cyc_len);
   (void) fprintf(stderr,"fastasplitn: status: Buffer size:     %d\n",gbl_bs);
   (void) fprintf(stderr,"fastasplitn: status: Input type:      %s\n",(gbl_q4 ? "fastq" : "fasta"));
   (void) fprintf(stderr,"fastasplitn: status: Output type:     %s\n",(gbl_q4 & !gbl_q2f ? "fastq" : "fasta"));
}

int main(int argc, char *argv[]){
char *bigstring;
char outname[200];
int streamcount=0;
unsigned long long stm_count  = 0; /* index for outupt streams */
unsigned long long cyc_count  = 0; /* index for blocks of output entries */
unsigned long long ent_count  = 0; /* count of entries */
unsigned long long rec_count  = 0; /* count of input records */
unsigned long long remainder  = 0; /* number of records to emit for this entry */
FILE *fin=NULL;
FILE *fout[MAXOUTFILES]; 
  
   process_command_line_args(argc,argv);

   bigstring=malloc(gbl_bs * sizeof(char));
   if(!bigstring)insane("fastasplitn: fatal error: could not allocate memory for input buffer\n");
  
  
   if(!gbl_infile || lcl_strcasecmp(gbl_infile,"-") == 0 || lcl_strcasecmp(gbl_infile,"stdin") == 0){
      fin=stdin;
   }
   else {
     fin=fopen(gbl_infile,"r");
   }
   if(!fin)insane("fastasplitn: fatal error: could not open input file\n");

   /* see if template is to use default or from a symbol */
  

   /* open all the output files */

   for(streamcount=0; streamcount<gbl_streams ;streamcount++){
      (void) sprintf(outname,gbl_template,streamcount+1);
      if(gbl_phase==0){
         (void) fprintf(stderr,"fastasplitn: status: Opening output file %s\n",outname);
         fout[streamcount] = fopen(outname,"w");
         if(fout[streamcount]==NULL){
            (void) fprintf(stderr,"fastasplitn: fatal error: Could not open output file %s\n",outname);
            exit(EXIT_FAILURE);
         }
      }
      else if(gbl_phase==streamcount+1){
         fout[streamcount] = stdout;
      }
   }

   /* now send fasta entries to each one */

   stm_count = 0;
   cyc_count = 0;
   ent_count = 0;
   /* the method used to check for input line truncation is not guaranteed to work by 
 *       the C language standards.  However, it would only fail if a compiler or library
 *             did something really crazy in its fgets() implementation. */
   bigstring[gbl_bs - 1]='\1';
   while( fgets(bigstring,gbl_bs,fin) != NULL){
      if(bigstring[gbl_bs-1]=='\0' && bigstring[gbl_bs-2]!='\n'){
         (void) fprintf(stderr,"fastasplitn: fatal error: line %lld is longer than the buffer, increase -bs\n",rec_count);
         exit(EXIT_FAILURE);
      }
      if(!ent_count){
         if(gbl_q4){
            if(*bigstring != '@')insane("fastasplitn: fatal error: input is not in fastq format");
         }
         else if( *bigstring != '>' )insane("fastasplitn: fatal error: input is not fasta format");
      }
      if(  (gbl_f2 && !(rec_count & 1)) ||
           (gbl_q4 && !(rec_count & 3)) ||
           bigstring[0] == '>'){
         ent_count++;
         cyc_count++;
         if(cyc_count > gbl_cyc_len){
            cyc_count=1;
            stm_count++;
            if(stm_count >= gbl_streams)stm_count=0;
         }
         if(gbl_f2){
            remainder = 2;
         }
         else if(gbl_q4){ 
            if(gbl_q2f){ 
               bigstring[0]='>';
               remainder = 2;
            }
            else { remainder = 4; }
         }
         else { 
            remainder = LLONG_MAX;
         }
      }
      rec_count++;
      if(remainder){
         if(gbl_phase==0 || gbl_phase==stm_count+1)(void) fprintf(fout[stm_count],"%s",bigstring);
         remainder--;
      }
      bigstring[gbl_bs-1]='\1';
   }
   if(gbl_phase==0)(void) fprintf(stderr,"fastasplitn: status: Processing completed on number of entries:  %lld\n",ent_count);
   exit(EXIT_SUCCESS);
}
