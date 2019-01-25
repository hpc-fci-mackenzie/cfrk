#ifndef _tipos_h
#define _tipos_h

//Macro performing 4^n
#define POW(k) (1U << 2*(k))

//Typedef
typedef unsigned short ushort;
typedef long int lint;
typedef unsigned int uint;

//0 disabled; 1 enable
const int DBG = 0;

struct seq
{
   char *header;
   char *seq;
   int len;
   struct seq *next;
};

struct read// Used to read sequences
{
   char *data;
   int *length;
   lint *start;
   int *Freq;
   struct read *next;
};

struct tmp_data// Used to auxiliate the struct read
{
   char *data;
   int length;
   struct tmp_data *next;
};

#endif
