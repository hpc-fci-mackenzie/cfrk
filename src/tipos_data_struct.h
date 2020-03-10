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
   char *read;
   char *data;
   int len;
};

struct chunk// Used to read sequences
{
   char *data;
   int *length;
   lint *start;
   int *n_combination;
   struct counter *counter;
   // struct read *reads;
   struct chunk *next;
};

struct read 
{
   
   
};

struct counter {
    int *index;
    int *frequence;
};

#endif
