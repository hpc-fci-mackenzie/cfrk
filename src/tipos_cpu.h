#ifndef _tipos_h
#define _tipos_h

//Macro performing 4^n
#define POW(k) (1U << 2*(k))

//Typedef
typedef unsigned short ushort;
typedef long int lint;
typedef unsigned int uint;

//0 disabled; 1 enable
// const int DBG = 0;

struct seq
{
   char *header;
   char *read;
   int *data;
   int len;
};

struct read// Used to read sequences
{
   int *data;
   int *length;
   lint *start;
   struct counter *counter;
   struct read *next;
};

struct counter {
    int **index;
    int *frequency;
};

#endif
