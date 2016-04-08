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

struct dic// Dictionary
{
   int *index;// Index of combination
   int *freq;// Frequency of repetition
   int *id_pos;// Index on the pos vector
   int *pos;// Where index appears
};

struct read// Used to read sequences
{
   short *data;
   int *length;
   int *start;
   struct read *next;
};

struct tmp_data// Used to auxiliate the struct read
{
   short *data;
   uint length;
   struct tmp_data *next;
};

#endif
