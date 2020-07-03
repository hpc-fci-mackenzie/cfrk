#ifndef _tipos_h
#define _tipos_h

//Macro performing 4^n
#define POW(k) (1U << 2*(k))

//Typedef
typedef unsigned short ushort;
typedef long int lint;
typedef unsigned int uint;

struct seq
{
    char *header;
    char *read;
    char *data;
    int len;
};

struct read// Used to read sequences
{
    char *data;
    int *length;
    lint *start;
    struct counter *counter;
    struct read *next;
};

struct counter
{
    char **index;
    char *frequency;
//    int **index;
//    int *frequency;
};

#endif
