
#include "Python.h"

//utilities
void print_read_count(PyObject* count, long total, int indent);
void print_read_dist(long *dist, size_t length, int width, int indent);


//structs
typedef struct {
    char *name, *seq, *qual;
} FastQSeq;

FastQSeq *FastQSeq_New(void);
size_t FastQSeq_Write(FastQSeq *s, FILE *f);
size_t FastQSeq_Read(FILE * f, FastQSeq **s);
void FastQSeq_Free(FastQSeq *s);
void FastQSeq_RemoveA(FastQSeq *self);
void FastQSeq_RemoveBarcode(FastQSeq *self, int barcode_length);

void get_barcode(const FastQSeq* s, const char* barcode_format, char* barcode);
void get_umi(const FastQSeq* s, const char* barcode_format, char* barcode);
int get_umi_length(const char* barcode_format);
long get_umi_long(const char* umi);


typedef struct SeqItem SeqItem;

struct SeqItem {
    FastQSeq *the_seq;
    SeqItem *next;
};

SeqItem *SeqItem_New(FastQSeq *seq);
void SeqItem_Free(SeqItem *self);
void SeqItem_Append(SeqItem* self, SeqItem* rhs);
//try and replace self with rhs, return 1 if they match
int SeqItem_TryMerge(SeqItem* self, FastQSeq* rhs);
//iterate through the list, return True if there is a next, otherwise
// next is unchanged
int SeqItem_Next(SeqItem **next);
size_t SeqItem_Length(SeqItem *s);
