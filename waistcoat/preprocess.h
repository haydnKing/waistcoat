
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
float FastQSeq_Distance(FastQSeq *lhs, FastQSeq *rhs);
void FastQSeq_RemoveA(FastQSeq *self);
float FastQSeq_Score(FastQSeq *self, int target_length);
void FastQSeq_Offset(FastQSeq *self, int offset);

typedef struct ConflictEl ConflictEl;

struct ConflictEl {
    FastQSeq *seq;
    ConflictEl *next;
};

//make a new conflict element
ConflictEl *ConflictEl_New(FastQSeq *seq);
//free a conflict element, its seq, and any subsequent elements
void ConflictEl_Free(ConflictEl *el);
//Add an element after self
void ConflictEl_Append(ConflictEl* self, ConflictEl* rhs);
//remove an element - the address is then filled by next
void ConflictEl_Remove(ConflictEl* self);

typedef struct Conflict Conflict;

struct Conflict {
    ConflictEl *first_element;
    Conflict *next;
};

//create a new conflict for seq and append it to the list
void Conflict_AppendNew(Conflict* self, FastQSeq *seq);
//Free the conflict, and its elements
void Conflict_Free(Conflict* self);
//Create a new conflict
Conflict *Conflict_New(FastQSeq* seq);
//resolve a conflict
FastQSeq *Conflict_Resolve(Conflict *self, int target_length);


