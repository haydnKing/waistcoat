
#include "Python.h"

typedef struct {
    char *name, *seq, *qual;
} FastQSeq;

size_t FastQSeq_Write(FastQSeq *s, FILE *f);
size_t FastQSeq_Read(FILE * f, FastQSeq *s);
void FastQSeq_Free(FastQSeq *s);

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

