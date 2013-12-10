//Speedy preprocessing

#include "preprocess.h"

//modules
PyObject *os, *os_path, *settings, *stats, *tempfile;

const float BOOST_FACTOR = 2.0;
const float MATCH_THRESHOLD = 0.04;
const size_t MIN_LENGTH = 15;

// ****************************************************************
// -------------------------- Functions from modules --------------
// ****************************************************************

// os.path.join
const char* os_path_join2(const char* lhs, const char* rhs)
{
    PyObject *join = PyObject_GetAttrString(os_path, "join");
    PyObject *str = PyObject_CallFunction(join, "ss", lhs, rhs);
    if(str == NULL)
        return NULL;

    const char* ret = PyString_AsString(str);
    Py_DECREF(join);
    Py_DECREF(str);
    
    return ret;
}

FILE *tempfile_mkstemp3(const char* dir, const char* prefix, const char* suffix,
        const char** filename)
{
    PyObject *mkstemp = PyObject_GetAttrString(tempfile, "mkstemp");
    PyObject *out = PyObject_CallFunction(mkstemp, "sss", suffix, prefix, dir);
    if (out == NULL)
    {
        return NULL;
    }

    int ret;
    PyArg_ParseTuple(out, "is", &ret, filename);

    Py_DECREF(mkstemp);
    Py_DECREF(out);

    return fdopen(ret, "w");
}

FILE *tempfile_mkstemp2(const char* dir, const char* prefix,
        const char** filename)
{
    return tempfile_mkstemp3(dir, prefix, "", filename);
}
FILE *tempfile_mkstemp1(const char* dir, const char** filename)
{
    return tempfile_mkstemp3(dir, "", "", filename);
}
FILE *tempfile_mkstemp0(const char** filename)
{
    return tempfile_mkstemp3("", "", "", filename);
}

int stats_addvalues(PyObject* values)
{
    PyObject* addValues = PyObject_GetAttrString(stats, "addValues");
    if(addValues == NULL) return -1;

    PyObject* ret = PyObject_CallFunctionObjArgs(addValues, values, NULL);
    if(ret == NULL) return -1;

    return 0;
}

// ****************************************************************
// -------------------------- Utility funcs --------------------------
// ****************************************************************

size_t write_lines(const char* out, size_t width, FILE *ofile)
{
    size_t len = strlen(out), pos = 0, written = 0, tmp = 0;

    while(pos < len)
    {
        tmp = fwrite(out + pos, sizeof(char), 
                (width < (len-pos)) ? width : (len-pos), ofile);
        fputc('\n', ofile);
        written += tmp + 1;
        pos += tmp;
    }
    return written;
}

void _extract(const FastQSeq* s, const char* barcode_format, char extract,
        char* barcode)
{
    int i = 0, j = 0;
    for(i = 0; i < strlen(barcode_format); i++)
    {
        if(barcode_format[i] == extract)
        {
            barcode[j] = s->seq[i];
            j++;
        }
    }
    barcode[j] = '\0';
}

void get_barcode(const FastQSeq* s, const char* barcode_format, char* barcode)
{
    _extract(s,barcode_format,'B',barcode);
}

void get_umi(const FastQSeq* s, const char* barcode_format, char* barcode)
{
    _extract(s,barcode_format,'U',barcode);
}

// ****************************************************************
// -------------------------- Structures --------------------------
// ****************************************************************

FastQSeq *FastQSeq_New(void)
{
    FastQSeq *ret = malloc(sizeof(FastQSeq));
    ret->name = NULL;
    ret->seq = NULL;
    ret->qual = NULL;
    return ret;
}

void FastQSeq_Free(FastQSeq *s)
{
    free(s->name);
    free(s->seq);
    free(s->qual);
    s->name = NULL;
    s->seq = NULL;
    s->qual = NULL;
}

size_t FastQSeq_Read(FILE * f, FastQSeq **s)
{
    if(feof(f)) {return 0;}
    size_t read = 0, buffsize = 1024, pos, len;
    char * buff = malloc(buffsize);

    //read the name
    if(fgets(buff, buffsize, f) == NULL) {return 0;}
    read += strlen(buff);
    if(buff[0] != '@') return 0;
    char* name = malloc(strlen(buff));
    strcpy(name, buff+1);

    //read the seq
    pos = 0;
    while(!feof(f) && !ferror(f))
    {
        if(fgets(buff+pos, buffsize-pos, f) == NULL) {return 0;};
        if(buff[pos] == '+') 
        {
            buff[pos] = '\0';
            break;
        }
        pos = strlen(buff);
        if((buffsize - pos) == 0 || (buffsize - pos) > buffsize)
        {
            free(buff);
            free(name);
            return 0;
        }
    }
    read += pos;
    len = pos;
    char *seq = malloc(len);
    strcpy(seq, buff);

    //read the quality
    pos = 0;
    while(!feof(f) && !ferror(f))
    {
        if(fgets(buff+pos, buffsize-pos, f) == NULL) {return 0;}
        pos = strlen(buff);
        if(pos >= len) break;
    }
    read += pos;
    char* qual = malloc(len);
    strcpy(qual, buff);

    free(buff);

    //remove trailing \n
    name[strlen(name)-1] = '\0';
    seq [strlen(seq )-1] = '\0';
    qual[strlen(qual)-1] = '\0';

    //create the object
    *s = FastQSeq_New();
    (*s)->name = name;
    (*s)->seq = seq;
    (*s)->qual = qual;
    return read;
}

size_t FastQSeq_Write(FastQSeq *s, FILE *f)
{
    size_t written = 0;
    written += fprintf(f, "@%s\n", s->name);
    written += write_lines(s->seq, 80, f);
    written += fprintf(f, "+\n");
    written += write_lines(s->qual, 80, f);

    return written;
}

float FastQSeq_Distance(FastQSeq *lhs, FastQSeq *rhs)
{
    int llen = strlen(lhs->seq), rlen = strlen(rhs->seq), i;
    int len = (llen > rlen) ? rlen : llen;
    if(len == 0){
        return -1;
    }
    float total = 0.0;
    for(i = 0; i < len; i++)
    {
        if(toupper(lhs->seq[i]) != toupper(rhs->seq[i]))
            total += 1.0;
    }

    return total / (float)len;
}

void FastQSeq_RemoveA(FastQSeq *self)
{
    int len = strlen(self->seq), i=0;
    for(i = len-1; i >=0; i--)
    {
        if(toupper(self->seq[i]) != 'A')
            break;
    }
    if(i < len-1)
    {
        self->seq[i+1] = '\0';
        self->qual[i+1] = '\0';
        self->seq = realloc(self->seq, i+1);
        self->qual= realloc(self->qual, i+1);
    }
}

float FastQSeq_Score(FastQSeq *self, int target_length)
{
    float score = 0;
    int len = strlen(self->seq), i;
    if(len == 0) return -1.0;
    for(i = 0; i < len; i++)
    {
        score += (float) self->qual[i];
    }
    score = score / (float) len;

    //give a boost to ones which are 28 long
    if(len == target_length){
        score = score * BOOST_FACTOR;
    }

    return score;
}

void FastQSeq_Offset(FastQSeq *self, int offset)
{
    self->seq = self->seq + offset;
    self->qual = self->qual + offset;
}

ConflictEl *ConflictEl_New(FastQSeq *seq)
{
    ConflictEl *ret = malloc(sizeof(ConflictEl));
    ret->seq = seq;
    ret->next = NULL;
    return ret;
}

void ConflictEl_Free(ConflictEl *el)
{
    FastQSeq_Free(el->seq);
    if(el->next != NULL)
        ConflictEl_Free(el->next);
    free(el);
}

void ConflictEl_Append(ConflictEl* self, ConflictEl* rhs)
{
    rhs->next = self->next;
    self->next = rhs;
}

void ConflictEl_Remove(ConflictEl* self)
{
    ConflictEl *n = self->next;
    FastQSeq_Free(self->seq);
    self->seq = n->seq;
    n->seq = NULL;
    self->next = n->next;
    n->next = NULL;
    free(n);
}

void Conflict_AppendNew(Conflict* self, FastQSeq *seq)
{
    Conflict *new = Conflict_New(seq);
    new->next = self->next;
    self->next = new;
}

void Conflict_Free(Conflict* self)
{
    ConflictEl_Free(self->first_element);
    self->first_element = NULL;
    free(self);
}

Conflict *Conflict_New(FastQSeq* seq)
{
    Conflict *ret = malloc(sizeof(Conflict));
    ret->first_element = ConflictEl_New(seq);
    ret->next = NULL;
    return ret;
}

FastQSeq *Conflict_Resolve(Conflict *self, int target_length)
{
    float max_score = 0.0, score;
    ConflictEl *el = self->first_element, *winner = NULL;
    while(el != NULL)
    {
        score = FastQSeq_Score(el->seq, target_length);
        if(score > max_score) 
        {
            max_score = score;
            winner = el;
        }
        el = el->next;
    }
    return winner->seq;
}

// ****************************************************************
// -------------------------- Module Functions --------------------
// ****************************************************************

// whether or not to output to console
PyObject* verbose;

PyObject* split_by_barcode(PyObject *self, PyObject *args)
{
    if(verbose == Py_True)
    {
        printf("Splitting Sequence by Barcode\n");
    }
    const char* in_file = NULL, * out_dir = NULL;
    PyObject* my_settings = NULL, * files, * barcodes, * count;
    int remove_input = 0;
    int ok = PyArg_ParseTuple(args, "sOs|i", &in_file, &my_settings, &out_dir,
            &remove_input);
    if(!ok){
        return NULL;
    }

    if(verbose == Py_True)
    {
        printf("Reading from \"%s\"\n", in_file);
    }

    //get sample names & barcodes
    int num_samples = 0;
    barcodes = PyObject_GetAttrString(my_settings, "barcodes");
    if (barcodes == NULL){ return NULL; }
    if (!PyDict_Check(barcodes)){
        PyErr_SetString(PyExc_TypeError, "settings.barcodes must be a dict");
        return NULL;
    }
    num_samples = PyDict_Size(barcodes);

    if(verbose == Py_True)
    {
        printf("Processing \"%d\" samples\n", num_samples);
    }

    //get barcode format
    PyObject *bfmt = PyObject_GetAttrString(my_settings, "barcode_format");
    if(bfmt == NULL) return NULL;
    const char* barcode_format = PyString_AsString(bfmt);
    Py_DECREF(bfmt);

    //open files
    files = PyDict_New();
    count = PyDict_New();
    const char* barcode_seqs[num_samples];
    PyObject* sample_names[num_samples];
    FILE *open_files[num_samples];
    Py_ssize_t pos = 0;
    int i = 0;
    PyObject *isample, *ibarcode;
    while(PyDict_Next(barcodes, &pos, &isample, &ibarcode))
    {
        //get the barcode
        barcode_seqs[i] = PyString_AsString(ibarcode);
        if(barcode_seqs[i] == NULL) return NULL;
        //get the sample name
        sample_names[i] = isample;
        Py_INCREF(isample);
        
        //open the temp file
        const char* filename;
        const char* ssample = PyString_AsString(isample);
        char sample_dot[strlen(ssample) + 2];
        sprintf(sample_dot, "%s.", ssample);
        open_files[i] = tempfile_mkstemp3(out_dir, sample_dot, ".barcode", 
                &filename);
        if(open_files[i] <= 0) return NULL;

        //store the filename
        PyObject* pFilename = PyString_FromString(filename);
        PyDict_SetItem(files, isample, pFilename);
        Py_DECREF(pFilename);

        i = i + 1;
        
        //initiate count for statistics
        PyObject* tZero = PyInt_FromLong(0);
        PyDict_SetItem(count, isample, tZero);
        Py_DECREF(tZero);

    }

    
    //open input file (possibly unziping)
    FILE *in = fopen(in_file, "rb");
    char barcode[strlen(barcode_format)+1];

    //for each seq
    FastQSeq * seq;
    PyObject * c;
    while(!feof(in) && !ferror(in)){
        int read = FastQSeq_Read(in, &seq);
        if(read == 0)
        {
            if(feof(in)) {break;}
            
            char err[64+strlen(in_file)];
            if(ferror(in))
            {
                sprintf(err, "Error reading from file \"%s\"", in_file);
                PyErr_SetString(PyExc_IOError, err);
                return NULL;
            }
            sprintf(err, "Error parsing file \"%s\"", in_file);
            PyErr_SetString(PyExc_IOError, err);
            return NULL;
        }
        
        //extract barcode
        //write it to the correct file

        get_barcode(seq, barcode_format, barcode);

        for(i=0; i < num_samples; i++)
        {
            if(strcmp(barcode, barcode_seqs[i]) == 0)
            {
                FastQSeq_Write(seq, open_files[i]);
                c = PyDict_GetItem(count, sample_names[i]);
                c = PyInt_FromLong(1 + PyInt_AsLong(c));
                PyDict_SetItem(count, sample_names[i], c);
                Py_DECREF(c);
                break;
            }
        }


        FastQSeq_Free(seq);
    }
    fclose(in);

    

    //close files
    for(i=0; i < num_samples; i++)
    {
        fclose(open_files[i]);
        Py_DECREF(sample_names[i]);
    }

    Py_DECREF(barcodes);


    if(!stats_addvalues(count))
    {
        Py_DECREF(count);
        return NULL;
    }
    Py_DECREF(count);

    if(remove_input)
    {
        if(!remove(in_file))
        {
            char e[64+strlen(in_file)];
            sprintf(e, "Could not remove input file \"%s\"", in_file);
            PyErr_SetString(PyExc_IOError, e);
            return NULL;
        }
    }

    return files;

}

PyObject *process_sample(PyObject* self, PyObject *args)
{
    //parse arguments
    PyObject *in_files = NULL, *my_settings = NULL, *ptemp;
    const char* outdir = NULL;
    int delete_input = 1;
    int ok = PyArg_ParseTuple(args, "O!Os|d", &PyDict_Type, &in_files,
            &my_settings, &outdir, &delete_input);
    if(!ok)
    {
        return NULL;
    }


    //extract barcode format
    ptemp = PyObject_GetAttrString(my_settings, "barcode_format");
    if(ptemp == NULL) return NULL;
    const char* barcode_format = PyString_AsString(ptemp);
    if(barcode_format == NULL) return NULL;
    Py_DECREF(ptemp);
    size_t barcode_length = strlen(barcode_format);


    //prepare output dict
    PyObject* out_files = PyDict_New();

    //for each sample and file
    Py_ssize_t pos = 0;
    PyObject *isample, *ifile;
    long count;
    while(PyDict_Next(in_files, &pos, &isample, &ifile))
    {
        const char* in_name = PyString_AsString(ifile), *out_name;


        //open input file
        FILE *in = fopen(in_name, "rb");

        //load in all seqs
        FastQSeq *seq = NULL;
        Conflict *conflict_first = NULL, *conflict_last, *conflict;
        float min, current;
        ConflictEl *target;
        while(!feof(in))
        {
            int read = FastQSeq_Read(in, &seq);
            if(read == 0)
            {
                if(feof(in)) {break;}
                
                char err[64+strlen(in_name)];
                if(ferror(in))
                {
                    sprintf(err, "Error reading from file \"%s\"", in_name);
                    PyErr_SetString(PyExc_IOError, err);
                    Py_DECREF(out_files);
                    return NULL;
                }
                sprintf(err, "Error parsing file \"%s\"", in_name);
                PyErr_SetString(PyExc_IOError, err);
                Py_DECREF(out_files);
                return NULL;
            }
            FastQSeq_RemoveA(seq);
            if((strlen(seq->seq)-barcode_length) < MIN_LENGTH)
            {
                continue;
            }

            conflict = conflict_first;
            //first case
            if(conflict_first == NULL)
            {
                conflict_first = Conflict_New(seq);
                conflict_last = conflict_first;
                continue;
            }
            //find the closest matching conflict
            target = NULL;
            min = 10.0 * MATCH_THRESHOLD;
            while(conflict != NULL)
            {
                current = FastQSeq_Distance(seq, conflict->first_element->seq);
                if(current < 0.0)
                {
                    conflict = conflict->next;
                    continue;
                }
                if(current < min)
                {
                    min = current;
                    target = conflict->first_element;
                }
                conflict = conflict->next;
            }

            //if there is a matching conflict
            if(target != NULL)
            {
                ConflictEl_Append(target, ConflictEl_New(seq));
            }
            else //if no conflict is close enough
            {
                Conflict_AppendNew(conflict_last, seq);
                conflict_last = conflict_last->next;
            }
        }

        //close input
        fclose(in);
        
        //create and open output file
        const char* ssample = PyString_AsString(isample);
        char sample_dot[strlen(ssample)+1];
        sprintf(sample_dot, "%s.", ssample);
        FILE *out = tempfile_mkstemp3(outdir, sample_dot, ".clean", &out_name);
        
        //store output file
        ptemp = PyString_FromString(out_name);
        PyDict_SetItem(out_files, isample, ptemp);
        Py_DECREF(ptemp);

        //resolve conflicts and save the winner
        conflict = conflict_first;
        count = 0;
        while(conflict != NULL)
        {
            seq = Conflict_Resolve(conflict, barcode_length + 28);
            FastQSeq_Offset(seq, barcode_length);
            FastQSeq_Write(seq, out);
            FastQSeq_Offset(seq, -barcode_length);
            count += 1;

            conflict_first = conflict->next;
            Conflict_Free(conflict);
            conflict = conflict_first;
        }

        //close output
        fclose(out);

        //delete input
        if(delete_input)
        {
            if(!remove(in_name))
            {
                char e[32+strlen(in_name)];
                sprintf(e, "Failed to remove file \"%s\"", in_name);
                PyErr_SetString(PyExc_IOError, e);
            }
        }
    }


    return out_files;
}

// Function table
static PyMethodDef
module_functions[] = {
    {"process_sample", process_sample, METH_VARARGS,
        "process_sample(files, my_settings, outdir, remove_input=True)\n"
            "  Clean samples and remove duplicates\n"
            "    files: dict mapping sample names to files\n"
            "    my_settings: settings.Settings object\n"
            "    outdir: directory to output to\n"
            "    remove_input: whether to remove the input files [True]"},
    {"split_by_barcode", split_by_barcode, METH_VARARGS,
        "split_by_barcode(filename, my_settings, outdir, remove_input=False)"
        " Split the fastq sequences found in filename into seperate files"
        " defined by my_settings, saving the files in outdir"},
    { NULL } //sentinel
};

void
initpreprocess(void)
{
    PyObject* module;

    module = Py_InitModule3("preprocess", module_functions,
            "Preprocess sequencing reads");

    //import other modules
    os = PyImport_ImportModule("os");
    if(os == NULL) return;
    os_path = PyImport_ImportModule("os.path");
    if(os_path == NULL) return;
    tempfile = PyImport_ImportModule("tempfile");
    if(tempfile == NULL) return;
    settings = PyImport_ImportModule("waistcoat.settings");
    if(settings == NULL) return;
    stats = PyImport_ImportModule("waistcoat.statistics");
    if(stats == NULL) return;

    //expose verbose
    verbose = PyBool_FromLong(1L);
    if(PyModule_AddObject(module, "verbose", verbose) < 0) return;

}




