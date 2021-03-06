//Speedy preprocessing

#include "preprocess.h"

#define VAL_A 0
#define VAL_T 1
#define VAL_C 2
#define VAL_G 3

//modules
PyObject *os, *os_path, *settings, *stats, *tempfile, *the_module;

const float BOOST_FACTOR = 2.0;
const float MATCH_THRESHOLD = 0.04;
const size_t MIN_LENGTH = 15;
const size_t LENGTH_DIST = 512;

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

int stats_addvalues(const char* name, PyObject* values)
{
    PyObject* addValues = PyObject_GetAttrString(stats, "addValues");
    if(addValues == NULL) return 0;

    PyObject* ret = PyObject_CallFunction(addValues, "sO", name, values);
    if(ret == NULL) return 0;

    return 1;
}

// ****************************************************************
// -------------------------- Utility funcs --------------------------
// ****************************************************************

int is_verbose(void)
{
    PyObject *v = PyObject_GetAttrString(the_module, "verbose");
    int r = (int) PyInt_AsLong(v);
    Py_DECREF(v);
    return r;
}

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
    _extract(s,barcode_format,'N',barcode);
}

int get_umi_length(const char* barcode_format)
{
    int l = 0, i;
    for(i=0; i < strlen(barcode_format); i++)
    {
        if(barcode_format[i] == 'N')
            l++;
    }
    return l;
}

long get_umi_long(const char* umi)
{
    int i;
    long r = 0;
    char c;
    for(i = 0; i < strlen(umi); i++)
    {
        c = toupper(umi[i]);
        r = r * 4;
        switch(c)
        {
            case 'A':
                r += VAL_A;
                break;
            case 'T':
                r += VAL_T;
                break;
            case 'C':
                r += VAL_C;
                break;
            case 'G':
                r += VAL_G;
                break;
        }
    }
    return r;
}


void print_read_count(PyObject* count, long total, int indent)
{
    Py_ssize_t pos = 0;
    int i;
    long current;
    PyObject *isample, *icount;
    while(PyDict_Next(count, &pos, &isample, &icount))
    {
        for(i=0; i < indent; i++)
            putchar('\t');
        current = PyInt_AsLong(icount);
        printf("%s : %3.1f%% : %ld\n", PyString_AsString(isample),
                (double)current / (double)total * 100.0, current);
    }
}

void print_read_dist(long *dist, size_t length, int width, int indent)
{
    //find max and start and end points
    int j;
    size_t i = 0, start = length, end=0;
    long max = 0L;
    for(i=0;i<length; i++)
    {
        if( (dist[i] > 0L) && (i < start) )
        {
            start = i;
        }
        if( (dist[i] > 0L) && (i > end) )
        {
            end = i;
        }
        if( dist[i] > max )
        {
            max = dist[i];
        }
    }

    double scale = (double)width / (double)max;
    for(i=start; i <= end; i++)
    {
        for(j=0; j < indent; j++) putchar('\t');
        printf("%3ld : ", i);
        for(j=0; j < (int) (scale * (double) dist[i]); j++) putchar('*');
        putchar('\n');
    }

    for(j=0; j < indent+1; j++) putchar('\t');
    printf("* = %.1f\n", 1.0 / scale);
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
    ret->low = ret->high = 0LL;
    return ret;
}

void FastQSeq_Free(FastQSeq *s)
{
    if(s == NULL) return;
    free(s->name);
    free(s->seq);
    free(s->qual);
    s->name = NULL;
    s->seq = NULL;
    s->qual = NULL;
    free(s);
}

size_t FastQSeq_Read(FILE * f, FastQSeq **s)
{
    if(feof(f)) {return 0;}
    size_t read = 0, buffsize = 1024, pos, len;
    char * buff = malloc(buffsize);

    //read the name - find the next line begining with @
    buff[0] = '0'; //something that's not '@'
    while(buff[0] != '@')
    {
        if(fgets(buff, buffsize, f) == NULL) {return 0;}
        read += strlen(buff);
    }
    char* name = malloc(strlen(buff));
    strcpy(name, buff+1); 
    name[strlen(name)-1] = '\0';

    //read the seq
    pos = 0;
    buff[0] = '\n';
    buff[1] = '\0';
    while(buff[pos] != '+')
    {
        //rewind over the newline
        pos = strlen(buff) - 1;
        if(fgets(buff+pos, buffsize-pos, f) == NULL) {return 0;};
    }
    //clean off the '+'
    buff[pos] = '\0';
    len = strlen(buff);
    read += len;
    char *seq = malloc(len+1);
    strcpy(seq, buff);

    //read the quality
    pos = 0;
    while(pos < len)
    {
      if(fgets(buff+pos, buffsize-pos, f) == NULL) {return 0;}
      //reverse over the \n
      pos = strlen(buff) - 1;
    }
    //null-terminate
    buff[len] = '\0';
    read += pos;
    char* qual = malloc(len+1);
    strcpy(qual, buff);

    free(buff);

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
    }
}

void FastQSeq_RemoveBarcode(FastQSeq *self, int barcode_length)
{
    size_t len = strlen(self->seq);
    if(barcode_length < 0) return;
    if(barcode_length > len) return;

    int i;
    for(i = 0; i < len - barcode_length; i++)
    {
        self->seq[i] = self->seq[i+barcode_length];
        self->qual[i] = self->qual[i+barcode_length];
    }
    self->seq[i] = '\0';
    self->qual[i]= '\0';
}

void FastQSeq_SetBits(FastQSeq *self)
{
    self->high = self->low = 0LL;
    size_t pos, 
           len = strlen(self->seq),
           max_high = (len > 32L) ? 32L : len,
           max_low  = (len > 64L) ? 64L : len;
    max_low = (max_low < 32L) ? 0L : max_low - 32L;

    for(pos = 0; pos < max_high; pos++)
    {
        self->high *= 4LL;
        switch(toupper(self->seq[pos]))
        {
            case 'A':
                self->high += VAL_A;
                break;
            case 'T':
                self->high += VAL_T;
                break;
            case 'C':
                self->high += VAL_C;
                break;
            case 'G':
                self->high += VAL_G;
                break;
        }
    }

    for(pos = 0; pos < max_low; pos++)
    {
        self->low *= 4LL;
        switch(toupper(self->seq[pos+32]))
        {
            case 'A':
                self->low += VAL_A;
                break;
            case 'T':
                self->low += VAL_T;
                break;
            case 'C':
                self->low += VAL_C;
                break;
            case 'G':
                self->low += VAL_G;
                break;
        }
    }

}

// -------------------------- SeqItem


SeqItem *SeqItem_New(FastQSeq *seq)
{
    SeqItem *r = malloc(sizeof(SeqItem));
    r->the_seq = seq;
    r->next = NULL;
    r->prev = NULL;
    return r;
}

void SeqItem_Free(SeqItem *self)
{
    if(self == NULL) return;
    FastQSeq_Free(self->the_seq);
    SeqItem_Free(self->next);
    free(self);
    return;
}

void SeqItem_Append(SeqItem* self, SeqItem* rhs)
{
    //links with self->next
    rhs->next = self->next;
    if(rhs->next != NULL)
    {
        rhs->next->prev = rhs;
    }
    //links with self
    rhs->prev = self;
    self->next = rhs;
}

void SeqItem_Prepend(SeqItem **head, SeqItem* self, SeqItem* rhs)
{
    //links with self->prev
    rhs->prev = self->prev;
    if(self->prev != NULL)
    {
        self->prev->next = rhs;
    }
    else
    {
        //we replaced head
        (*head) = rhs;
    }
    //links with self
    rhs->next = self;
    self->prev = rhs;
}

int SeqItem_Next(SeqItem **next)
{
    if((*next) == NULL) return 0;
    if((*next)->next == NULL)
        return 0;
    (*next) = (*next)->next;
    return 1;
}

size_t SeqItem_Length(SeqItem *s)
{
    if(s == NULL) return 0;
    size_t r = 1;
    while(SeqItem_Next(&s))
        r++;
    return r;
}

// ****************************************************************
// -------------------------- Module Functions --------------------
// ****************************************************************


PyObject* split_by_barcode(PyObject *self, PyObject *args)
{
    if(is_verbose())
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

    if(is_verbose())
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

    if(is_verbose())
    {
        printf("Processing \"%d\" samples\n", num_samples);
    }

    //get barcode format
    PyObject *bfmt = PyObject_GetAttrString(my_settings, "barcode_format");
    if(bfmt == NULL) return NULL;
    const char* barcode_format = PyString_AsString(bfmt);
    Py_DECREF(bfmt);


    //open output files
    files = PyDict_New();
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
        const char* filename = NULL;
        const char* ssample = PyString_AsString(isample);
        char sample_dot[strlen(ssample) + 2];
        sprintf(sample_dot, "%s.", ssample);
        open_files[i] = tempfile_mkstemp3(out_dir, sample_dot, ".barcode", 
                &filename);
        if(open_files[i] <= 0)
        {
            Py_DECREF(barcodes);
            Py_DECREF(files);
            return NULL;
        }

        //store the filename
        PyObject* pFilename = PyString_FromString(filename);
        PyDict_SetItem(files, isample, pFilename);
        Py_DECREF(pFilename);

        i = i + 1;
    }
    Py_DECREF(barcodes);

    //open input file
    FILE *in = fopen(in_file, "rb");
    if(in == NULL)
    {
        char e[32+strlen(in_file)];
        sprintf(e, "could not open input file \"%s\"", in_file);
        PyErr_SetString(PyExc_IOError, e);
        Py_DECREF(files);
        return NULL;
    }


    char barcode[strlen(barcode_format)+1];

    //for each seq
    FastQSeq * seq = NULL;
    long total = 0;
    long ccount[num_samples];
    for(i=0; i<num_samples; i++)
        ccount[i] = 0L;

    while(FastQSeq_Read(in, &seq))
    {
        //extract barcode
        //write it to the correct file

        get_barcode(seq, barcode_format, barcode);

        for(i=0; i < num_samples; i++)
        {
            if(strcmp(barcode, barcode_seqs[i]) == 0)
            {
                FastQSeq_Write(seq, open_files[i]);
                ccount[i] += 1;
                total += 1;
                break;
            }
        }

        FastQSeq_Free(seq);
        seq = NULL;
    }
    //close files and fill in count and filename
    count = PyDict_New();
    for(i=0; i < num_samples; i++)
    {
        fclose(open_files[i]);
        PyObject *c = PyInt_FromLong(ccount[i]);
        PyDict_SetItem(count, sample_names[i], c);
        
        PyObject *f = PyDict_GetItem(files, sample_names[i]);
        PyObject *t = PyTuple_Pack(2, f, c);
        PyDict_SetItem(files, sample_names[i], t);

        Py_DECREF(t);
        Py_DECREF(c);
        Py_DECREF(sample_names[i]);
    }
    //check for error
    if(ferror(in))
    {
        char err[64+strlen(in_file)];
        sprintf(err, "Error reading from file \"%s\"", in_file);
        PyErr_SetString(PyExc_IOError, err);
        Py_DECREF(count);
        Py_DECREF(files);
        return NULL;
    }
    fclose(in);

    if(is_verbose())
    {
        printf("Found %ld reads:-\n", total);
        print_read_count(count, total, 1);
    }

    //store statistics
    if(!stats_addvalues("split_by_barcode",count))
    {
        Py_DECREF(count);
        Py_DECREF(files);
        return NULL;
    }
    Py_DECREF(count);

    //delete input file if requested
    if(remove_input)
    {
        if(remove(in_file))
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
    if(is_verbose())
    {
        printf("Processing Samples\n");
    }
    //parse arguments
    PyObject *in_files = NULL, *my_settings = NULL, *ptemp;
    const char* out_dir = NULL;
    int remove_input = 1, i;
    int ok = PyArg_ParseTuple(args, "O!Os|i", &PyDict_Type, &in_files,
            &my_settings, &out_dir, &remove_input);
    if(!ok)
    {
        return NULL;
    }
    //check dict
    Py_ssize_t pos = 0;
    PyObject *isample = NULL, *ifile = NULL;
    while(PyDict_Next(in_files, &pos, &isample, &ifile))
    {
        if(!PyString_Check(isample))
        {
            PyErr_SetString(PyExc_ValueError, "Sample name must be a string");
            return NULL;
        }
        const char *fname = NULL;
        long length = 1;
        int ok = PyArg_ParseTuple(ifile, "sl", &fname, &length);
        if(!ok)
        {
            return NULL;
        }
    }

    long length_dist[LENGTH_DIST];
    for(i = 0; i < LENGTH_DIST; i++)
        length_dist[i] = 0;



    //extract barcode format
    ptemp = PyObject_GetAttrString(my_settings, "barcode_format");
    if(ptemp == NULL) return NULL;
    const char* barcode_format = PyString_AsString(ptemp);
    if(barcode_format == NULL) return NULL;
    Py_DECREF(ptemp);
    size_t barcode_length = strlen(barcode_format);
    long umi_length = get_umi_length(barcode_format);
    long num_umi = pow(4, umi_length);


    //prepare output dict
    PyObject* out_files = PyDict_New();

    //for each sample and file
    pos = 0;
    long count, total = 0;
    PyObject *PyCount = PyDict_New();
    SeqItem **UMI = malloc(num_umi * sizeof(SeqItem*));
    for(i=0; i < num_umi; i++)
        UMI[i] = NULL;
    while(PyDict_Next(in_files, &pos, &isample, &ifile))
    {
        const char *in_name = NULL, *out_name = NULL;
        long length = 1;
        int ok = PyArg_ParseTuple(ifile, "sl", &in_name, &length);
        if(!ok)
        {
            Py_DECREF(out_files);
            Py_DECREF(PyCount);
            return NULL;
        }

        if(is_verbose())
        {
            const char *sample_name = PyString_AsString(isample);
            printf("Sample \"%s\"\n", sample_name);
            printf("\tReading from \"%s\"\n", in_name);
        }

        //open input file
        FILE *in = fopen(in_name, "rb");
        if(in == NULL)
        {
            char e[32+strlen(in_name)];
            sprintf(e, "Could not open file \"%s\"", in_name);
            PyErr_SetString(PyExc_IOError, e);
            Py_DECREF(PyCount);
            Py_DECREF(out_files);
            return NULL;
        }

        //load in all seqs
        FastQSeq *seq = NULL;
        SeqItem *item;
        char umi[umi_length];
        long lumi = 0;
        count = 0;
        while(FastQSeq_Read(in, &seq))
        {
            count += 1;
            if(count % 1000 == 0)
            {
                if(is_verbose())
                {
                    printf("\rReading %ld/%ld (%3.1f%%)         ", 
                        count, length,
                        100.0 * (float)((double)count / (double)length));
                    fflush(stdout);
                }
            }
            FastQSeq_RemoveA(seq);
            if((strlen(seq->seq)-barcode_length) < MIN_LENGTH)
            {
                FastQSeq_Free(seq);
                seq = NULL;
                continue;
            }
            get_umi(seq, barcode_format, umi);
            //Trim the barcode
            FastQSeq_RemoveBarcode(seq, barcode_length);
            //set seq->high and seq->low
            FastQSeq_SetBits(seq);


            //check sequences with the same UMI
            lumi = get_umi_long(umi);
            item = UMI[lumi];

            //first case
            if(item == NULL)
            {
                UMI[lumi] = SeqItem_New(seq);
                continue;
            }
            //for each seq with the same UMI
            while(1)
            {
                //compare high
                if(seq->high == item->the_seq->high)
                {
                    if(seq->low == item->the_seq->low)
                    {
                        //we have a match
                        FastQSeq_Free(seq);
                        seq = NULL;
                        break;
                    }
                    else if(seq->low > item->the_seq->low)
                    {
                        //greater than
                        SeqItem_Prepend(UMI+lumi, item, SeqItem_New(seq));
                        break;
                    }

                }
                else if(seq->high > item->the_seq->high)
                {
                    //greater than
                    SeqItem *new = SeqItem_New(seq);
                    SeqItem_Prepend(UMI+lumi, item, new);
                    break;
                }
                // if less than, we need to compare with the next seq
                //move to the next seq
                if(item->next != NULL)
                    item = item->next;
                else
                {
                    //seq is less than the last item
                    SeqItem_Append(item, SeqItem_New(seq));
                    break;
                }
            }
        }
        //check for error
        if(ferror(in))
        {
            char err[64+strlen(in_name)];
            sprintf(err, "Error reading from file \"%s\"", in_name);
            PyErr_SetString(PyExc_IOError, err);
            Py_DECREF(PyCount);
            Py_DECREF(out_files);
            return NULL;
        }

        //close input
        fclose(in);
        
        //create and open output file
        const char* ssample = PyString_AsString(isample);
        char sample_dot[strlen(ssample)+1];
        sprintf(sample_dot, "%s.", ssample);
        FILE *out = tempfile_mkstemp3(out_dir, sample_dot, ".clean", &out_name);
        
        //store output file
        ptemp = PyString_FromString(out_name);
        PyDict_SetItem(out_files, isample, ptemp);
        Py_DECREF(ptemp);

        if(is_verbose())
        {
            printf("\r                                                       \r");
            printf("\t\tWriting to \"%s\"\n", out_name);
        }

        //save each item
        count = 0;
        ok = 0;
        for(i=0; i < num_umi; i++)
        {
            item = UMI[i];
            while(item != NULL)
            {
                //write
                FastQSeq_Write(item->the_seq, out);

                //statistics
                int len = strlen(item->the_seq->seq);
                if(len < LENGTH_DIST)
                    length_dist[len] += 1;
                count += 1;
                total += 1;

                item = item->next;
            }
            //free all items in this UMI
            SeqItem_Free(UMI[i]);
            UMI[i] = NULL;
        }

        //close output
        fclose(out);

        //delete input
        if(remove_input)
        {
            if(remove(in_name))
            {
                char e[32+strlen(in_name)];
                sprintf(e, "Failed to remove file \"%s\"", in_name);
                PyErr_SetString(PyExc_IOError, e);
                Py_DECREF(out_files);
                Py_DECREF(PyCount);
                return NULL;
            }
        }

        //store the count
        ptemp = PyInt_FromLong(count);
        PyDict_SetItem(PyCount, isample, ptemp);
        Py_DECREF(ptemp);

        if(is_verbose())
        {
            printf("\t\tWritten %ld reads\n", count);
        }
    }
    free(UMI);
    
    //print summary
    if(is_verbose())
    {
        printf("Found %ld reads :-\n", total);
        print_read_count(PyCount, total, 1);
        printf("Length Distribution :-\n");
        print_read_dist(length_dist, LENGTH_DIST, 50, 1);
    }

    //save statistics
    if(!stats_addvalues("clean", PyCount))
    {
        Py_DECREF(out_files);
        return NULL;
    }


    return out_files;
}

PyObject* run(PyObject *self, PyObject *args)
{
    PyObject *my_settings = NULL, 
             *files1      = NULL, 
             *files2      = NULL, 
             *the_args    = NULL;
    const char *in_file   = NULL, 
               *out_dir    = NULL;
    int remove_input = 0;
    int ok = PyArg_ParseTuple(args, "sOs|i", &in_file, &my_settings, &out_dir,
            &remove_input);
    if(!ok)
    {
        return NULL;
    }

    //test that in_file exists?

    the_args = Py_BuildValue("sOsi", in_file, my_settings, out_dir, remove_input);
    //call split_by_barcode
    files1 = split_by_barcode(self, the_args);
    Py_DECREF(the_args);
    if(files1 == NULL) return NULL;

    //call process_sample -- always remove_input for internal tempfiles
    the_args = Py_BuildValue("OOsi", files1, my_settings, out_dir, 1);
    Py_DECREF(files1);
    files2 = process_sample(self, the_args);
    Py_DECREF(the_args);
    if(files2 == NULL) return NULL;

    return files2;
}

// Function table
static PyMethodDef
module_functions[] = {
    {"run", run, METH_VARARGS,
        "run(in_file, my_settings, out_dir, remove_input=True)\n"
            "  Run the preprocess pipeline\n"
            "   in_file: input file (fastQ format)\n"
            "   my_settings: Settings object\n"
            "   out_dir: directory to write output and temp files\n"
            "   remove_input: whether or not to remove in_file\n"
            "Returns:\n"
            "   dictionary mapping sample name to file name"},
    {"process_sample", process_sample, METH_VARARGS,
        "process_sample(files, my_settings, out_dir, remove_input=True)\n"
            "  Clean samples and remove duplicates\n"
            "    files: dict mapping sample names to files\n"
            "    my_settings: settings.Settings object\n"
            "    out_dir: directory to output to\n"
            "    remove_input: whether to remove the input files [True]"},
    {"split_by_barcode", split_by_barcode, METH_VARARGS,
        "split_by_barcode(filename, my_settings, out_dir, remove_input=False)"
        " Split the fastq sequences found in filename into seperate files"
        " defined by my_settings, saving the files in out_dir"},
    { NULL } //sentinel
};

void
initpreprocess(void)
{
    the_module = Py_InitModule3("preprocess", module_functions,
            "Preprocess sequencing reads");

    //import other modules
    os = PyImport_ImportModule("os");
    if(os == NULL) return;
    os_path = PyImport_ImportModule("os.path");
    if(os_path == NULL) return;
    tempfile = PyImport_ImportModule("tempfile");
    if(tempfile == NULL) return;
    settings = PyImport_ImportModule("settings");
    if(settings == NULL) 
    {
        PyErr_Clear();
        settings = PyImport_ImportModule("waistcoat.settings");
        if(settings == NULL) return;
    }
    stats = PyImport_ImportModule("statistics");
    if(stats == NULL)
    {
        PyErr_Clear();
        stats = PyImport_ImportModule("waistcoat.statistics");
        if(stats == NULL) return;
    }

    //expose verbose
    PyObject *verbose = PyBool_FromLong(1L);
    if(PyModule_AddObject(the_module, "verbose", verbose) < 0) return;
    Py_DECREF(verbose);

}




