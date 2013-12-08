//Speedy preprocessing

#include "Python.h"

//modules
PyObject *os, *os_path, *settings, *stats;


// whether or not to output to console
PyObject* verbose;



// Function table
static PyMethodDef
module_functions[] = {
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
    settings = PyImport_ImportModule("settings");
    if(settings == NULL) return;
    stats = PyImport_ImportModule("statistics");
    if(stats == NULL) return;


    //expose verbose
    verbose = PyBool_FromLong(1L);
    if(!PyModule_AddObject(module, "verbose", verbose)) return;


}




