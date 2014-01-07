waistcoat
=========

A pipeline for mapping RNA-seq data to a genome during Ribosome Profiling

RNA-seq data in FASTQ format (optionally gziped) is processed using the
following method:
* Reads are split depending on their barcode and processed separately
* Non-unique reads are discarded
* Barcodes and UMIs are removed
* Reads which map to any of a specified list of Bowtie indexes are discarded
  (e.g. rRNA)
* Reads are mapped to a Bowtie index of the genome using Tophat

Settings for configuring the barcode format and Tophat settings are given in a
JSON encoded settings file, a template for which is given below.

### Installation

Waistcoat is simplest to use in place, but the C-extension must first be built
by calling `python setup.by build_ext --inplace`.

You can also run the test suite using `python setup.py test`.

### Example Settings

```JSON
{
"_comment": "This is an example settings file for waistcoat. Comments begin with underscores",

"_barcode_format": "REQUIRED. Format of the barcode - B = Barcode character N = Random (UMI) character",
"barcode_format" : "BBBNNNNBB",

"_barcodes": "REQUIRED. A dictionary mapping sample names to barcodes",
"barcodes" : {
	"sample 1": "ACCTA",
	"sample 2": "GCGAT"
	},

"_discard": "OPTIONAL. Sequences which map to indexes listed here will be discarded",
"discard": [
	"path/to/discard_index_base1",
	"discard_index_base2"
	],

"_discard_settings": "OPTIONAL. Default settings for tophat when discarding",
"discard_settings": {
	"_comment": "settings go here as key value pairs, e.g. this sets --max-insertion-length 5",
	"max_insertion_length": 5
},

"_discard_index_base1_settings": "OPTIONAL. Override default discard settings for a specific index by setting [index_name]_settings",
"discard_index_base1_settings": {
	"_comment": "override settings here"
	},

"_target": "REQUIRED. Index to perform final map against", 
"target": "final_ref",

"_target_settings": "OPTIONAL. Settings for final mapping",
"target_map_settings": {
	"_comment": "Settings go here"
	}
}
```
