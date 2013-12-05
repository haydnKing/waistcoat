import csv, os, os.path, pysam

recording = True
#for length recording
max_len = 128

#data
names = []
data = {}
lengths = {}

def clear():
	"""clear recorded statistics"""
	data.clear()
	lengths.clear()
	del names[0:len(names)]

def setUp(sample_names):
	"""set sample_names"""
	clear()
	for sample in sample_names:
		data[sample] = []


def addValues(name, new_data):
	"""Add a measurement"""
	if not recording:
		return
	if not isinstance(name, basestring):
		raise TypeError("name must be a string, not {}".format(type(name)))
	
	if sorted(new_data.keys()) != sorted(data.keys()):
		raise KeyError("new_data must have all the same keys as data, {} != {}"
				.format(list(sorted(new_data.keys())), list(sorted(data.keys()))))

	for sample, value in new_data.iteritems():
		data[sample].append(value)
	names.append(name)

def collectFinalStats(sample, samfile):
	"""Collect statistics from output samfile"""
	l = [0,] * max_len
	for alg in pysam.Samfile(samfile, 'rb'):
		if alg.alen >= 0 and alg.alen < max_len:
			l[alg.alen] += 1
	
	lengths[sample] = l

def prettyString():
	"""return a pretty representation of the statistics"""
	col_widths = []
	ret = []

	#find the maximum width of each column
	col_widths.append(max(len(str(x)) for x in ['sample',] + data.keys()))
	for i,name in enumerate(names):
		col_widths.append(max(len(str(x)) for x in 
				[names[i],] + [s[i] for s in data.values()]))

	hdr = ' '.join(["\t{{:<{}s}}".format(col_widths[0]),] + 
			["{{:>{}s}}".format(w) for w in col_widths[1:]])
	line = ' '.join(["\t{{:<{}s}}".format(col_widths[0]),] + 
			["{{:>{}d}}".format(w) for w in col_widths[1:]])

	ret.append('Pipeline Reads:')
	ret.append(hdr.format('sample', *names))
	for sample, values in data.iteritems():
		ret.append(line.format(sample, *values))

	l_start = max_len
	l_end = 0
	for sample,l in lengths.iteritems():
		s = next((i for i,x in enumerate(l) if x))
		if s < l_start: l_start = s
		e = len(l) - next((i for i,x in enumerate(reversed(l)) if x))
		if e > l_end: l_end = e

	maximum = max(sum(x[i] for x in lengths.values()) for i in range(max_len))
	if maximum:
		scale = 50.0 / float(maximum)
		ret.append('Read Lengths:')
		fmt = "\t{{:>{}d}}: {{}}".format(len(str(l_end)))
		for i in range(l_start, l_end):
			value = sum(x[i] for x in lengths.itervalues())
			ret.append(fmt.format(i, '*'*int(scale*value)))

		ret.append('\t\t* = {:.1f} reads'.format(1.0 / scale))

	return '\n'.join(ret)

def write(directory):
	"""Write statistics to a directory"""
	if os.path.isfile(directory):
		raise ValueError(
			"\'{}\' should be a directory, not a file".format(directory))
	if not os.path.exists(directory):
		os.mkdir(directory)
	
	with open(os.path.join(directory, 'pipeline.csv'), 'wb') as csvfile:
		out = csv.writer(csvfile, delimiter=' ')
		out.writerow(['sample'] + names)
		for sample, values in data.iteritems():
			out.writerow([sample,] + values)

	with open(os.path.join(directory, 'lengths.csv'), 'wb') as csvfile:
		out = csv.writer(csvfile, delimiter=' ')
		out.writerow(['sample'] + range(0,max_len))
		for sample, l in lengths.iteritems():
			out.writerow([sample,] + l)

