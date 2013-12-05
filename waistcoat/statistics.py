import csv, os, os.path

recording = True

#data
names = []
data = {}

def clear():
	"""clear recorded statistics"""
	data.clear()
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


def prettyString():
	"""return a pretty representation of the statistics"""
	col_widths = []

	#find the maximum width of each column
	col_widths.append(max(len(str(x)) for x in ['sample',] + data.keys()))
	for i,name in enumerate(names):
		col_widths.append(max(len(str(x)) for x in 
				[names[i],] + [s[i] for s in data.values()]))

	hdr = ' '.join(["{{:<{}s}}".format(col_widths[0]),] + 
			["{{:>{}s}}".format(w) for w in col_widths[1:]])
	line = ' '.join(["{{:<{}s}}".format(col_widths[0]),] + 
			["{{:>{}d}}".format(w) for w in col_widths[1:]])

	ret = [hdr.format('sample', *names),]
	for sample, values in data.iteritems():
		ret.append(line.format(sample, *values))

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



