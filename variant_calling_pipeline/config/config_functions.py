
# These functions parse the file names and determine if the files are samples or blanks

# These are the only functions that may need to be changed for each run
def parse_samples(fname):
    # Change the number of split elements returned depending on the sample name
    return '_'.join(fname.split('_')[0:2])

def parse_blanks(fname):
    # Change the number of split elements returned depending on the sample name
    return '_'.join(fname.split('_')[0:2])


# These functions don't need to be changed.
def parse_sampleID(fname, prefix_list):
    if parse_samples(fname).startswith(tuple(prefix_list)):
        return parse_samples(fname)
    else:
        return parse_blanks(fname)




