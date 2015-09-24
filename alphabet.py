"""
This module was created by Mikael Boden.
"""

class Alphabet(object):
    """ Defines an immutable biological alphabet (e.g. the alphabet for DNA is AGCT) 
    that can be used to create sequences (see sequence.py).
    We use alphabets to define "tuple" tables, where entries are keyed by combinations
    of symbols of an alphabet (see class TupleStore below). 
    Alphabets are used to define probability distributions for stochastic events
    (see prob.py). """
    
    def __init__(self, symbolString):
        """ Construct an alphabet from a string of symbols. Lower case characters 
        will be converted to upper case, repeated characters are ignored.
        Example of constructing the DNA alphabet:
        >>> alpha = Alphabet('ACGTttga')
        >>> alpha.symbols
        ('A', 'C', 'G', 'T') """
        # Add each symbol to the symbols list, one at a time, and ignore doubles (could use "set" here...)
        _symbols = [] # create a temporary list
        for s in symbolString:
            if not str(s).upper()[0] in _symbols:
                _symbols.append(str(s).upper()[0])
        _symbols.sort() # we put them in alphabetical (one canonical) order
        # OK done extracting, put them in place
        self.symbols = tuple(_symbols); # create the immutable tuple from the extracted list
        self.length = len(self.symbols)
        self.annotations = {}

    def __str__(self):
        return str(self.symbols)
    
    def __len__(self):
        return len(self.symbols)
    
    def __iter__(self):
        return self.symbols.__iter__()
    
    def __getitem__(self, ndx):
        """ Retrieve the symbol(s) at the specified index (or slice of indices) """
        return self.symbols[ndx]
    
    def __contains__(self, sym):
        """ Check if the given symbol is a member of the alphabet. """
        return sym in self.symbols
    
    def index(self, sym):
        """ Retrieve the index of the given symbol in the alphabet. """
        # If the symbol is valid, use the tuple's index function
        if sym in self.symbols:
            syms = self.symbols
            return syms.index(sym)
        else:
            raise RuntimeError('Symbol %s is not indexed by alphabet %s' % (sym, str(self.symbols)))
        
    def __eq__(self, rhs):
        """ Test if the rhs alphabet is equal to ours. """
        if rhs == None:
            return False
        if len(rhs) != len(self):
            return False
        # OK we know they're same size...
        for sym in self.symbols:
            if not sym in rhs:
                return False
        return True

    def isSubsetOf(self, alpha2):
        """ Test if this alphabet is a subset of alpha2. """
        for sym in self.symbols:
            if not alpha2.isValidSymbol(sym):
                return False
        return True
    
    def isSupersetOf(self, alpha2):
        """ Test if this alphabet is a superset of alpha2. """
        return alpha2.isSubsetOf(self)
    
    def annotateSym(self, label, sym, value):
        try:
            lookup = self.annotations[label]
        except KeyError:
            lookup = self.annotations[label] = {}
        lookup[sym] = value
            
    def annotateAll(self, label, symdictOrFilename):
        if isinstance(symdictOrFilename, str): # we assume it is a filename
            fh = open(symdictOrFilename)
            string = fh.read()
            d = {}
            for line in string.splitlines():
                if len(line.strip()) == 0:
                    continue
                sections = line.split()
                symstr, value = sections[0:2]
                for sym in symstr:
                    d[sym] = value
            fh.close()
        else: # we assume it is a dictionary 
            d = symdictOrFilename
        for sym in d:
            self.annotateSym(label, sym, d[sym])
        
    def getAnnotation(self, label, sym):
        try:
            lookup = self.annotations[label]
            return lookup[sym]
        except KeyError:
            return None


""" Below we declare alphabets that are going to be available when 
this module is imported """
Bool_Alphabet = Alphabet('TF')
DNA_Alphabet = Alphabet('ACGT')
DNA_Alphabet_wN = Alphabet('ACGTN')
RNA_Alphabet = Alphabet('ACGU')
Protein_Alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
Protein_Alphabet_wX = Protein_wX = Alphabet('ACDEFGHIKLMNPQRSTVWYX')
Protein_Alphabet_wSTOP = Alphabet('ACDEFGHIKLMNPQRSTVWY*')
DSSP_Alphabet = Alphabet('GHITEBSC')
DSSP3_Alphabet = Alphabet('HEC')

predefAlphabets = {'DNA': DNA_Alphabet,
                   'RNA': RNA_Alphabet,
                   'DNAwN': Alphabet('ACGTN'),
                   'RNAwN': Alphabet('ACGUN'),
                   'Protein': Protein_Alphabet,
                   'ProteinwX': Protein_wX}
# The preferred order in which a predefined alphabet is assigned to a sequence 
# (e.g., we'd want to assign DNA to 'AGCT', even though Protein is also valid)
preferredOrder = ['DNA', 'RNA', 'DNAwN', 'RNAwN', 'Protein', 'ProteinwX']
# Useful annotations
DNA_Alphabet.annotateAll('html-color', {'A':'green','C':'orange','G':'red','T':'#66bbff'})
RNA_Alphabet.annotateAll('html-color', {'A':'green','C':'orange','G':'red','U':'#66bbff'})
Protein_Alphabet.annotateAll('html-color', {'G':'orange','P':'orange','S':'orange','T':'orange','H':'red','K':'red','R':'red','F':'#66bbff','Y':'#66bbff','W':'#66bbff','I':'green','L':'green','M':'green','V':'green'})


