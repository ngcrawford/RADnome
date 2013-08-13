import sys
import os
import time
import re

class bioFile(object):
    # fills info about file and supplies polymorphic routines for getting next hdr, next record, etc.
    def __init__(self, fileName, keepFileOpen=False):
        self.inf = [0, 0, 0] # inf[0] is type, [1] is estimate or 0 if not estimated, [2] is fh if keepFileOpen True

        self.fileName = fileName
        self.open() # open in read binary
        self.firstRecOfs = 0 # true for all 3 but SAM file with @ headers
        self.singleLineSeq = False; # True if file is fastA format and single line per sequence
        self.fileFormat = self.getBioFileFormat() # set vars associated with file type, ie fasta, fastq, or SAM

        self.nextHdr = self._hdrDefFromType(self.inf[0]) # nextHdr used by estimateRecords
        self.estRecs = self._estimateRecords()
        self.inf[1]  = self.estRecs
        self.getRecRaw = self._recRawDefFromType(self.inf[0])
       
        # close file or keep open for more operations
        if keepFileOpen and self.inf[0] > 0:
            self.inf[2] = self.fh
            self.seekToFirstRec()
        else:
            self.close()
        
    # fasta has header line with '>', then one or more sequence lines, possible blank line, this then repeats
    # fastq has header line with '@' followed by nucleotide line, then line with '+' then quality line (4 lines per record)
    # SAM file starts with 0 or more lines beginning with '@', then 1 line per record with at least 11 tabbed delimited fields
    # Given this: use very, very simple tests to see if we have one of these three (and which) or none
    def getBioFileFormat(self):
        format = 0
        ln = self.fh.readline()
        if ln:
            if ln[0] == '>': # we're just presuming fasta with this
                format = self.inf[0] = 1
                self.singleLineSeq = self.fastaSingleSeqLine() # returns True if single line for sequence in fastA
            else:
                firstRecSAM = 0
                atLines = 0
                if ln[0] == '@': # either fastq or SAM headers to skip
                    while ln and ln[0] == '@':
                        atLines += 1
                        firstRecSAM = self.fh.tell()
                        ln = self.fh.readline()
                if atLines == 1 and ln.find('\t') == -1: # could be fastq. read next line to check
                    ln = self.fh.readline()
                    if ln and ln[0] == '+': # 1st line '@', 3rd '+' we'll take it as fastq
                        format = self.inf[0] = 2
                if format == 0: # check for SAM by seeing if there are multiple tabs in the line
                    ary = ln.split('\t')
                    if len(ary) >= 11:
                        format = self.inf[0] = 3
                        self.firstRecOfs = firstRecSAM
        return format

    def open(self):
        self.fh = open(self.fileName, 'rb')
        
    def close(self):
        self.fh.close()
        self.inf[2] = 0
        self.fh = 0
        
    def seekToFirstRec(self):
        if self.fh == 0: # need to reopen file
            self.open()
        try:
            self.fh.seek(self.firstRecOfs)
        except IOError:
            self.open()
            self.fh.seek(self.firstRecOfs)

    def getNextHdr(self):
        return self.nextHdr(self.fh)

    def getRecByOfs(self, ofs): # go to offset in file and read a record
        self.fh.seek( self.indexDict[key] )
        return self.getRecord()
        
    def format(self):
        typ = self.inf[0]
        if typ == 1: return 'FASTA'
        if typ == 2: return 'FASTQ'
        if typ == 3: return 'SAM'
        return 'unknown'

    def ext(self):
        typ = self.inf[0]
        if typ == 1: return 'fa'
        if typ == 2: return 'fq'
        if typ == 3: return 'sam'
        return 'txt'

    def fastaSingleSeqLine(self): # for fasta, look through 100 records and see if these are 1 line per sequence
        self.seekToFirstRec()
        numRecs = 0
        single = True
        while numRecs < 100:
            line = self.fh.readline()
            if not line: break
            numRecs += 1
            if line[0] != '>':
                single = False # if every other line doesn't begin with '>' we don't have single line sequneces
                break
            self.fh.readline() # read putative single sequence line

        self.seekToFirstRec()
        self._lastFastaHdr = None # used with getRec routine
        return single
        
    # check 100 records at a time until at least 1 second elapsed (or file exhausted)
    def _estimateRecords(self):
        fileSize = os.stat(self.fileName).st_size
        if fileSize < 160 or self.nextHdr == None:
            return 0
            
        estimate = True # False if we read entire file
        recs = 0
        self.seekToFirstRec()
        estimateStartTime = time.time()
        while 1:
            ln, ofs = self.nextHdr(self.fh)
            if not ln:
                estimate = False
                break
            recs += 1
            if (recs % 100)==0 and (time.time() - estimateStartTime) >= 1:
                break
        if estimate:
            avgRec = self.fh.tell() / recs 
            recs = int( (fileSize+avgRec-1) / avgRec )
            
        self.seekToFirstRec()
        return recs
        
    # returns a routine that given a file handle returns next hdr and its file offset
    def _hdrDefFromType(self, type):
        if   type == 1: return self._nextFastA
        elif type == 2: return self._nextFastQ
        elif type == 3: return self._nextSAM
        else:
            return None

    @staticmethod           
    def _nextFastA(fh):
        while 1 : # first read should get '>' header line except if multi-line sequences or blank lines between records
            headerOfs = fh.tell()
            header = fh.readline()
            if not header:
                return None, headerOfs
            elif header[0] == '>':
                header = header[1:]
                break
        fh.readline()  # move past sequence line (potentially only the first) after header
        return header, headerOfs
        
    @staticmethod
    def _nextSAM(fh): # we have moved past any @ headers so each line has record header and sequence info
        headerOfs = fh.tell()
        header = fh.readline()
        return header, headerOfs
        
    @staticmethod
    def _nextFastQ(fh): # absolutely depends upon 4 lines per record
        headerOfs = fh.tell()
        header = fh.readline()
        if not header or header[0] != '@':
            return None, headerOfs
            
        fh.readline() # move past sequence line
        fh.readline() # move past '+' line that indicates quality codes on next line
        fh.readline() # move past quality code line, setting up for next header record to be read
        
        header = header[1:]
        return header, headerOfs

    # return a routine that reads the record from the current ofs, leaves ofs ready for next record to be read
    # these routines only return an array of lines read, the non-Raw routines break record into fields
    def _recRawDefFromType(self, type):
        if type == 1:
            if self.singleLineSeq:
                return self._getFastARecRaw
            else: # need to use slightly more complicated version
                return _getFastAMultiSeqRecRaw
        elif type == 2: return self._getFastQRecRaw
        elif type == 3: return self._getSAMRecRaw
        else:
            return None
            
    @staticmethod
    def _getSAMRecRaw(fh):
        lst = []
        line = fh.readline()
        if line:
            lst.append(line);
            return lst
        else:
            return None
       
    @staticmethod
    def _getFastQRecRaw(fh):
        lst = []
        line = fh.readline()
        if not line or line[0] != '@':
            return None
            
        lst.append(line)              # @ line 
        lst.append( fh.readline() );  # nucleotide line
        lst.append( fh.readline() );  # + line
        lst.append( fh.readline() );  # quality code line
        return lst
    
    @staticmethod
    def _getFastARecRaw(fh):
        lst = []
        line = fh.readline()
        if not line or line[0] != '>':
            return None
            
        lst.append(line)
        lst.append( fh.readline() )
        return lst

    def _getFastAMultiSeqRecRaw(self, fh):
        lst = []
        if self._lastFastaHdr: # last hdr read getting all lines of sequence
            line = self._lastFastaHdr
        else: # first time through very first line is header
            line = fh.readline()

        if not line or line[0] != '>':
            return None
            
        lst.append(line) # add header line to list
        
        # add sequence lines to lst until we get to next header line or finish the file
        while 1:
            line = fh.readline()
            if not line: break
                
            if line[0] != '>':
                lst.append(line)
            else: # cache header line for next call to this routine
                self._lastFastaHdr = line
                break
            
        return lst
        
    def getRecord(self):
        rec = self.getRecRaw(self.fh)
        if len(rec) > 0 and self.inf[0] < 3:
            rec[0] = rec[0][1:]
        for r in xrange(len(rec)):
            rec[r] = rec[r].rstrip()
        
        return rec
        
# end class bioFile(object)
        
# subclass of bioFile that understands Illumina style Header info.
# this is used to index the file. once indexed, the file can be merged
# with another IlluminaFile of same format (ie both fasta, fastq or SAM)
# the presumption is that the files are paired-end /1 and /2 files
# the merged file is called merged.<ext> (where <ext> is fa, fq, or sam)
# unpaired records are written into unpaired.<ext>
class IlluminaFile(bioFile):
    def __init__(self, filename):
        self.numRecs = 0; self.indexRunTime = 0; self.mergeRunTime = 0
        self.indexDict = {}; self.dupDict = {}; self.machMap = {}
        
        super(IlluminaFile, self).__init__(filename, True) # True keeps file handle open for use

    def __getitem__(self, key):
        if type(key) == int:
            return self._keyName(key) # look up ixth (i.e key as ix) key in dict and return it
        elif key not in self.indexDict:
            return None
        else:
            return self.getRecordByOfs(self.indexDict[key])

    def _keyName(self, ixKey):
        for k in self.indexDict:
            if ixKey <= 0:
                return k
            ixKey -= 1
        return None
            
    def getRecordByOfs(self, ofs):
        self.fh.seek(ofs)
        return self.getRecord()

    def getRecord(self):
        rec = bioFile.getRecord(self)
        self.header = rec[0]
        return rec
        
    def index(self):
        self.indexDict = {}; self.dupDict = {}; self.machMap = {}
        self.numRecs, self.indexRunTime = self._indexRecordsByHdr(self.indexDict, self.dupDict)
        return self.numRecs, self.indexRunTime

    def mergeWith(self, pairFileObj): # IlluminaFile object: /1 if self is /2 of paired ends, /2 file if self /1
        def isIndexed(ob):
            return len(ob.indexDict) > 0 and ob.indexRunTime > 0

        indexedFile = None; mergeFile = None   
        if not isinstance(pairFileObj, IlluminaFile):
            return 0, "Merge file must be IlluminaFile object"
        elif self.fileFormat != pairFileObj.fileFormat:
            return 0, "Can't merge " + pairFileObj.format() + " file into " + self.format() + " file."
        elif isIndexed(self):
            indexedFile = self
            mergeFile = pairFileObj
        elif isIndexed(pairFileObj):
            indexedFile = pairFileObj
            mergeFile = self
        else:
            return 0, "One of the files must be indexed"

        # open the 2 output files - one for merged record pairs, the other for unpaired records
        thisDir = os.path.split(self.fileName)[0]
        outputName = os.path.join(thisDir, "merge." + self.ext())
        fhMerge = open(outputName, "wb")
        
        unpairedName = os.path.join(thisDir, "unpaired." + self.ext())
        fhUnpaired = open(unpairedName, "wb")
        
        return self._mergeFiles(mergeFile, indexedFile, fhMerge, fhUnpaired)
        
    # encapsulates Illumina machine name and lane into small alpha code of record ID
    # and this with tile#, x-coord, y-coord becomes record ID used as index's key
    def _indexRecordsByHdr(self, indexDict, dupDict):
        machMap = self.machMap
        nextHdr = self.nextHdr # local variable fastest lookup (self.nextHdr from parent class bioFile)
        makeRecID = self.ilHdrToRecordID

        numRecs = 0L
        self.seekToFirstRec()
        fh = self.fh # if file close, seekToFirstRec will open it, cache this after it is called
        startTime = time.time()
        while 1:
            hdr, hdrOfs = nextHdr(fh) # if FastA or FastQ, hdr line returned without '>' or '@'
            if not hdr: break
            
            numRecs += 1
            key = makeRecID(hdr, machMap)

            if key not in indexDict:
                indexDict[ key ] = hdrOfs
            elif key not in dupDict: # only storing first duplicate now, could store list
                dupDict[ key ] = hdrOfs
                
        return numRecs, time.time() - startTime

    def _mergeFiles(self, mergeWith, indexedFile, fhOut, fhOut2): # wrapper class validates the arguments
        def SAMfileWithHeaders(ob): return ob.inf[0] == 3 and ob.firstRecOfs > 0
        
        keyMap  = indexedFile.indexDict # local variables fastest lexical resolution, use them inside loop
        machMap = indexedFile.machMap
        nextMergeRec = mergeWith.getRecRaw
        getIndexFileRec = indexedFile.getRecRaw # we'll position file pointer, then call this
        makeRecID = mergeWith.ilHdrToRecordID

        mergeWith.seekToFirstRec() # this will reopen file if necessary
        indexedFile.seekToFirstRec() # this will reopen file if necessary
        
        fhM = mergeWith.fh
        fhI = indexedFile.fh
        startTime = time.time()
        
        if SAMfileWithHeaders(mergeWith): # SAM file with '@' header records, write them first
            mergeWith.seek(0)
            while 1:
                line = fhM.readline()
                if not line: break
                if line[0] == '@':
                    fhOut.write(line) # write this header line into merge file
        
        numRecs = 0L; numPaired = 0L; unpairedM = 0L; unpairedI = 0L
        hdrChr = self.inf[0] < 3 # fasta and fastq files need hdr char stripped before recID made from hdr
        mergeWith.seekToFirstRec()
        
        # loop through mergeWith files records, looking up indexedFile's record with same recID
        # and writing both of them out to merge.<ext> file of fhOut
        while 1:
            rec = nextMergeRec(fhM) # return list of record lines from file to Merge With
            if not rec: break
            
            numRecs += 1
            # make recID key from hdr of mergeWith file using indexedFile's machMap so same alpha codes used
            if not hdrChr: # SAM file doesn't have '>' or '@'
                key = makeRecID(rec[0], machMap)
            else:
                key = makeRecID(rec[0][1:], machMap)
                
            if key in keyMap: # found the paired record to this one in indexedFile, write both records out
                numPaired += 2
                fhI.seek( keyMap.pop(key) ) # position to matching record in indexedFile and delete this key
                recI = getIndexFileRec(fhI)
                
                for r in xrange(len(rec)):  fhOut.write( rec[r] )
                for r in xrange(len(recI)): fhOut.write( recI[r] )
            else: # no pair with which to merge rec, write to "unpaired.<ext>" file
                unpairedM += 1
                for r in xrange(len(rec)): fhOut2.write( rec[r] )
 
        fhOut.close() # finished with merge file, flush to disk

        while keyMap: # flush unmatched indexedFile records to "unpaired.<ext>" file
            unpairedI += 1
            kvLst = keyMap.popitem() # return key value pair in no particular order
            fhI.seek( kvLst[1] )
            recI = getIndexFileRec(fhI)
            if not recI: break
            for r in xrange(len(recI)): fhOut2.write( recI[r] )
        fhOut2.close() # finished with unpaired files  
        
        mergeWith.mergeRunTime = time.time() - startTime
        indexedFile.mergeRunTime = mergeWith.mergeRunTime
        
        return (numRecs, numPaired, unpairedM, unpairedI, self.mergeRunTime)
        
    # take apart an Illumina header e.g. "HWIST640_0162:3:1101:9294:69735#TCTTGTCATC/1"
    # and turn into recordID string e.g.                "C1101:9294:69735" (C for 3rd unique hdr:lane seen)
    # updating the machMap as necessary, so same machine:lane maps to same alpha encoding
    re_digits = re.compile(r"(\d+)") # make a class variable
    @staticmethod
    def ilHdrToRecordID(hdrStr, machMap, codeMap=None):
    
        # given an Illumina machine string (preferably with lane) e.g. "EAS54_67:6",
        # if it exists in machMap, return the code, else generate new code and add it to both maps
        # codes start with upper case alphas (A-Z), then lower (a-z)
        # thus unique short character string is created for each unique machine id, lane pair
        def ilmachToCode(machStr, machMap, codeMap=None):
            if machStr in machMap:
                machCode = machMap[ machStr ]
            else:
                machCode = numToAlpha( len(machMap) ) 
                machMap[ machStr ] = machCode # map e.g. "EAS54_67:6" to "C"
                if codeMap != None:
                    codeMap[ machCode ] = machStr # reverse map e.g. "C" to "EAS54_67:6"
            return machCode
        
        # return alpha for num 0='A', 25='Z' 26='a', 52='AA' 53='AB' 103='Az' 104='BA'
        # Only 1 or 2 chars up to 2755 machine/lane pairs: 2755='zz', 2756='AAA'
        def numToAlpha(num):
            rmdr = num % 52
            if rmdr < 26:
                rmdrChr = chr( rmdr + ord('A') ) # 0 - 25 is 'A'-'Z'
            else: # rmdr mod 52 so has to be < 52
                rmdrChr = chr( rmdr-26 + ord('a') ) # 27 - 51 is 'a'-'z'

            togo = num-rmdr
            if togo <= 0:
                alpha = rmdrChr
            else:
                alpha = numToAlpha( int(togo/52)-1 ) + rmdrChr
            return alpha
    
        aParts = hdrStr.split(':') # prefix char, if any, removed beforehand
        #assert(len(aParts) > 4) # ill formatted hdr handling needs to be defined

        # encode colon separated unique machine ID and lane ID
        recordID = ilmachToCode(aParts[0]+':'+aParts[1], machMap, codeMap)

        # add tile# and x-coord (no sep before tile# since machCode is alphabetic)
        recordID += aParts[2] + ':' + aParts[3] + ':' 
        
        # add y-coord with anything after run of digits thrown away
        m = IlluminaFile.re_digits.match(aParts[4]) # n.b.: compiled r.e takes 3 sec off 2.38MM rec test
        if m: recordID += m.group()
        else: recordID += aParts[4]

        return recordID
#end class IlluminaFile(bioFile)    
        