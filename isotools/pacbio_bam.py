from pysam import AlignmentFile,AlignedSegment,AlignmentHeader
import sys
import gzip
from tqdm import tqdm
import hashlib 
import struct
from Bio import bgzf
import logging
logger=logging.getLogger('isotools')

log.setLevel(logging.INFO)
log_format=logging.Formatter('%(levelname)s: [%(asctime)s] %(name)s: %(message)s')
#log_file=logging.FileHandler('logfile.txt')
log_stream=logging.StreamHandler()
log_stream.setFormatter(log_format)
log.handlers=[]
log.addHandler(log_stream)

type_lookup={'c':'b','C':'B', 's':'h','S':'H' }
#sam tag type characters
#A  Printable character
#i  Signed integer16
#f  Single-precision floating number
#Z  Printable string, including space
#H  Byte array in the Hex format17
#B  Integer or numeric array


base_lookup=' ACMGRSVTWYHKDBN'

#fn='/project/42/pacbio/hecatos_isoseq/01-ccs/Control_190315_A_m54070_190315_132401_ccs.bam'
class PacbioBam():
    #this class uses the pac bio bam index to assess transcripts by name
    # not optimized for high performance, but it works
    def __init__(self, fn):
        self.fn=fn
        self.is_open=False
        self.nreads=None
        self.header=None
        self.idx_fn=None
        self._open()
        self.read_header()
        self._close()
        try:
            self.read_index(fn+'.pbi')
        except FileNotFoundError:
            print('no pacbio index found, only sequential access to data')
        


    def read_index(self, idx_fn):
        # https://pacbiofileformats.readthedocs.io/en/3.0/PacBioBamIndex.html#pbi-header
        #idx_fn=self.fn+'.pbi'
        with gzip.open(idx_fn, "r") as idx:
            #header=char[4] magic, uint32 version, uint16 bpi_flags, uint32 n_reads, char[18] reserved
            magic, version, flags, nreads, _=struct.unpack('<4sIHI18s',idx.read(32)) #32 byte header
            assert magic==b'PBI\x01'
            self.nreads=nreads
            # read basic data section
            #read group is calculated as follwos:
            #id_string=hashlib.md5(bytes(movieName + "//" + readType, 'utf-8')).digest[:8] #(e.g. fa272339)
            #id_int=int.from_bytes(bytes.fromhex('fa272339'),'big')
            rgid=struct.unpack('<{}I'.format(nreads), idx.read(4*nreads)) #important if data from several runs
            start=struct.unpack('<{}I'.format(nreads), idx.read(4*nreads)) #always 0
            end=struct.unpack('<{}I'.format(nreads), idx.read(4*nreads))
            hnr=struct.unpack('<{}I'.format(nreads), idx.read(4*nreads))
            rq=struct.unpack('<{}f'.format(nreads), idx.read(4*nreads)) 
            cx=struct.unpack('<{}B'.format(nreads), idx.read(nreads)) #always 12 (or 0?)
            offset=struct.unpack('<{}Q'.format(nreads), idx.read(8*nreads)) #bgzf_tell
            log.debug(f'read pbi index: rgid={rgid[:3]}, start={start[:3]}, end={end[:3]}, hnr={hnr[:3]}, rq={rq[:3]}, cx={cx[:3]},offset={offset[:3]}')
            self.offset=dict(zip(hnr, offset))
            self.idx_fn=idx_fn
        #return nreads, rgid, dict(zip(hnr, offset))

    def read_header(self):
        #self.open()
        self._fh.seek(0)
        magic, hlen=struct.unpack('<4sI',self._fh.read(8)) 
        assert magic==b'BAM\x01'
        header=dict()
        for section in self._fh.read(hlen).decode('utf-8')[1:].split('@'):
            for elem in section.split('\n')[:-1]:
                k,*val = elem.split('\t')
                header.setdefault(k,list()).append(dict(v.split(':') for v in val))
        header['HD']=header['HD'][0]
        n_ref=int.from_bytes(self._fh.read(4), 'little')
        self.refs=list()
        for _ in range(n_ref):
            ln=int.from_bytes(self._fh.read(4), 'little')
            self.refs.append(struct.unpack('<{}sI'.format(ln),self._fh.read(4+ln)) )
        #header[]
        #self.header=header
        self.header=AlignmentHeader.from_dict(header)
        #return hstr, refs, fh.tell()


    def fetch(self, query_list=None):
        if not self.is_open:
            raise ValueError('I/O operation on a closed file')
        if query_list is None:
            yield self.read_record()
        else:
            for query_name in query_list:
                try:
                    offset=self.offset[query_name]
                except KeyError:
                    log.warning(f'no record with name {query_name} in index')
                    continue
                except TypeError:                     
                    raise ValueError('sequential access only with index')    
                _=self._fh.seek(offset)
                yield self.read_record()
        

    def read_record(self):
        offset=self._fh.tell()
        #block_size=int.from_bytes(fh.read(4), 'little')
        #data=fh.read(block_size)
        block_size, refid, pos, lrn, mapq, bam_bin, n_cigar, flag, l_seq, next_refid, next_pos, tlen=struct.unpack('<3i2B3H4i', self._fh.read(36))
        end=offset+block_size+4
        read_name=struct.unpack('<{}s'.format(lrn),self._fh.read(lrn) )
        cigar=_make_cigar(struct.unpack('<{}I'.format( n_cigar ),self._fh.read(n_cigar*4)))
        seq=_make_seq(struct.unpack('<{}B'.format( (l_seq+1)//2) ,self._fh.read((l_seq+1)//2 )))
        #log.debug(f'reading qual at {self._fh.tell()}')
        #qual=struct.unpack('<1B',self._fh.read(1 ))
        log.debug(f'now at {self._fh.tell()}')

        qual=struct.unpack('<{}B'.format( l_seq),self._fh.read(l_seq ))
        if qual[0]==255: #pacbio seems to fill empty quality strings with 255
            qual=[9] # 9=ord('*')-33 #this is how it is supposed to be 
        #todo: in order to read correctly formated bam qual, one would need to read one byte, check if it is 9 and if not read the rest.
        
        log.debug(f'read {qual}\nnow at {self._fh.tell()}')
        
        #parse the tags
        tags=list()
        # type_lookup={'c':'b','C':'B', 's':'h','S':'H' }
        # sam tag type characters
        # A  Printable character
        # i  Signed integer16
        # f  Single-precision floating number
        # Z  Printable string, including space
        # H  Byte array in the Hex format e.g ‘1AE301’ represents the byte array [0x1a, 0xe3, 0x1]
        # B  Integer or numeric array
        # Array types: cCsSiIf u_int8,int8,u_int16,int16,u_int32,int32,f
        while(self._fh.tell()<end): #next tag
            tag_name=self._fh.read(2 ).decode('utf-8')
            vtype=self._fh.read(1).decode('utf-8')
            if vtype=='B': #integer array
                vtype=self._fh.read(1).decode('utf-8')
                n=struct.unpack('<i',self._fh.read(4) )[0]
            else: 
                n=1
            if vtype in ('Z', 'H'): #byte string
                b=self._fh.read(1) #todo is there something like buffer=read_until(b'\00')
                tagval=bytearray()
                while b != b'\00': #todo: handle H to be 2 byte/char (what encoding? is it relevant?)
                    tagval+=(b)                
                    b=self._fh.read(1)
                val=tagval.decode('utf-8')
            else:
                if vtype in type_lookup:
                    vtype=type_lookup[vtype]
                fstring='<{}{}'.format(n,vtype)
                val=struct.unpack(fstring,self._fh.read(struct.calcsize(fstring)) )
                #if n==1:
                #    val=val[0]
                val=','.join(str(v) for v in val)#conversion list to string
                #tags[tag_name]=val
            tags.append(f'{tag_name}:{vtype}:{val}')
            #tags.append((tag_name,val))
        read={'name':read_name[0][:-1].decode('utf-8'), 'flag':str(flag), "ref_name":'*' if refid<0 else self.refs[refid], "ref_pos":str(pos+1), "map_quality":str(mapq), "cigar":cigar,
             "next_ref_name":'*' if next_refid<0 else self.refs[next_refid], "next_ref_pos":str(next_pos+1), "length":str(tlen), "seq":seq, "qual":''.join((chr(c+33) for c in qual)), "tags":tags}
        #return read
        try:
            return AlignedSegment.from_dict(read,self.header) 
        except TypeError:
            log.warning('Caught type error when transforming the alignment record')
            return read
        #return(read_name[0][:-1].decode('utf-8'), refid, pos, cigar, seq, qual, tags)

    

    def _open(self):
        if self.is_open:            
            self._close()
        idx_fn=self.fn+'.pbi'
        if idx_fn != self.idx_fn:
            self.read_index(idx_fn)

        self._fh=bgzf.open(self.fn,'rb')
        self.is_open=True
        self.read_header()
        return self
        
        
    def _close(self):
        if self.is_open:
            self._fh.close()
            self.is_open=False

    @classmethod
    def open(cls, fn, use_pbi=True):
        if use_pbi:
            try:
                align=cls(fn)._open()
            except FileNotFoundError: #cannot find pbindex, fall back to AlignmentFile
                pass
            else:                
                return(align)
        return AlignmentFile(fn, "rb")

    def __exit__(self, *args):
        self._close()
        return True
    
    def __enter__(self):
        self._open()
        return self
 

def _make_cigar( data):
    #todo: implement according to bam definition
    # this dummy returns '*' (e.g. N/A)
    return '*'
def _make_seq(data):        
        seq="".join(c for  b in data for c in (base_lookup[ b >> 4], base_lookup[b & 0x0F] ) )
        return seq.rstrip()



class Tabix_file:
    #just a proof of principle implementation of the tabix index (currently not usefull)
    #the intention is to find out the number of records for EAT estimations
    #https://samtools.github.io/hts-specs/tabix.pdf
    #apperantly it does not hold the number of records, nor means to compute it ... :(
    # the only option would be to make use of the file size...
    def __init__(self, fn):
        self.fn=fn

    def read_index(self):
        idx_fn=self.fn+'.tbi'
        with gzip.open(idx_fn, "r") as idx:
            magic, n_ref, format_id, col_seq, col_bg, col_end,cmt_char, skip, l_nm=struct.unpack('<4s8i',idx.read(4*9))
            assert magic==b'TBI\x01', 'not a tabix index'
            self.refnames=[n.decode('utf-8') for n in struct.unpack('<{}s'.format(l_nm),idx.read(l_nm))[0].split(b'\x00')]
            total=0
            for ref_idx in range(n_ref):
                n_bin=struct.unpack('<i', idx.read(4))[0]
                for bin_nr in range(n_bin):
                    bin_id, n_chunk=struct.unpack('<Ii', idx.read(8))
                    total+=n_chunk
                    for chunk_nr in range(n_chunk):
                        cnk_start, cnk_end=struct.unpack('<2Q', idx.read(16))
                nintv,=struct.unpack('<i', idx.read(4))
                ioff=struct.unpack('<{}Q'.format(nintv), idx.read(8*nintv))
            try:
                no_coor=struct.unpack('<Q', idx.read(8))
            except struct.error: #at eof: unpack requires a buffer of 8 bytes
                no_coor=0
    



if __name__=='__main__':
    #demo how to use: 
    nargs = len(sys.argv) - 1
    if nargs>0:
        #1) open bamfile and index
        pbbam=PacbioBam(sys.argv[1])
    else:   
        raise ValueError('please specify bam file')
    if nargs>1:
        #2a) access specific transcript
        print('\t'.join(pbbam.fetch(sys.argv[2]) ))
    else:
        #2b) access next transcript
        #print('\t'.join(pbbam.fetch()))
        print(pbbam.nreads)

