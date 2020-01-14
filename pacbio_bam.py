from pysam import AlignmentFile
import sys
import gzip
from tqdm import tqdm
import hashlib 
import struct
from Bio import bgzf

#fn='/project/42/pacbio/hecatos_isoseq/01-ccs/Control_190315_A_m54070_190315_132401_ccs.bam'
class PacbioBam():
    base_lookup='=ACMGRSVTWYHKDBN'
    type_lookup={'c':'b','C':'B', 's':'h','S':'H' }
    def __init__(self, fn):
        self.fn=fn
        self.is_open=False
        self.nreads=None
        self.header=None
        self.read_header()
        try:
            self.read_index()
        except FileNotFoundError:
            print('no pacbio index found, only sequential access to data')
        


    def read_index(self):
        # https://pacbiofileformats.readthedocs.io/en/3.0/PacBioBamIndex.html#pbi-header
        fn=self.fn+'.pbi'
        with gzip.open(fn, "r") as idx:
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
            cx=struct.unpack('<{}B'.format(nreads), idx.read(nreads)) #always 12
            offset=struct.unpack('<{}Q'.format(nreads), idx.read(8*nreads)) #bgzf_tell
            self.offset=dict(zip(hnr, offset))
        #return nreads, rgid, dict(zip(hnr, offset))

    def read_header(self):
        self.open()
        self._fh.seek(0)
        magic, hlen=struct.unpack('<4sI',self._fh.read(8)) 
        assert magic==b'BAM\x01'
        self.header=self._fh.read(hlen).decode('utf-8')
        n_ref=int.from_bytes(self._fh.read(4), 'little')
        self.refs=list()
        for i in range(n_ref):
            ln=int.from_bytes(self._fh.read(4), 'little')
            self.refs.append(struct.unpack('<{}sI'.format(ln),self._fh.read(4+ln)) )
        #return hstr, refs, fh.tell()

    def fetch_zmw(self, zmw=None):
        self.open()
        if zmw is not None:
            try:
                offset=self.offset[zmw]
            except KeyError:
                raise ValueError('ZMW {} not in index'.format(zmw))
            except TypeError:
                raise ValueError('no index')    
            _=self._fh.seek(offset)
        else: 
            #try to read from where we are
            offset=self._fh.tell()

        #block_size=int.from_bytes(fh.read(4), 'little')
        #data=fh.read(block_size)
        block_size, refid, pos, lrn, mapq, bam_bin, n_cigar, flag, l_seq, next_refid, next_pos, tlen=struct.unpack('<3i2B3H4i', self._fh.read(36))
        end=offset+block_size+4
        read_name=struct.unpack('<{}s'.format(lrn),self._fh.read(lrn) )
        cigar=self._make_cigar(struct.unpack('<{}I'.format( n_cigar ),self._fh.read(n_cigar*4)))
        seq=self._make_seq(struct.unpack('<{}B'.format( (l_seq+1)//2) ,self._fh.read((l_seq+1)//2 )))
        qual=struct.unpack('<{}B'.format( l_seq),self._fh.read(l_seq ))
        tags=dict()
        #fh.seek(7014)
        while(self._fh.tell()<end): #read the tags
            tag_name=self._fh.read(2 ).decode('utf-8')
            vtype=self._fh.read(1).decode('utf-8')
            if vtype=='B':
                vtype=self._fh.read(1).decode('utf-8')
                n=struct.unpack('<i',self._fh.read(4) )[0]
            else: 
                n=1
            if vtype in ('Z', 'H'): #array
                b=self._fh.read(1)
                tagval=bytearray()
                while b != b'\00': #todo: handle H to be 2 byte/char (what encoding?)
                    tagval+=(b)                
                    b=self._fh.read(1)
                tags[tag_name]=tagval.decode('utf-8')
            else:
                if vtype in self.type_lookup:
                    vtype=self.type_lookup[vtype]
                fstring='<{}{}'.format(n,vtype)
                val=struct.unpack(fstring,self._fh.read(struct.calcsize(fstring)) )
                if n==1:
                    val=val[0]
                tags[tag_name]=val
        return(read_name[0][:-1].decode('utf-8'), refid, pos, cigar, seq, qual, tags)


    def _make_cigar(self, data):
        #todo: implement according to bam definition
        return data

    def _make_seq(self, data):        
        return("".join(c for  b in data for c in (self.base_lookup[ b >> 4], self.base_lookup[b & 0x0F] )))


    def open(self):
        if not self.is_open:
            self._fh=bgzf.open(self.fn,'rb')
            self.is_open=True
        
    def close(self):
        if self.is_open:
            self._fh.close()
            self.is_open=False

  
class Tabix_file:
    #just a proof of principle implementation of the tabix index
    #https://samtools.github.io/hts-specs/tabix.pdf
    #apperantly it does not hold the number of records, nor means to compute it ... :(
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
        #2a) access specific ZMW
        print('\t'.join(pbbam.fetch_zmw(sys.argv[2]) ))
    else:
        #2b) access next ZMW
        print('\t'.join(pbbam.fetch_zmw()))

