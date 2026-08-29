"""Canberra CNF (CAM) reader — adapted from becquerel/xylib with stricter
block validation (header entry must point to a block whose first bytes match
the section id + 0x20)."""
import struct, datetime
import numpy as np

def _le(b,i,n): return int.from_bytes(b[i:i+n],'little')

def _date(b,i):
    d=_le(b,i,8); t=(d/1e7)-3506716800
    return datetime.datetime.utcfromtimestamp(t)

def _time(b,i):
    d=_le(b,i,8); d=(pow(2,64)-1) & ~d
    return d*1e-7

def _pdp11(b,i):
    sign=-1 if (b[i+1]&0x80) else 1
    exb=((b[i+1]&0x7F)<<1)+((b[i]&0x80)>>7)
    if exb==0: return 0.0 if sign==1 else float('nan')
    h=b[i+2]/256.**3 + b[i+3]/256.**2 + (128+(b[i]&0x7F))/256.
    return sign*h*2.**(exb-128)

def _ecal(b,i):
    c=[_pdp11(b,i+2*4+28+4*k) for k in range(4)]
    return None if c[1]==0.0 else c

def read_cnf(fn):
    b=open(fn,'rb').read()
    offs={0:0,1:0,2:0,5:0}
    for i in range(112,min(len(b),128*1024),48):
        off=_le(b,i+10,4)
        b0,b1,b2=b[i],b[i+1],b[i+2]
        if not ((b1==0x20 and b2==0x01) or b1==0 or b2==0): continue
        if b0 in offs and offs[b0]==0 and 0<off<len(b):
            # validate block start
            if off+1<len(b) and b[off]==b0 and b[off+1]==0x20:
                offs[b0]=off
    oa,os_,oc5=offs[0],offs[1],offs[5]
    if not oa or not oc5: raise ValueError('sections not found '+str(offs))
    out={}
    if os_:
        g=lambda a,n: b[os_+a:os_+a+n].decode('latin-1').strip('\x00').strip()
        out['sample_name']=g(48,64); out['sample_id']=g(112,64)
        out['sample_desc']=g(0x036E,256)
    off1=_le(b,oa+34,2); off2=_le(b,oa+36,2)
    opha=oa+48+128
    nch=256*_le(b,opha+10,2)
    odate=oa+48+off2+1
    out['start']=_date(b,odate)
    out['realtime']=_time(b,odate+8)
    out['livetime']=_time(b,odate+16)
    ocal=oa+48+32+off1
    cc=_ecal(b,ocal) or _ecal(b,ocal-off1)
    out['ecal']=cc
    cnt=np.frombuffer(b,dtype='<u4',count=nch,offset=oc5+512).astype(float).copy()
    for i in range(2):
        if int(cnt[i])==int(out['realtime']) or int(cnt[i])==int(out['livetime']): cnt[i]=0
    out['counts']=cnt; out['nch']=nch
    return out

def energy(ch, cc):
    return cc[0]+cc[1]*ch+cc[2]*ch**2+cc[3]*ch**3
