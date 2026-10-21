"""Find mm axis mapping of an E1300 chart using OCR of tick labels (bottom mm axis, right mm axis)."""
import sys, subprocess, numpy as np
from PIL import Image
def ocr_boxes(img, scale=3):
    im = img.resize((img.width*scale, img.height*scale), Image.LANCZOS)
    im.save('/tmp/claude-0/-home-claude/fa34c777-b303-5da7-acc5-8854c649bf21/scratchpad/e1300/_o.png')
    out = subprocess.run(['tesseract','/tmp/claude-0/-home-claude/fa34c777-b303-5da7-acc5-8854c649bf21/scratchpad/e1300/_o.png','-','--psm','11','tsv'],capture_output=True,text=True).stdout
    res=[]
    for line in out.splitlines()[1:]:
        p=line.split('\t')
        if len(p)<12: continue
        txt=p[11].strip().replace(',','')
        if txt.isdigit():
            l,t,w,h=map(int,p[6:10]); res.append((int(txt),(l+w/2)/scale,(t+h/2)/scale))
    return res
def frame(fn):
    im=np.array(Image.open(fn).convert('L'))<110
    H,W=im.shape
    cols=im.sum(0); rows=im.sum(1)
    return im
def mapping(fn):
    img=Image.open(fn).convert('L'); W,H=img.size
    arr=np.array(img)<110
    # plot frame: longest vertical dark lines
    colsum=arr.sum(0); left=int(np.argmax(colsum[:W//3])); right=int(W//2+np.argmax(colsum[W//2:]))
    rowsum=arr.sum(1); bottom=int(H//2+np.argmax(rowsum[H//2:int(H*0.85)])); top=int(np.argmax(rowsum[:H//4]))
    # bottom labels: strip below bottom frame
    b=ocr_boxes(img.crop((0,bottom+3,W,bottom+30)))
    bx=[(v,x) for v,x,y in b if v%500==0]
    # right labels: strip right of frame
    r=ocr_boxes(img.crop((right+3,top-10,W,bottom+12)))
    ry=[(v,y+top-10) for v,x,y in r if v%500==0]
    return dict(left=left,right=right,top=top,bottom=bottom,bx=bx,ry=ry)
if __name__=='__main__':
    for fn in sys.argv[1:]:
        m=mapping(fn); print(fn,m)
