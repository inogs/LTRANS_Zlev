import math 
import numpy as np
import math
import os,sys
import glob
from collections import Counter

def is_non_zero_file(fpath):  
     return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def runcommand(command):
        print(command)
        os.system(command)

def strmask(number,string,optstring='.'):
        if number==0 : outstring=optstring
        else:
          if(number<1):
            outstring='+'
          else: 
            outstring=string
        return outstring


def printval(var,mul=1.0,leadingzeros=0):
        (nz,ny,nx)=getnznynx(var) 
        var=var*mul
        if(leadingzeros==0):space=0
        else: space=1
        lenitem=len(str(int(np.max(var))))+leadingzeros+space
        for k in range(0,max(1,nz)):

          print('-'*nx*(lenitem+2))
          for j in range(ny-1,-1,-1):
              for i in range(0,nx): 
                print(emptyfillednumberstring(cropfloat(out[j,i],leadingzeros),lenitem+1),)
              print(' ')
        print('-'*nx*(lenitem+2))
        print(' ')

def printmask(m_rho,m_u=None,m_v=None):
        (nz,ny,nx)=getnznynx(m_rho)
        nxu=nx-1
        if(m_u != None):(nzu,nyu,nxu)=getnznynx(m_u)
        if(m_v != None):(nzv,nyv,nxv)=getnznynx(m_v)
        print('nz=',nz)
        if(nz==0):
                print('----'*nx)
                for j in range(ny-1,-1,-1):
                  if(m_v != None):  
                    if(j<nyv): 
                      for i in range(0,nxv): 
                        print(strmask(m_v[j,i],'v'),' ',)
                      print(' ')
                  for i in range(0,nx):
                      if(i<nxu): 
                        print(strmask(m_rho[j,i],'X'),)
                        if(m_u != None):  strmask(m_u[j,i],'u'),
                      else:       print(strmask(m_rho[j,i],'X','.'))
        else:
                for k in range(0,nz):
                  print('----'*int((nx-1)/2),)
                  print(' '+str(k)+' ',)
                  print('----'*int(nx-(nx-1)/2-1))
                  for j in range(ny-1,-1,-1):
                    if(j<nyv): 
                      for i in range(0,nxv): 
                        print(strmask(m_v[k,j,i],'v'),' ',)
                      print(' ')
                    for i in range(0,nx): 
                      if(i<nxu):  print(strmask(m_rho[k,j,i],'X'),strmask(m_u[k,j,i],'u'),)
                      else:       print(strmask(m_rho[k,j,i],'X'))
        print('----'*nx)
        print(' '       )

def getnznynx(ARRAY):
   if(len(np.shape(ARRAY))==2):
        nz=0
        ny=np.shape(ARRAY)[0]
        nx=np.shape(ARRAY)[1]
   else:
        nz=np.shape(ARRAY)[0]
        ny=np.shape(ARRAY)[1]
        nx=np.shape(ARRAY)[2]
   return (nz,ny,nx)


def getgridparamsfromSTDOUT(directory,identifier,dirout='',
           f_Eta='Eta.',f_RHOA='RHOAnoma.',
           f_U='U.',f_V='V.',f_W='W.',f_S='S.',f_T='T.',
           f_KPPdiffS='KPPdiffS.',f_EXFuwind='EXFuwind.',f_EXFvwind='EXFvwind.',f_EXFIwind='EXFiwind.',plotdir='./',fileextension=''):

        FilesNamesContain=[f_Eta,f_RHOA,f_U,f_V,f_W,f_S,f_T,f_KPPdiffS,f_EXFuwind,f_EXFvwind,f_EXFIwind]
        Freqcy_List=np.zeros((1000),dtype=float)
        VarNam_List=np.empty((1000),dtype="S15")
        VarNam_List[:]=''
        try:
                LIST = glob.glob(directory+'/STDOUT*')
        except:
                LIST = glob.glob(directory)
        LIST.sort()
        if(len(LIST)==0):
          print('No STDOUT file found in ',directory)
          resultsdirectory=directory.replace('run_','results_')
          try: LIST = glob.glob(resultsdirectory+'/STDOUT*')
          except: LIST = glob.glob(directory)
          if(len(LIST)==0):
            print('No STDOUT file found in ',resultsdirectory)
        elif(len(LIST)!=1):
          print('found more than 1 STDOUT file in ',directory,':',LIST)
          print('trying with ',LIST[0])
        if(True):       
          if(is_non_zero_file(plotdir+'/setup_'+identifier+'.py')):
              print('File '+plotdir+'/setup_'+identifier+'.py ALREADY EXISTS you first have to delete it if you want to proceed')
              quit()
          os.system('mkdir '+plotdir)
          fw=open(plotdir+'/setup_'+identifier+'.py','w')
          fw.write("dirin='"+directory+"' \n")
          fw.write("identifier='"+identifier+"' \n")
          file=LIST[0]
          f = open(file, 'r')
          filecontent = f.readlines()
          waitforend=False
          waitforfrequency=False
          waitingforprecision=False
          delXYZ=''
          PRECISION=''
          nXYZ=''
          countline=0
          words=['','']
          nXYZ_dimsFacets=''
          for line in filecontent:
             countline+=1
             if(countline>5000):break
             oldwords=words
             words = line.split()
             #foundarrow=False
             for i in range(1,len(words)):
               #if(foundarrow or words[i-1]=='>'):foundarrow=True
               #else:continue
               word=words[i]
               if(word[0]=='#' or word[min(len(word)-1,1)]=='#' or word[min(len(word)-1,2)]=='#'): break         
               if(waitforend):
                  #print('waitforend at line ',countline,' among::: ',line,' ; checking ',word)
                  if(word[1:4]=='END'):
                        delXYZ=delXYZ[:-1]+'] \n'
                        #print(delXYZ)
                        waitforend=False
                        break
                  elif(words[i-1]=='>'): 
                        #print('++++ ',''.join(words[i:]))
                        string=''.join(words[i:])
                        delXYZ=delXYZ+string 
                        break
               elif(waitforfrequency):
                   if(word[0:9]=='frequency'):
                        
                        for p1,char in enumerate(word[9:]):
                           if char=='(':
                              for p2,charend in enumerate(word[9+p1+1:]):
                                   if charend==')':
                                       num=int(word[9+p1+1:9+p1+1+p2])
                              Freqcy_List[num]=words[i+2][:-1]
                   if(word[0:8]=='fileName' or word[0:8]=='filename' or word[0:8]=='FileName'):
                        for p1,char in enumerate(word[8:]):
                           if char=='(':
                               for p2,charend in enumerate(word[8+p1+1:]):
                                   if charend==')':
                                       num=int(word[8+p1+1:8+p1+1+p2])
                               VarNam_List[num]=words[i+2][:-1]
                   if(word[1:4]=='END' or word[1:5]=='&END'):
                        for i in range(1,len(VarNam_List)):
                          if(VarNam_List[i]==''):break
                          #print('output frequency(',VarNam_List[i],')=',Freqcy_List[i])
                        for num,n1 in enumerate(VarNam_List):
                          name1=n1.decode('ascii')
                          for char in name1:
                             if char in "',.": name1=name1.replace(char,"")
                          for name2 in FilesNamesContain:
                             for char in name2:
                               if char in "',.": name2=name2.replace(char,'')
                             if(name1==name2):
                                 FreqOutput=Freqcy_List[num]
                                 waitforfrequency=False
                                 print('taking as output frequency the one of ',name1,'=',FreqOutput)
                                 break
               else:
                  if(word[0:5]=='delX=' or word[0:5]=='delY='):
                        string=(''.join(words[i:]))[:-1]
                        delXYZ=delXYZ+string+' \n'
                        break
                  elif(word[0:5]=='delZ='):
                        string=words[i]+'['+(''.join(words[i+1:]))
                        delXYZ=delXYZ+string
                        waitforend=True
                        #print(string,' :: RUN')
                        break
                  elif(word[0:15]=='writeBinaryPrec'):
                        string=(''.join(words[i:]))[:-1]
                        print('searching for PRECISION ',PRECISION)
                        if(string[16:17]!='/*' and string[16]!='#' and string[16:17]!='//'):
                          try:
                            intprec=int(string[16:-1])
                            PRECISION="writeBinaryPrec="+str(intprec)+' \n'
                          except:
                            waitingforprecision=True
                        else:
                            waitingforprecision=True
                        break
                  elif(waitingforprecision):
                        print(word)
                        print(words)
                        try:
                          intprec=int(word)
                          PRECISION="writeBinaryPrec="+str(intprec)+' \n'
                          waitingforprecision=False
                        except:
                          continue
                  elif(word=='Nx' or word=='Ny'):
                        string=(''.join(words[i:i+3]))
                        string='n'+string[1]+'rho_in'+string[2:] #  'Nx -> nxrho_in'
                        nXYZ=nXYZ+string+' \n'
                        #print('nXYZ=',nXYZ)
                        break
                  elif(word=='dimsFacets'):
                        nx_dimsFacets=words[i+2:i+3][0][:-1]
                        ny_dimsFacets=words[i+3:i+4][0][:-1]
                        nXYZ_dimsFacets='nxrho_in='+nx_dimsFacets+' \nnyrho_in='+ny_dimsFacets+' \n'
                        break
                  elif(word=='&DIAGNOSTICS_LIST'):
#                  elif(len(words)>=i+3):
#                      if(words[i][0:min(len(words[i]),8)]=='fileName' and words[i+2][0:min(len(words[i+2]),5)]=="'Eta'"):
                        waitforfrequency=True
                        print('lets wait for frequency, words=',words)
                        linestartwait=countline
                        break
          if(len(nXYZ_dimsFacets)>0) : nXYZ=nXYZ_dimsFacets
          print('----------')
          fw.write(delXYZ)
          print(delXYZ)
          fw.write(PRECISION)
          print(PRECISION)
          nXYZ=nXYZ+'nzrho_in=len(delZ) \n'+'nzrho_in_cut=nzrho_in \n'
          fw.write(nXYZ)
          print(nXYZ)
          fw.write('OUTPUTfrequency = '+str(FreqOutput)+' \n')
          GridFilesNamesContain=['XC','YC','hFacC','XG','YG']
          GridFilesNamesDOESNOTContain=['DXC','DYC','DXG','DYG']
          EXCLUDELIST=[]
          for name in GridFilesNamesDOESNOTContain:
            LIST = glob.glob(directory+'*'+name+'*')
            LIST.sort()
            EXCLUDELIST=EXCLUDELIST+LIST
          g_filenames="g_filenames =["
          for name in GridFilesNamesContain:
            LIST = glob.glob(directory+'*'+name+'*')
            LIST.sort()
            if(len(LIST)==1):
                g_filenames=g_filenames+"'"+LIST[0][len(directory):]+"',"
                print('for file ',name,' kept only ',LIST[0][len(directory):])
            else:
                print('for file ',name,' kept :',)
                for n in range(0,len(LIST)):
                  found=False
                  for exclude in EXCLUDELIST:
                    if(LIST[n]==exclude):found=True
                  if(not found):
                    g_filenames=g_filenames+"'"+LIST[n][len(directory):]+"',"
                    print(LIST[0][len(directory):],)
                    break
                print(' DONE.')
                #if(not found):print('kept ',LIST[n][len(directory):],' among ',LIST,' searching for ',directory+'*'+name+'*')
                #else:
                #   print('error ',LIST,' searching for ',directory+'*'+name+'*')
                #   quit()
          GridFilesNamesContain=['atimetry','epth*.data','athimetry','athymetry']
          found=False
          for name in GridFilesNamesContain:
            LIST = glob.glob(directory+'*'+name+'*')
            LIST.sort()
            print(name,LIST)
            if(len(LIST)==1):
                if(found):print('ERROR found 2 options for bathymetry/depth files')
                g_filenames=g_filenames+"'"+str(LIST[0][len(directory):])+"'] \n"
                found=True
          g_filenames=g_filenames+'nfXC   =         0  \n'
          g_filenames=g_filenames+'nfYC   =         1  \n'
          g_filenames=g_filenames+'nfhFacC=         2  \n'
          g_filenames=g_filenames+'nfXG   =         3  \n'
          g_filenames=g_filenames+'nfYG   =         4  \n'
          g_filenames=g_filenames+'nfDepth=         5  \n'
          print(g_filenames)
          fw.write(g_filenames)
          fw.write('numgridfiles=6 \n')
          fw.write('reslon=float(nxrho_in)/delX \n')
          fw.write('reslat=float(nyrho_in)/delY \n')
          fw.write('f_newvarname=[ "zeta","rho", "u", "v" ,"w", "salt","temp", "AKs", "uwind", "vwind"] \n')
          maxnumfiles=11
          FoundFile=[False]*maxnumfiles
          FilesNames=['']*maxnumfiles
          FilesPrefix=['']*maxnumfiles
          numfieldfiles=0
          f_filenames='f_filenames = '
          count=0
          found=0
          print(FilesNamesContain)

          for numf,name in enumerate(FilesNamesContain):
            LIST = glob.glob(directory+name+'*')
            LIST.sort()
            if(len(LIST)==1):
                item=0
                keep=0
            else:
                length=100
                keep=-1
                for i,item in enumerate(LIST):
                   if(len(item)<length):
                      keep=i
                      length=len(item)
                print('WARNING found ',len(LIST),' files for variable',name)
                if(keep>=0):
                  #for i,item in enumerate(LIST):print('     -',item[len(directory):])
                  print('selected file is:',LIST[keep][len(directory):])
            if(keep>=0):
                fname=LIST[keep][len(directory):]
                NUMDIGITS=findOccurencesnumber(fname,['0','1','2','3','4','5','6','7','8','9'])
                extension=fname[reversefindOccurences(fname,'.'):]
                if(NUMDIGITS==0):
                        FILENUM=0
                        FILESTEP=0
                        TDIM=0
                        NUMDAYS=1.0
                else:
                        TDIM=1
                        FILENUM=int(fname[-len(extension)-NUMDIGITS:-len(extension)])
                        FILESTEP=9999999999
                        NUMFILES=0
                        for k in range(keep+1,len(LIST)):
                            otherfname=LIST[k][len(directory):]
                            if(findOccurencesnumber(otherfname,['0','1','2','3','4','5','6','7','8','9'])==NUMDIGITS 
                                and otherfname[reversefindOccurences(otherfname,'.'):]==extension):
                              NUMFILES=NUMFILES+1
                              FILEDIFF=int(otherfname[-len(extension)-NUMDIGITS:-len(extension)])-FILENUM
                              if(FILEDIFF<FILESTEP):FILESTEP=FILEDIFF
                        if(FILESTEP==9999999999):FILESTEP=0
                        NUMDAYS=int(NUMFILES*abs(float(FreqOutput)/86400.0))
                FilesNames[count]=directory+fname
                if(NUMDIGITS>0):ONEFORDOT=1
                else:ONEFORDOT=0
                fname=fname[:-len(extension)-NUMDIGITS+1-ONEFORDOT]
                FilesPrefix[count]=fname
                print(FilesNames[count],FilesPrefix[count],fname,extension,NUMDIGITS)
                f_filenames=f_filenames+"'"+fname+"',"
                found=found+1
                numfieldfiles=numfieldfiles+1
                FoundFile[count]=True
            count=count+1
          for count in range(0,maxnumfiles):
                if(not FoundFile[count]):
                     fname=str(FilesNamesContain[count])
                     print('prefix ',fname,' added to LTRANS data file even if file is not present in hydro directory')
                     f_filenames=f_filenames+"'"+fname+"',"
          f_filenames=f_filenames[:-1]+' \n'
          print(f_filenames)
          fw.write(f_filenames)
          print('f_filesextension="'+extension+'"')
          fw.write('f_filesextension="'+extension+'"'+'\n')
          fw.write('hydrofileextension="'+extension+'"'+'\n')
          if(FoundFile[2]):fw.write('hydroUfile="'+FilesNames[2]+'" \n')
          if(FoundFile[3]):fw.write('hydroVfile="'+FilesNames[3]+'" \n')
          if(FoundFile[8]):fw.write('windUfile="'+FilesNames[8]+'" \n')
          if(FoundFile[9]):fw.write('windVfile="'+FilesNames[9]+'" \n')
          fw.write('freqforcing='+str(int(float(FreqOutput)))+' \n')
          fw.write('recordnum=1 \n')
          print('numfieldfiles='+str(numfieldfiles))
          fw.write('numfieldfiles='+str(numfieldfiles)+'\n')
          fw.write('numzerosinfieldfile=0 \n')
          fw.write('fieldfilesmergedintime=False \n')
          fw.write('timerange=[0]   # uncomment if multiples files \n')
          fw.close()
          #------------------------------------
          if(is_non_zero_file(dirout+'LTRANS_'+identifier+'.data')):
              print('File '+dirout+'LTRANS_'+identifier+'.data ALREADY EXISTS')
              print('you may delete it AS WELL AS '+plotdir+'/setup_'+identifier+'.py if you want to proceed')
              quit()
          else:
              print('creating file '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('cp LTRANS_model_file.data '+dirout+'LTRANS_'+identifier+'.data')
          print('import '+plotdir+'/setup_'+identifier+'.py')
          filename=plotdir+'/setup_'+identifier+'.py'
          with open(filename, "rb") as source_file:
            code = compile(source_file.read(), filename, "exec")
          exec(code, globals(), None)
          runcommand('sed -i s/NUMUSLEVELS/'+str(nzrho_in)+'/ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i s/NUMWSLEVELS/'+str(nzrho_in+1)+'/ '+dirout+'LTRANS_'+identifier+'.data')
          listsigles=['READZETA','READDENS','READU','READV','READW','READSALT','READTEMP','READAKS','READuWIND','READvWIND','READiWIND']
          for count in range(0,maxnumfiles):
           if(FoundFile[count]):
             runcommand('sed -i s/'+listsigles[count]+'/.TRUE./ '+dirout+'LTRANS_'+identifier+'.data')
           else:
             runcommand('sed -i s/'+listsigles[count]+'/.FALSE./ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i s/IDENTIFIER/'+str(identifier)+'/ '+dirout+'LTRANS_'+identifier+'.data')
        
          runcommand('sed -i s/MITgcmFIELDSDIRECTORY/'+str(directory.replace('/','!?'))+'/ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/!?/\//g" '+dirout+'LTRANS_'+identifier+'.data')
          print('*')
          runcommand('sed -i "s/PREFIXZETA/'+FilesPrefix[0].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXDENS/'+FilesPrefix[1].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXUVEL/'+FilesPrefix[2].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXVVEL/'+FilesPrefix[3].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXWVEL/'+FilesPrefix[4].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXSALT/'+FilesPrefix[5].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXTEMP/'+FilesPrefix[6].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXAKS/'+FilesPrefix[7].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXUWIND/'+FilesPrefix[8].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXVWIND/'+FilesPrefix[9].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i "s/PREFIXIWIND/'+FilesPrefix[10].replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          #runcommand('sed -i "s/PREFIXLIST/'+(str(f_filenames)[13:-1]).replace("'",'?!')+'/" '+dirout+'LTRANS_'+identifier+'.data')
          print('***')
          runcommand('sed -i s/EXTENSION/?!'+f_filesextension+'?!/ '+dirout+'LTRANS_'+identifier+'.data')
          print('*****')
          runcommand('sed -i "s/?!/'+"'/g"+'"'+" "+dirout+"LTRANS_"+identifier+'.data')
          runcommand('sed -i "s/OUTPUTFREQ/'+str(int(abs(float(FreqOutput))))+'/g" '+dirout+"LTRANS_"+identifier+'.data')
          runcommand('sed -i s/NUMDAYS/'+str(NUMDAYS)+'/ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i s/TDIM/'+str(TDIM)+'/ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i s/NUMDIGITS/'+str(NUMDIGITS)+'/ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i s/FILENUM/'+str(FILENUM)+'/ '+dirout+'LTRANS_'+identifier+'.data')
          runcommand('sed -i s/FILESTEP/'+str(FILESTEP)+'/ '+dirout+'LTRANS_'+identifier+'.data')
def reversefindOccurences(s, ch):
        for i in range(len(s)-1,-1,-1):
                if (s[i]==ch):break
        return i

def findOccurencesnumber(s, ch):
        count=0
        for i in range(len(s)-1,-1,-1):
          for char in ch:
                if (s[i]==char):count+=1
        return count

def zerosfillednumberstring(number,stringsize):  
        outstring=str(number)
        len_num=len(outstring)
        if(len_num<stringsize):
                for count in range(len_num,stringsize): outstring='0'+outstring
        else:
                outstring=outstring[0:stringsize]
        return outstring

