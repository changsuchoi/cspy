import os
import glob

def calibfitslist(caliblist):
    oklist,nolist=[],[]
    currentdir=os.getcwd()
    currentdir_split=currentdir.split('/')
    for i, j in enumerate(caliblist):
        # sample='Calib-DOAO-NGC3367-20180330-111526-B-60.fits'    
        # print(i,j)
        filename = j.split('.')
        components=j.split('.')[0].split('-')
        if len(components) !=7 : 
            print (i,j)
            nolist.append(j)
        if  components[0] != 'Calib' : 
            os.system('rename Calibrated Calib '+j)
            print (i,j,'rename filename Calib-')
            nolist.append(j)
        if components[1] not in currentdir_split:
            print (i,j,'Different OBS')
            nolist.append(j)
        if components[2] not in currentdir_split:
            print (i,j,'Different Object')
            nolist.append(j)
        if len(components[3]) != 8 or  components[3][:2] !='20' :
            print(i,j,'Date value is strange')
            nolist.append(j)        
        if len(components[4]) != 6 : 
            print(i,j,'Time value is strange')  
            nolist.append(j)
        else: oklist.append(j)  
    return oklist, nolist    


if __name__=="__main__": 


    caliblist=glob.glob("Calib*.fits")
    oklist, nolist=calibfitslist(caliblist)

def calibfitslist(calibstr):
    caliblist=glob.glob("Calib*.fits")
    oklist, nolist=calibfitslist(caliblist)
    print('oklist', len(oklist),'nolist', len(nolist))
    return oklist, nolist
   