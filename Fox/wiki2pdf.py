import urllib
import urlparse
import re
import os
from time import sleep

class MoinSpider:
   def __init__(self,site='objcryst.sourceforge.net',
                exclude=["RecentChanges","action=",
                         "FindPage","TitleIndex","WordIndex",
                         "Help","template","Template","MoinMoin",
                         "UserPreferences","WikiSandBox",
                         "ScriptAlias","ScriptAlias"]):
      self.u=urllib.URLopener()
      self.u.addheader(('USER_AGENT', 'Mozilla/4.0'))
      self.base='href=\"/Fox/'
      self.suffix="?action=print"
      self.site=site
      self.pages=[] # list of pairs [relative URL, page content]
      self.d={} # dictionnary with keys=relative URL, value= short filename for the downloaded page
      self.exclude=exclude
      self.nbFail=0 # pages which failed to load
      self.img=set()
   def Weave(self, lnk='/Fox/FoxWiki',nbtry=3):
      """ Download recursively all pages, starting from one relative URL.
      """
      if self.d.has_key(lnk): # we already got that page !
         return
      self.d[lnk]="wiki_%i.html"%(1000+len(self.d))
      url="http://"+self.site+lnk+self.suffix #:TODO: use urlparse !
      print
      print "Getting page: %s"%url
      print "         -> %s"%(self.d[lnk])
      nb=nbtry
      cont=True
      while(nb>0):
         try:
            p=self.u.open(url)
            page=p.read()
            nb=-1
         except IOError:
            nb-=1
            print "IOError..... retry #%i"%(nbtry-nb)
            sleep(1)
      if nb==0:
         print "Failed to load page, after %i trials:"%nbtry,lnk
         self.nbFail+=1
         return
      if re.search("This page does not exist yet",page)!=None:
         print "  -> Page has not been written yet !"
         self.d[lnk]="http://"+self.site+lnk # Link directly to site
         return
      self.pages.append([lnk,page])
      for m in re.finditer(r"href\=\"(.*?)\"",page):
         newlink=m.group()
         if len(newlink)>=len(self.base):
            if newlink[:len(self.base)]==self.base:
               keep=True
               for x in self.exclude:
                  if re.search(x,newlink)!= None:
                     keep=False
                     break
               if keep:
                  #print "    ->%s"%newlink
                  newlink=newlink[6:-1]# [6:-1] -> exlude ' href=" ' and the end ' " '
                  newlink=re.split('#',newlink)[0] # exclude anchors
                  self.Weave(newlink) 
               #else:
               #   print "    ->%s ?  NO"%newlink
            
   def WeaveStatic(self, pagelist,nbtry=3):
      """ Alternative to weave: download a pre-selected list of pages
      """
      for lnk in pagelist:
         self.d[lnk]="wiki_%i.html"%(1000+len(self.d))
         url="http://"+self.site+lnk+self.suffix #:TODO: use urlparse !
         print "Getting page: %s         -> %s"%(url,self.d[lnk])
         nb=nbtry
         cont=True
         while(nb>0):
            try:
               p=self.u.open(url)
               page=p.read()
               nb=-1
            except IOError:
               nb-=1
               print "IOError..... retry #%i"%(nbtry-nb)
               sleep(1)
            if nb==0:
               print "Failed to load page, after %i trials:"%nbtry,lnk
            if re.search("This page does not exist yet",page)!=None:
               print "  -> Page has not been written yet !"
               self.d[lnk]="http://"+self.site+lnk # Link directly to site
               nb=0
            else:
               self.pages.append([lnk,page])
         
   def Pages2Html(self,d="wikihtml"):
      #TODO : remove links to non-written pages
      if not os.path.exists(d):
         os.mkdir(d)
      #this is necessary so that urls that contain other (smaller) urls
      #are replaced first
      ks=self.d.keys()
      ks.sort(reverse=True)
      
      for p in self.pages:
         for m in re.finditer(r"img .*? src\=\"(.*?)\"",p[1]):
            print re.findall(r"src\=\"(.*?)\"",m.group())
            url=re.findall(r"src\=\"(.*?)\"",m.group())[0]
            up=urlparse.urlparse(url)
            print url
            up0,up1,up2,up3,up4,up5=up[0],up[1],up[2],up[3],up[4],up[5]
            if up4 != '':
               name=re.split('=',up4).pop()
            else:
               name=re.split('/',up2).pop()
            if name not in self.img:#download image once
               self.img.add(name)
               if up0=='':
                  up0='http'
               if up1=='':
                  up1=self.site
               urlimg=urlparse.urlunparse((up0,up1,up2,up3,up4,up5))
               print "   %s -> %s"%(urlimg,name)
               nbTry=3
               nb=nbTry
               while nb>0:
                  try:
                     urllib.urlretrieve(urlimg,d+"/"+name)
                     nb=-1
                  except IOError:
                     nb-=1
                     print "IOError..... retry #%i to get %s"%(nbTry-nb,name)
                     sleep(1)
               if nb==0:
                      print "Failed to load image, after %i trials: %s"%(nbtry,name)
               else: # KLUDGE png->png cause htmldoc chokes on these...
                  if name[-4:]==".png":
                    print "convert %s %s"%(d+"/"+name,d+"/"+name[:-3]+"jpg")
                    os.system("convert %s %s"%(d+"/"+name,d+"/"+name[:-3]+"jpg"))
                    os.system("rm -f %s"%(d+"/"+name))
            p[1]=p[1].replace(url,name)
         for k in ks:# change to local url
            if k!=self.d[k]:
               p[1]=p[1].replace(k,self.d[k])
            # Change src field of img from "wiki_1002.html?action=AttachFile&amp;do=get&amp;target=toto.jpg" to "toto.jpg"
            p[1]=p[1].replace("%s?action=AttachFile&amp;do=get&amp;target="%k,"")
         p[1]=p[1].replace(".png",".jpg")
         f=open(d+"/"+self.d[p[0]],'w')
         f.write(p[1])
   def Html2pdf(self,d="wikihtml"):
      os.system("mogrify -resize '600x>' wikihtml/*.jpg")
      #os.system("htmldoc --jpeg=85 --webpage %s/*.html --linkcolor blue -f wiki.pdf"%d)
      os.system("htmldoc --jpeg=85 --webpage %s/*.html --linkcolor blue --size a4 --format pdf14 --links --book --toclevels 3 --left 1.5cm --right 1.5cm --top 1.5cm --bottom 1.5cm --footer Dc1 -f wiki.pdf"%d)   
      #os.system("rm -f wikihtml/*")

#m=MoinSpider(site="objcryst.sourceforge.net")
m=MoinSpider(site="vincefn.net")

m.WeaveStatic(["/Fox/FoxWiki",
               "/Fox/BiblioReferences",
               "/Fox/Screenshots",
               "/Fox/Download",
               "/Fox/Install",
               "/Fox/Install/Linux",
               "/Fox/Install/MacOSX",
               "/Fox/Install/Windows",
               "/Fox/Tutorials",
               "/Fox/Tutorials/Cimetidine",
               "/Fox/Tutorials/PbSO4",
               "/Fox/Manual",
               "/Fox/Manual/Crystal",
               "/Fox/Manual/Powder",
               "/Fox/Manual/Indexing",
               "/Fox/Manual/ProfileFitting",
               "/Fox/Manual/SingleCrystalDiffraction",
               "/Fox/Manual/Algorithm",
               "/Fox/Manual/CIF",
               "/Fox/Manual/Tips",
               "/Fox/MailingList",
               "/Fox/Changelog",
               "/Fox/Compile/Linux",
               "/Fox/Compile/MacOSX",
               "/Fox/Compile/Windows",
               #"/Fox/BiblioStructures",
               #"/Fox/VincentFavreNicolin"
               ])

m.Pages2Html()
m.Html2pdf()

