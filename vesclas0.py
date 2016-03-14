'''
    VesClas v. 0.01 - detects and measures round donut-shaped microvesicle
    features from TSEM images of extracellular miocrovesicles. 
    Copyright (C) 2015  Petr Cizmar

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import cv2, cv, sys, numpy as np
import math as m
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from scipy.optimize import fmin
import string

## TWEAK THESE PARAMETERS TO SUIT YOUR IMAGES ##
DEBUG = 0

R_MAX=200 #default 100
BLUR=2.5 #2.5
DCT_ZERO_RAD=500 #default 500
MIN_GRADIENT=8 #default 8
MIN_BELOW_MAX=0.9 #default 0.9 (= minumim must be lower than 90% of maximum)
SCANSTEP=7 #default 7
RAY_STEP=0.5 #default 0.5
MAX_CIRCLE_ERROR=0.8
CROP_BOTTOM = 74 # default 74

MASK_MIN_DIFF=40 # default 40
MASK_OVERDRIVE=10

#SCANNING='no'
SCANNING='yes'

debname='image' 

#debug functions
#{{{
global debimg
def dbank(direction):
  global debimg
  global debimg_bank
  if direction == 'push':
    debimg_bank = debimg.copy()
  if direction == 'pull':
    debimg = debimg_bank

dcounter=0
def dimage():
  global dcounter
#  cv2.imshow('image', debimg) #Show the image
  fname = "vesistat_seq_%04d.png" % dcounter
  cv2.imwrite(fname, debimg)
  dcounter = dcounter + 1
##############################################################  #}}}

def interploated_value(img, x, y):#{{{
  a = img[m.floor(y),m.floor(x)]
  b = img[m.floor(y),m.ceil(x)]
  c = img[m.ceil(y),m.ceil(x)]
  d = img[m.ceil(y),m.floor(x)]
  dx = x - m.floor(x)
  dy = y - m.floor(y)
  e = a + dx * (b-a)
  f = d + dx * (c-d)
  g = e + dy * (f-e)
  # debug #print "int:",x,y,":",a,b,c,d,"->",g
  return g
###############################################################}}}
def cast_ray(imgf, x, y, phi):#{{{
  img_maxy = len(imgf)
  img_maxx = len(imgf[0])
  step = RAY_STEP
  dx = step*m.cos(phi*m.pi/180)
  dy = step*m.sin(phi*m.pi/180)
  ray = []
  xx = []
  yy = []
  tt = []
  v = interploated_value(imgf, x, y)
  max_found = 0
  min_found = 0
  gone_up = 0
  for i in range(1, R_MAX):
    v_old = v;
    xi = x + dx*i
    yi = y + dy*i
    if xi<2 or yi<2 or xi>(img_maxx-2) or yi>(img_maxy-2): 
      break
    xx = xx + [xi]
    yy = yy + [yi]
    tt = tt + [i]
    v = interploated_value(imgf, xi, yi)
    ray = ray + [v]
    up = 0
    if v > v_old:
      up = 1
      gone_up = 1
    if not max_found and not up and gone_up:
      max_found = 1
      v_max = i-1
    if max_found and not min_found and up and v < MIN_BELOW_MAX * ray[v_max]:
      min_found = 1
      v_min = i-1
# display plot window
#  plt.figure(1)
#  plt.plot(tt, ray)
#  if max_found: plt.plot(tt[v_max-1],ray[v_max-1],'rx')
#  if min_found: plt.plot(tt[v_min-1],ray[v_min-1],'gx')
#  plt.show(block=False)

  vmax = -1
  vmin = -1
  xmax = -1
  xmin = -1
  ymax = -1
  ymin = -1
  if min_found:
    trsmin = 0.2*(ray[v_max]-ray[v_min])+ray[v_min]
    trsmax = 0.8*(ray[v_max]-ray[v_min])+ray[v_min]
    for i in reversed(range(0,len(ray))):
      v = ray[i]
      if tt[i] < tt[v_min] and v >= trsmin and vmin == -1:  
        vmin = v
        xmin = xx[i]
        ymin = yy[i]
      if tt[i] < tt[v_min] and v >= trsmax and vmax == -1:
        vmax = v
        xmax = xx[i]
        ymax = yy[i]
        break
    stp = (vmax-vmin) / m.sqrt(m.pow(xmax-xmin,2)+m.pow(ymax-ymin,2))  
    #print 'stp:', stp
    if stp < MIN_GRADIENT:
      return -1,-1,-1
  return xmin,ymin,vmin
###############################################################}}}
def circ_err(x, xx, yy):#{{{
  cx = x[0]
  cy = x[1]
  r = x[2]
  l = len(xx)
  sq_dev = 0
  for i in range(0,l):
    sq_dev = sq_dev + (m.sqrt((cx - xx[i])**2 + (cy - yy[i])**2) - r)**2
  return sq_dev
###############################################################}}}
def fit_circle(xx, yy):#{{{
  cx = 0
  cy = 0
  l = len(xx)
# probe center is the mean (center of mass) of the found border points
  print xx
  cx = sum(xx) / len (xx)
  cy = sum(yy) / len (yy)
# probe radius is the mean distance from the probe center
  d = 0
  for i in range(0, len(xx)):
    d = d + (cx - xx[i])**2 + (cy - yy[i])**2
  d = m.sqrt(d / l)
# call fitting function
  print 'probe:', cx, cy, l
  xopt = fmin(circ_err,(cx,cy,d), (xx,yy), xtol=0.01)
  err = circ_err(xopt, xx, yy) / len(xx) # average quadratic deviation due to different border-point counts
  print 'result:', xopt, 'err: ', err
  return xopt[0], xopt[1], xopt[2], err # cx, cy, r
###############################################################}}}
def find_border(imgf, x, y, anglestep):#{{{
  bxx = []
  byy = []
  for phi in range(00,360,anglestep):
    bx, by, bv = cast_ray(imgf, x, y, phi)
    if bv > -1:
      bxx = bxx + [bx]
      byy = byy + [by]
      if DEBUG > 2:
        cv2.circle(debimg, (int(round(bx)),int(round(by))), 2, (0,255,0))
        dimage()
  return bxx, byy
###############################################################}}}
def explore_particle(imgf, x, y):#{{{
  message={}
  bxx, byy = find_border(imgf, x, y, 60) # first try
  if len(bxx) < 5:
    message['error']='Too few border points found.'
    return message
  x, y, r, err = fit_circle(bxx, byy)
  bxx, byy = find_border(imgf, x, y, 10) # corrective iteration
  if len(bxx) < 27:
    message['error']='Too few border points found.'
    return message
  x, y, r, err = fit_circle(bxx, byy)
  if err > MAX_CIRCLE_ERROR:
    message['error']='Too low roundness'
    return message
  if DEBUG > 2:
    dbank('pull')
    cv2.circle(debimg, (int(round(x)),int(round(y))), int(round(r)), (0,0,255), 2)
    dimage()
    dbank('push')
  message['center'] = (x,y)
  message['radius'] = r
  message['roundness_coeff'] = err
  message['border_points'] = len(bxx)
  return message

###############################################################}}}
def on_mouse(event, x, y, flags, par):#{{{
  if event == cv.CV_EVENT_LBUTTONUP:
    imgf = par['imgf']
    if DEBUG > 2:  
      cv2.circle(debimg, (x,y), 3, (0,120,255))
      dimage()
    explore_particle(imgf, x, y)
###############################################################}}}
def scan_image_brute_force(img):#{{{
  sizex = len(img[0])
  sizey = len(img)
  results=[]
  for y in range(SCANSTEP,sizey,SCANSTEP):
    for x in range(SCANSTEP,sizex,SCANSTEP):
      print x,y
      if DEBUG > 2:
        cv2.circle(debimg, (x,y), 3, (0,120,255))
        dimage()
      result = explore_particle(imgf, x, y)
      if not result.has_key('error'):
        results = results + [result]
  return results
###############################################################}}}
def scan_image_with_mask(img, mask):#{{{
  sizex = len(img[0])
  sizey = len(img)
  results=[]
  for y in range(SCANSTEP,sizey,SCANSTEP):
    for x in range(SCANSTEP,sizex,SCANSTEP):
      if mask[y,x] == 1:
        hit_vesicle = 0
        for r in results: # search in previous results
          c = r['center']
          r = r['radius']
          if (x-c[0])**2 + (y-c[1])**2 <= r**2: 
            print "HIT DUPL:", x, y
            hit_vesicle = 1
#        print x,y
        if hit_vesicle: continue    
        if DEBUG > 2:
          cv2.circle(debimg, (x,y), 3, (0,120,255))
          dimage()
        result = explore_particle(imgf, x, y)
        if not result.has_key('error'):
          results = results + [result]
  return results
###############################################################}}}
def calculate_mask(img):#{{{
# Threshold
  trsh = np.empty(img.shape, img.dtype)
  cv2.threshold(img, 1.3*img.mean(), 1, cv2.THRESH_BINARY, trsh)

  min_diff = MASK_MIN_DIFF

  mask = np.zeros(img.shape, 'uint8')
  sblx = cv2.Sobel(img, -1, 1, 0)
  print sblx.max()
  print sblx.min()
  thsblxl = np.empty(sblx.shape, sblx.dtype)
  thsblxr = np.empty(sblx.shape, sblx.dtype)
  cv2.threshold(sblx, min_diff, 1, cv2.THRESH_BINARY, thsblxl)
  cv2.threshold(-sblx, min_diff, 1, cv2.THRESH_BINARY, thsblxr)
  xcandidatess = np.uint8(thsblxl) & np.uint8(trsh)
  xcandidatesk = np.uint8(thsblxr) & np.uint8(trsh)

  sbly = cv2.Sobel(img, -1, 0, 1)
  thsblyl = np.empty(sbly.shape, sbly.dtype)
  thsblyr = np.empty(sbly.shape, sbly.dtype)
  cv2.threshold(sbly, min_diff, 1, cv2.THRESH_BINARY, thsblyl)
  cv2.threshold(-sbly, min_diff, 1, cv2.THRESH_BINARY, thsblyr)
  ycandidatess = np.uint8(thsblyl) & np.uint8(trsh)
  ycandidatesk = np.uint8(thsblyr) & np.uint8(trsh)

  candidates = (xcandidatess | ycandidatess | xcandidatesk | ycandidatesk)
  kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(2,2))
  candidates = cv2.erode(candidates,kernel)
#search for contours
  contours, hierarchy = cv2.findContours(candidates, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
  candidates = cv2.cvtColor(candidates*255, cv2.COLOR_GRAY2BGR)
  print "Contour:", contours[0]
  if DEBUG > 2: dbank('push')
  for co in contours:
    rect = cv2.boundingRect(co)
    print "rect:", rect
    xl = rect[0]-MASK_OVERDRIVE
    if xl < 0: xl = 0
    xh = rect[0]+rect[2]+MASK_OVERDRIVE
    if xh >= img.shape[1]: xh = img.shape[1]-1
    yl = rect[1]-MASK_OVERDRIVE
    if yl < 0: yl = 0
    yh = rect[1]+rect[3]+MASK_OVERDRIVE
    if yh >= img.shape[0]: yh = img.shape[0]-1
    mask[yl:yh,xl:xh] = 1
    if DEBUG > 2: # currently broken
      cv2.rectangle(debimg,(rect[0]-MASK_OVERDRIVE,rect[1]-MASK_OVERDRIVE),(rect[0]+rect[2]+MASK_OVERDRIVE,rect[1]+rect[3]+MASK_OVERDRIVE), (255,0,0), -2)
      cv2.drawContours(debimg, contours, -1, (0,0,255), 1)
      dimage()
  if DEBUG > 2: # currently broken
    dimage()
    dbank('pull')
  return mask
###############################################################}}}
def prepare_image(img):#{{{
#crop off the bottom
  img = img[:-CROP_BOTTOM,:]
#convert to floats
  img = np.float32(img)
#Discrete Cosine Transform
  cimg = img
  cv2.dct(img, cimg, 0);
# Zero the top-left corner - round
  dct_zero_radius = DCT_ZERO_RAD 
  for j in range(0, cimg.shape[1]):
     for i in range(0, cimg.shape[0]):
        if i**2 + j**2 <= dct_zero_radius:
           cimg[i,j] = 0
# Inverse DCT
  evbg = cimg
  cv2.idct(cimg, evbg)
# Gaussian Blur
  blr = cv2.GaussianBlur( evbg, (0,0), BLUR) 
  blr = (blr - blr.min())
  blr = blr / blr.max() * 255
  #blr = cv2.adaptiveThreshold(blr, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY,5, 0)
  return blr
###############################################################}}}
def get_pixel_size_from_tiff(filename):# {{{
  pixel_size = 1
  unit = "px"
  file = open(filename, "r")
  for line in file:
    if string.find(line, "Pixel Size = ") > -1:
      _,pixel_size,unit = string.rsplit(line, " ", 2)
  file.close()
  return string.atof(pixel_size), string.rstrip(unit,"\n\r")# }}}
########### MAIN
if __name__ == "__main__":#{{{
  global debimg
  global debimg_bank

  print('''
  VesClas v. 0.01  Copyright (C) 2015  Petr Cizmar

    This program comes with ABSOLUTELY NO WARRANTY.
    This is free software, and you are welcome to redistribute it
    under certain conditions; see source code for details.

  ''')


  if len(sys.argv) < 2:
    print( "Usage: %s [image-file name]" % sys.argv[0] )
    sys.exit(127)
  
  pixel_size, unit = get_pixel_size_from_tiff(sys.argv[1]) 
  file_base,file_ext = string.rsplit(sys.argv[1],".",1) # get the file-name base
# read the image
  src = cv2.imread(sys.argv[1], cv2.CV_LOAD_IMAGE_GRAYSCALE)
# prepare image
  imgf = prepare_image(src)
# display it


  img = np.uint8(imgf)
#  imgd = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
  srcd = cv2.cvtColor(src, cv2.COLOR_GRAY2BGR)
  cv2.imwrite("bgblr.tiff", img)
# set the mouse callback function
  debimg = srcd
  dimage()
  dbank('push')
#calculate mask

  
  if SCANNING=='yes':
    mask = calculate_mask(imgf)

    if DEBUG > 2:
      maske = np.zeros(debimg.shape, debimg.dtype)
      maskd = cv2.cvtColor(mask, cv2.COLOR_GRAY2BGR)
      maske[0:maskd.shape[0],:] = maskd
      print "MASK x IMG:", debimg.shape, debimg.dtype, maskd.shape, maskd.dtype
      debimg = cv2.addWeighted(debimg, 1, 1-maske, -50, 0)
      dimage()
      dbank('push')

  par={'imgf':imgf}
  if not SCANNING=='yes':
    cv.SetMouseCallback('image', on_mouse, par)
    cv2.waitKey()
  else:
    #results = scan_image_brute_force(imgf)
    results = scan_image_with_mask(imgf, mask)
    srcd = cv2.cvtColor(src, cv2.COLOR_GRAY2BGR)
    result_file = open(file_base + '_result' + ".txt", 'w')
    result_file.write("# index, center_x [px], center_y [px], radius [%s], sigma_divergence [%s], valid_border_points\n" % (unit, unit))
    for i,res in enumerate(results):
      print res
      cx = int(round(res['center'][0]))
      cy = int(round(res['center'][1]))
      cv2.circle(srcd, (cx, cy), int(round(res['radius'])), (0,0,255), 2)
      txt = "%d" % i
      tsize = cv2.getTextSize(txt, cv.CV_FONT_HERSHEY_COMPLEX_SMALL, 0.8, 1)
      cv2.putText(srcd, txt, (cx-tsize[0][0]/2,cy+tsize[0][1]/2), cv.CV_FONT_HERSHEY_COMPLEX_SMALL, 0.8 , (255,255,255))
      result_file.write("%d, %0.3f, %0.3f, %0.3f, %0.3f, %d\n" % (i,
      res['center'][0], res['center'][1], res['radius'] * pixel_size,
      m.sqrt(res['roundness_coeff']) * pixel_size, res['border_points']))

    result_file.close()
    if DEBUG > 0:
      cv2.imshow('Result', srcd)
    #result_image = file_base + '_result.' + file_ext # if everything works OK.
    result_image = file_base + '_result.png' # in case tiff lib does not save properly.
    cv2.imwrite(result_image, srcd)

    cv2.waitKey()

#}}}
# vim: smartindent:expandtab:shiftwidth=2:softtabstop=2

