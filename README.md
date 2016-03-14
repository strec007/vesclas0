# Vesclas v. 0.01

&copy; 2015 Petr Cizmar @ PTB Germany

This is a **very preliminary** version of the *VesClas* software. It **only**
can detect some **donut-shaped** features. As of March 14th 2016, development of
a new version is undergo and will soon be released. The new version can detect a
lot more features and is also significantly faster than this version. 

## Files

* *README.md* is this file,
* *presentation* directory contains a short presentation about this technique
  includind two videos.
* *LICENSE* contains the license, and
* *example_tsem_image.tiff* is an example of a TSEM image of extracellular
  microvesicles
* vesclas0.py is the source code.


## Requirements
This program requires:

* python v .2

and following python libraries:

* OpenCV2
* matplotlib
* numpy
* scipy

## Invokation

Invoke the program as follows:

    python vesclas0.py example_tsem_image.tiff

## Teaking of parameters

At the beginning of the source code, there is a parameter section. If the
program cannot detect anything in your images, you might need to tweak some
parameters. The original setting looks like this:

```
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
```

where

* *R_MAX* is the maximum feature radius. Whatever is bigger, will be discarded.
* *BLUR* is the sigma of the Gaussian blurring. If you have a very noisy image,
  you might need to add a bit here.
* *DCT_ZERO_RAD* is the radius of the area in direct-cosine-transform image to
  be zeroed. It removes the low-frequency info from the image. Do not change,
  unless you need to.
* *MIN_GRADIENT* is a minimum derivative of the border transition.
* *MIN_BELOW_MAX* is a lower limit of the height of the border transition.
* *SCANSTEP* is the length of the scanning step throughout the entire image.
* *RAY_STEP* is the length of the step for the ray test.
* *MAX_CIRCLE_ERROR* is maximum allowed circle-fitting error. Decrease if you
  want to keep only very round features.
* *CROP_BOTTOM* is the height of the info-bar to be cropped off.
* *MASK_MIN_DIFF* is a minimum derivative to unmask.
* *MASK_OVERDRIVE* is a size of the square area to extend each unmasked pixel.


## Credits

This program has been written by Petr Cizmar at Physikalisch-Technische
Bundesanstalt (PTB), Braunschweig, Germany as a part of the work on the MetVes
EMRP project (<http://www.metves.eu>). 

## License

Copyright (c) 2015, Petr Cizmar.    
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

