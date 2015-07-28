 ** If gimp can save xpm:
If size is too big: change number of colors to 32 may help size
(gimp: image>mode>indexed, use normal dithering)
(need to save the images in xpm file first).

If gimp can't save xpm:
First save in ppm
Use linux ppmtoxpm a.ppm > a.xpm
Do it for bunches of images to save time

foreach f ( `ls` )
	set o = `echo $f | awk -F. '{print $1}'`
	ppmtoxpm $f > $o"_xpm.h"
end


For transparency:
1 Select region to be transparent
2 Colors > To Alpha
3 Notice: jpg doesn't support transparency but xpm does.

# For icons:
# After saving in xpm, set size to 256x256

# mv *.xpm ../src/images/
# Then change the type to static const char *
# May want to change name too (remove cpm)
# Then move it to _xpm.h  (then it's same name as the object)





Alternatively:


1. Save original png in the "original" directory (256x256)
2. Convert to xpm: http://www.office-converter.com/Convert-to-XPM
3. Change the type to static const char *
4. Then move it to _xpm.h  (then it's same name as the object)
