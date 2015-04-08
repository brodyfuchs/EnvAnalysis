from numpy import ma,array,linspace,log,cos,sin,pi,zeros,exp,arange
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.ticker import FixedLocator, AutoLocator, ScalarFormatter
import matplotlib.transforms as transforms
import matplotlib.axis as maxis
import matplotlib.artist as artist
from matplotlib.projections import register_projection
from matplotlib.pyplot import rcParams,figure,show,draw
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

from thermodynamics import *

from UserDict import UserDict
from datetime import datetime
import os,sys

class SkewXTick(maxis.XTick):
    #Copyright (c) 2008 Ryan May
    def draw(self, renderer):
        if not self.get_visible(): return
        renderer.open_group(self.__name__)

        if self.gridOn:
            self.gridline.draw(renderer)
        if self.tick1On:
            self.tick1line.draw(renderer)
        if self.tick2On:
            self.tick2line.draw(renderer)

        if self.label1On:
            self.label1.draw(renderer)
        if self.label2On:
            self.label2.draw(renderer)

        renderer.close_group(self.__name__)

    def set_clip_path(self, clippath, transform=None):
        artist.Artist.set_clip_path(self, clippath, transform)
        self.tick1line.set_clip_path(clippath, transform)
        self.tick2line.set_clip_path(clippath, transform)
        self.gridline.set_clip_path(clippath, transform)
    set_clip_path.__doc__ = artist.Artist.set_clip_path.__doc__

class SkewXAxis(maxis.XAxis):
    #Copyright (c) 2008 Ryan May
    def _get_tick(self, major):
        return SkewXTick(self.axes, 0, '', major=major)

    def draw(self, renderer, *args, **kwargs):
        'Draw the axis lines, grid lines, tick lines and labels'
        ticklabelBoxes = []
        ticklabelBoxes2 = []

        if not self.get_visible(): return
        renderer.open_group(__name__)
        interval = self.get_view_interval()
        for tick, loc, label in self.iter_ticks():
            if tick is None: continue
            if transforms.interval_contains(interval, loc):
                tick.set_label1(label)
                tick.set_label2(label)
            tick.update_position(loc)
            tick.draw(renderer)
            if tick.label1On and tick.label1.get_visible():
                extent = tick.label1.get_window_extent(renderer)
                ticklabelBoxes.append(extent)
            if tick.label2On and tick.label2.get_visible():
                extent = tick.label2.get_window_extent(renderer)
                ticklabelBoxes2.append(extent)

        # scale up the axis label box to also find the neighbors, not
        # just the tick labels that actually overlap note we need a
        # *copy* of the axis label box because we don't wan't to scale
        # the actual bbox

        self._update_label_position(ticklabelBoxes, ticklabelBoxes2)

        self.label.draw(renderer)

        self._update_offset_text_position(ticklabelBoxes, ticklabelBoxes2)
        self.offsetText.set_text( self.major.formatter.get_offset() )
        self.offsetText.draw(renderer)

class SkewXAxes(Axes):
    #Copyright (c) 2008 Ryan May
    # The projection must specify a name.  This will be used be the
    # user to select the projection, i.e. ``subplot(111,
    # projection='skewx')``.
    name = 'skewx'

    def _init_axis(self):
        #Taken from Axes and modified to use our modified X-axis
        "move this out of __init__ because non-separable axes don't use it"
        self.xaxis = SkewXAxis(self)
        self.yaxis = maxis.YAxis(self)
        self._update_transScale()

    def draw(self, *args):
        '''
        draw() is overridden here to allow the data transform to be updated
        before calling the Axes.draw() method.  This allows resizes to be
        properly handled without registering callbacks.  The amount of
        work done here is kept to a minimum.
        '''
        self._update_data_transform()
        Axes.draw(self, *args)

    def _update_data_transform(self):
        '''
        This separates out the creating of the data transform so that
        it alone is updated at draw time.
        '''
        # This transforms x in pixel space to be x + the offset in y from
        # the lower left corner - producing an x-axis sloped 45 degrees
        # down, or x-axis grid lines sloped 45 degrees to the right
        self.transProjection.set(transforms.Affine2D(
            array([[1, 1, -self.bbox.ymin], [0, 1, 0], [0, 0, 1]])))

        # Full data transform
        self.transData.set(self._transDataNonskew + self.transProjection)

    def _set_lim_and_transforms(self):
        """
        This is called once when the plot is created to set up all the
        transforms for the data, text and grids.
        """
        #Get the standard transform setup from the Axes base class
        Axes._set_lim_and_transforms(self)

        #Save the unskewed data transform for our own use when regenerating
        #the data transform. The user might want this as well
        self._transDataNonskew = self.transData

        #Create a wrapper for the data transform, so that any object that
        #grabs this transform will see an updated version when we change it
        self.transData = transforms.TransformWrapper(
            transforms.IdentityTransform())

        #Create a wrapper for the proj. transform, so that any object that
        #grabs this transform will see an updated version when we change it
        self.transProjection = transforms.TransformWrapper(
            transforms.IdentityTransform())
        self._update_data_transform()

    def get_xaxis_transform(self, which='grid'):
        """
        Get the transformation used for drawing x-axis labels, ticks
        and gridlines.  The x-direction is in data coordinates and the
        y-direction is in axis coordinates.

        We override here so that the x-axis gridlines get properly
        transformed for the skewed plot.
        """
        return self._xaxis_transform + self.transProjection

    # Disable panning until we find a way to handle the problem with
    # the projection
    def start_pan(self, x, y, button):
        pass

    def end_pan(self):
        pass

    def drag_pan(self, button, key, x, y):
        pass

    def other_housekeeping(self,pmin=100.,mixratio=array([])):
	# Added by Thomas Chubb
	self.yaxis.grid(True,ls='-',color='y',lw=0.5)

	# Plot x-grid lines instead of using xaxis.grid().
	# This is because xaxis.grid only plots skew lines 
	# that intersect the lower x axis... a possible fix
	# would be to use twinx() to plot upper and lower
	# grid lines but I think that's messy too.
	for TT in linspace(-100,100,21):
	    self.plot([TT,TT],[1050,pmin],color='y',lw=0.5)

	# self.set_ylabel('Pressure (hPa)')
	self.set_xlabel('Temperature ($^{\circ}$C)',fontsize=12)
        self.set_ylabel('Pressure (hPa)',fontsize=12)
        self.yaxis.labelpad = -2
	self.set_yticks(linspace(100,1000,10))
	self.yaxis.set_major_formatter(ScalarFormatter())
	self.set_xlim(-40,45)
	self.set_ylim(1050.,pmin)
	self.spines['right'].set_visible(False)
	self.get_yaxis().set_tick_params(which="both",size=0)
	self.get_xaxis().set_tick_params(which="both",size=0)

    def set_xticklocs(self,xticklocs):
	# Added by Thomas Chubb
	self.set_xticks(xticklocs)

    def add_dry_adiabats(self,T0,P,do_labels=True,**kwargs):
	# Added by Thomas Chubb
	P0=1000.
	T=array([ (st+273.15)*(P/P0)**(Rs_da/Cp_da)-273.15 for st in T0 ])
	labelt=[ (st+273.15)*1**(Rs_da/Cp_da) for st in T0 ]
	if kwargs.has_key('color'): 
	    col=kwargs['color']
	else: 
	    col='k'
	for tt,ll in zip(T,labelt):
	    self.plot(tt,P,**kwargs)
	    if do_labels:
		if (tt[8]>-50) and (tt[8]<40):# 20 is default, how far to label the dry adibats
		    self.text(tt[8],P[8]+30,'%d'%(ll),fontsize=8,\
			    ha='center',va='bottom',rotation=-40,color=col)#, bbox={'facecolor':'w','edgecolor':'w'})
# BF LAST PART IS THE BOX BEHIND THE TEXT I THINK

	return T
    

    def add_moist_adiabats(self,T0,P,do_labels=True,**kwargs):
	# Added by Thomas Chubb
	T=array([lift_wet(st,P) for st in T0])
	if kwargs.has_key('color'): 
	    col=kwargs['color']
	else: 
	    col='k'
	for tt in T:
	    self.plot(tt,P,**kwargs)
	    # if (tt[-1]>-60) and (tt[-1]<-10):
	    if do_labels:
		self.text(tt[-1],P[-1],'%d'%tt[0],ha='center',va='bottom',\
			fontsize=8, bbox={'facecolor':'w','edgecolor':'w'},color=col)

    def add_mixratio_isopleths(self,w,P,do_labels=True,**kwargs):
	# Added by Thomas Chubb
	e=array([P*ww/(.622+ww) for ww in w])
	T = 243.5/(17.67/log(e/6.112) - 1)
	if kwargs.has_key('color'): 
	    col=kwargs['color']
	else: 
	    col='k'

	for tt,mr in zip(T,w):
	    self.plot(tt,P.flatten(),**kwargs)
	    if do_labels:
		if (tt[-1]>-45) and (tt[-1]<20):
		    if mr*1000<1.:
			fmt="%4.1f"
		    else:
			fmt="%d"
		    self.text(tt[-1],P[-1],fmt%(mr*1000),\
			    color=col, fontsize=8,ha='center',va='bottom',\
			    bbox={'facecolor':'w','edgecolor':'w'})

# Now register the projection with matplotlib so the user can select
# it.
register_projection(SkewXAxes)

class Sounding(UserDict):
    # Copyright (c) 2013 Thomas Chubb 
    """Utilities to read, write and plot sounding data quickly and without fuss
    
    INPUTS:
    filename:   If creating a sounding from a file, the full file name. The 
		format of this file is quite pedantic and needs to conform 
		to the format given by the University of Wyoming soundings 
		(see weather.uwyo.edu/upperair/sounding.html) 
    data: 	Soundings can be made from atmospheric data. This should be 
		in the form of a python dict with (at minimum) the following 
		fields:

		TEMP: dry-bulb temperature (Deg C)
		DWPT: dew point temperature (Deg C)
		PRES: pressure (hPa)
		SKNT: wind speed (knots)
		WDIR: wind direction (deg)

		The following fields are also used, but not required by the 
		plot_skewt routine:

		HGHT (m)
		RELH (%)
		MIXR (g/kg)
		THTA (K)
		THTE (K)
		THTV (K)
    """


    def __init__(self,filename=None,data=None):
	UserDict.__init__(self)

	if data is None:
	    self.data={}
	    self.readfile(filename)
	else:
	    self.data=data
	    self['SoundingDate']=""

    def plot_skewt(self, imagename = None, title = None, text = None, mixdepth = 100, pres_s = 1000.0,
        cape = None, ncape = None, **kwargs):
	"""A wrapper for plotting the skewt diagram for a Sounding instance."""
	
	self.make_skewt_axes()
	self.add_profile(pres_s = pres_s, **kwargs)
	parcel=self.surface_parcel(mixdepth = mixdepth, pres_s = pres_s)
#	print parcel
	self.lift_parcel(*parcel,cape=cape,ncape=ncape,**kwargs)
        # here, could calc interp_parcel_path() and calc_cape() to calculate CAPE;
        # the only problem is that lift_parcel() labels CAPE and NCAPE on the sounding,
        # and lift_parcel is already called. One possible solution is to have the string
        # that's created in lift_parcel be defined outside lift_parcel()

	if isinstance(title, str):
	    self.skewxaxis.set_title(title,fontsize=14)
	else:
	    self.skewxaxis.set_title("%s: %s"%(self["StationNumber"],self['SoundingDate']),fontsize=14)

        #if text is not None:
            #self.plt.figtext(0.1,0.7,text,fontsize=25)
        #if cape is not None:
        #    self.fig.text(0.1,0.71,"CAPE: %0.0f J kg$^{-1}$\n"%cape)
        #if ncape is not None:
        #    self.fig.text(0.1,0.69,"NCAPE: %0.2f m s$^{-2}$\n"%ncape)

	if imagename is not None:
	    print("saving figure")
	    self.fig.savefig(imagename,dpi=100,bbox_inches='tight')


    def add_profile(self,bloc=0.5, pres_s = 1000.0, **kwargs):
	"""Add a new profile to the SkewT plot.

	This is abstracted from plot_skewt to enable the plotting of
	multiple profiles on a single axis, by updating the data attribute.
	For example:

	S=SkewT.Sounding(data={})
	S.make_skewt_axes()
	S.readfile("../examples/94975.2013062800.txt")
	S.add_profile(color="b",bloc=0.5)
	S.readfile("../examples/94975.2013070900.txt")
	S.add_profile(color="r",bloc=1.)

	Use the kwarg 'bloc' to set the alignment of the wind barbs from the centerline
        (useful if plotting multiple profiles on the one axis)


	Modified 25/07/2013: enforce masking of input data for this 
	function (does not affect the data attribute).
	"""

	try: pres = ma.masked_invalid(self.data['pres'])
	except KeyError: raise KeyError, "Temperature in hPa (PRES) is required!"

	try: tc=ma.masked_invalid(self.data['temp'])
	except KeyError: raise KeyError, "Temperature in C (TEMP) is required!"

	try: dwpt=ma.masked_invalid(self.data['dwpt'])
	except KeyError:
	    print "Warning: No DWPT available"
	    dwpt=ma.masked_array(zeros(pres.shape),mask=False)

	try:
	    sknt=self.data['sknt']
	    drct=self.data['drct']
	    rdir = (270.-drct)*(pi/180.)
	    uu = ma.masked_invalid(sknt*cos(rdir))
	    vv = ma.masked_invalid(sknt*sin(rdir))
	except KeyError:
	    print "Warning: No SKNT/DRCT available"
	    uu=ma.masked_array(zeros(pres.shape),mask=True)
	    vv=ma.masked_array(zeros(pres.shape),mask=True)

	vp = pres <= pres_s # valid pressures to plot
#	pres_plot = pres[valid_pres]
#	tc_plot = tc[valid_pres]
#	dwpt_plot = dwpt[valid_pres]
# BF OLD WAY
#	tcprof=self.skewxaxis.plot(tc, pres, zorder=5,**kwargs)
#	dpprof=self.skewxaxis.plot(dwpt, pres, zorder=5,**kwargs)
# NEW WAY TO ONLY PLOT ABOVE SURFACE PARCEL
	tcprof=self.skewxaxis.plot(tc[vp], pres[vp], zorder=5,**kwargs)
	dpprof=self.skewxaxis.plot(dwpt[vp], pres[vp], zorder=5,**kwargs)


	# this line should no longer cause an exception
	nbarbs=(~uu.mask).sum()

	skip=max(1,int(nbarbs/32))

	if kwargs.has_key('color'): bcol=kwargs['color']
	else: bcol='k'

	if kwargs.has_key('alpha'): balph=kwargs['alpha']
	else: balph=1.

#	self.wbax.barbs((zeros(pres.shape)+bloc)[::skip]-0.5, pres[::skip],\
#		uu[::skip], vv[::skip],\
#		length=6,color=bcol,alpha=balph,lw=0.5)

	self.wbax.barbs((zeros(pres[vp].shape)+bloc)[::skip]-0.5, pres[vp][::skip],\
		uu[vp][::skip], vv[vp][::skip],\
		length=6,color=bcol,alpha=balph,lw=0.5)


	self.skewxaxis.other_housekeeping()

	return tcprof

    def make_skewt_axes(self,pmax=1050.,pmin=100.):
	
	self.fig = figure(figsize=(7,8))
	self.fig.clf()
	
	rcParams.update({\
		'font.size':10,\
		})

	self.skewxaxis=self.fig.add_axes([.085,.1,.83,.8], projection='skewx')
	self.skewxaxis.set_yscale('log')

	xticklocs=arange(-80,45,10)
	T0 = xticklocs

	P=linspace(pmax,pmin,37)

	w = array([0.0001,0.0004,0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032])
	self.skewxaxis.add_mixratio_isopleths(w,P[P>=550],color='purple',ls='--',alpha=1.,lw=0.5)
	self.skewxaxis.add_dry_adiabats(linspace(250,440,20)-273.15,P,color='red',ls='--',alpha=1.,lw=0.5)
	self.skewxaxis.add_moist_adiabats(linspace(8,32,7),P[P>=200],color='g',ls='--',alpha=1.,lw=0.5)
	self.skewxaxis.other_housekeeping()

	self.wbax=self.fig.add_axes([0.815,0.1,0.1,0.8],sharey=self.skewxaxis,frameon=False)
	self.wbax.xaxis.set_ticks([],[])
	self.wbax.yaxis.grid(True,ls='-',color='y',lw=0.5)
	for tick in self.wbax.yaxis.get_major_ticks():
	    # tick.label1On = False
	    pass
	self.wbax.get_yaxis().set_tick_params(size=0,color='y')
	self.wbax.set_xlim(-1.5,1.5)
	self.wbax.get_yaxis().set_visible(False)

	# Set up standard atmosphere height scale on 
	# LHS of plot. It's jus
	majorLocatorKM   = MultipleLocator(2)
	majorLocatorKFT  = MultipleLocator(5)
	minorLocator     = MultipleLocator(1)

	self.kmhax=self.fig.add_axes([0.915,0.1,1e-6,0.8],frameon=True)
	self.kmhax.xaxis.set_ticks([],[])
	self.kmhax.spines['left'].set_color('k')
	self.kmhax.spines['right'].set_visible(False)
	self.kmhax.tick_params(axis='y', colors='k',labelsize=10)
	self.kmhax.set_ylim(0,16.18)
	self.kmhax.set_title("km/kft",fontsize=12)
	self.kmhax.get_yaxis().set_tick_params(which="both",direction='out')
	self.kmhax.yaxis.set_major_locator(majorLocatorKM)
	self.kmhax.yaxis.set_minor_locator(minorLocator)

	self.fthax=self.kmhax.twinx()
	self.fthax.xaxis.set_ticks([],[])
	self.fthax.tick_params(axis='y', colors='k',labelsize=10)
	self.fthax.set_ylim(0,53.084)
	self.fthax.get_yaxis().set_tick_params(which="both",direction='out')
	self.fthax.yaxis.set_major_locator(majorLocatorKFT)
	self.fthax.yaxis.set_minor_locator(minorLocator)

    def readfile(self,fname):
	#--------------------------------------------------------------------
	# This *should* be a convenient way to read a uwyo sounding
	#--------------------------------------------------------------------
	fid=open(fname)
	lines=fid.readlines()
	nlines=len(lines)
	ndata=nlines-34
	output={}

	fields=lines[3].split()
	units=lines[4].split()

	# First line for WRF profiles differs from the UWYO soundings
	header=lines[0]
	if header[:5]=='00000':
	    # WRF profile
	    self['StationNumber']='-99999'
	    self['Longitude']=float(header.split()[5].strip(","))
	    self['Latitude']=float(header.split()[6])
	    self['SoundingDate']=header.split()[-1]
	else:
	    self['StationNumber']=header[:5]
	    dstr=(' ').join(header.split()[-4:])
	    self['SoundingDate']=datetime.strptime(dstr,"%HZ %d %b %Y").strftime("%Y-%m-%d_%H:%M:%S") 

	for ff in fields:
	    output[ff.lower()]=zeros((nlines-34))-999.

	lhi=[1, 9,16,23,30,37,46,53,58,65,72]
	rhi=[7,14,21,28,35,42,49,56,63,70,77]

	lcounter=5
	for line,idx in zip(lines[6:],range(ndata)):
	    lcounter+=1

	    try: output[fields[0].lower()][idx]=float(line[lhi[0]:rhi[0]])
	    except ValueError: break
	    
	    for ii in range(1,len(rhi)):
		try: 
		    # Debug only:
		    # print fields[ii].lower(), float(line[lhi[ii]:rhi[ii]].strip())
		    output[fields[ii].lower()][idx]=float(line[lhi[ii]:rhi[ii]].strip())
		except ValueError: 
		    pass

	for field in fields:
	    ff=field.lower()
	    self.data[ff]=ma.masked_values(output[ff],-999.)

	return None

    def lift_parcel(self,startp,startt,startdp,cape=None,ncape=None,**kwargs):
	"""Lift a parcel to discover certain properties.
	
	INPUTS:
	startp:  Pressure (hPa)
	startt:  Temperature (C)
	startdp: Dew Point Temperature (C)
	"""
	from numpy import interp

	#print startp, startt, startdp
	assert startt>startdp,"Not a valid parcel. Check Td<Tc"
	Pres=linspace(startp,100,100)

	# Lift the dry parcel
	T_dry=(startt+273.15)*(Pres/startp)**(Rs_da/Cp_da)-273.15 

#	print T_dry

	# Mixing ratio isopleth
	starte=SatVap(startdp)
	startw=MixRatio(starte,startp*100)
	e=Pres*startw/(.622+startw)
	T_iso=243.5/(17.67/log(e/6.112)-1)

	# Solve for the intersection of these lines (LCL).
	# interp requires the x argument (argument 2)
	# to be ascending in order!
	P_lcl=interp(0,T_iso-T_dry,Pres)
	T_lcl=interp(P_lcl,Pres[::-1],T_dry[::-1])

	col=[.6,.6,.6]

	# Plot traces below LCL
	self.skewxaxis.plot(T_dry[Pres>=P_lcl],Pres[Pres>=P_lcl],color=col,**kwargs)
	self.skewxaxis.plot(T_iso[Pres>=P_lcl],Pres[Pres>=P_lcl],color=col,**kwargs)
	self.skewxaxis.plot(T_lcl,P_lcl,ls='',marker='o',mec=col,mfc=col)

	# Now lift a wet parcel from the intersection point
	preswet=linspace(P_lcl,100)
	tempwet=lift_wet(T_lcl,preswet)

	# Plot trace above LCL
	self.skewxaxis.plot(tempwet,preswet,color=col,**kwargs)

	# Add text to sounding
	dtext ="Parcel:\n"
	dtext+="Ps:  %6.1fhPa\n"%startp
	dtext+="Ts:    %4.1fC\n"%startt
	dtext+="Ds:    %4.1fC\n"%startdp
	dtext+="Plcl:%6.1fhPa\n"%P_lcl
	dtext+="Tlcl:  %4.1fC\n"%T_lcl
        if cape is not None:
            dtext+="CAPE: %0.0f J kg$^{-1}$\n"%cape
        if ncape is not None:
            dtext+="NCAPE: %0.2f m s$^{-2}$\n"%ncape

        # instead of putting text on figure here, could return dtext, and annotate within
        # the plot_skewt() funtion. There, you could add CAPE and NCAPE to dtext
	self.fig.text(0.1,0.895,dtext,fontname="monospace",va='top',size=12)

	return Pres, T_dry, T_iso, preswet, tempwet

    def surface_parcel(self,mixdepth=100, pres_s = 1000.0): # added pres_s so can just change for diff regions
	"""Returns parameters for a parcel initialised by:
	1. Surface pressure (i.e. pressure of lowest level)
	2. Surface temperature determined from max(theta) of lowest <mixdepth> mbar
	3. Dew point temperature representative of lowest <mixdepth> mbar

	Inputs:
	mixdepth (mbar): depth to average mixing ratio over
	"""

#	print 'mixdepth: ', mixdepth
	pres=self.data["pres"]
	temp=self.data["temp"]
	dwpt=self.data["dwpt"]

	# parcel pressure is surface pressure
	# pres_s=pres[0]

	# identify the layers for averaging
	# layers=pres>pres[0]-mixdepth
	# BF ABOVE IS OLD WAY, DOING IT NEW WAY
	
	#layers = pres > (pres_s - mixdepth)
        # BB problem with above is if model data
        # goes down to 1000 mb, you're averaging (see below) over all data down 
        # to 1000 mb, even if surface pressure is higher (say 850 mb for Colorado)
        # BB added this line below to get rid of layers below surface (if there are some)
        # want to average over all layers from top of mixed layer to surface
        mdtop = pres_s - mixdepth # top of mixed layer (hPa)
	layers = (pres<=pres_s) & (pres >= mdtop) # layers above surface but below top of mixed layer
        
	# get max theta over mixheight to give
	# parcel temperature
        # BB why max? Should we use MEAN theta here??
	thta_mix=Theta(temp[layers]+273.15,pres[layers]*100.).max()
        # TempK inverts the theta function to get actual temperature
	temp_s=TempK(thta_mix,pres_s*100)-273.15

	# average mixing ratio over mixheight
	vpres=SatVap(dwpt) # get vapor pressure at all levels from dwpt
	mixr=MixRatio(vpres,pres*100) # get mixing ratio from vapor pressure
	mixr_mix=mixr[layers].mean() # averaged mix-layer mixing ratio
	vpres_s=MixR2VaporPress(mixr_mix,pres_s*100) # surface mixing ratio

	# surface dew point temp
        # empirical formala to get dewpoint from vapor pressure (Bolton, 1980)
	dwpt_s=VapPres2Dwpt(vpres_s)

	return pres_s,temp_s,dwpt_s

def lift_wet(startt,pres):
    #--------------------------------------------------------------------
    # Lift a parcel moist adiabatically from startp to endp.
    # Init temp is startt in C, pressure levels are in hPa    
    #--------------------------------------------------------------------

    temp=startt
    t_out=zeros(pres.shape);t_out[0]=startt
    for ii in range(pres.shape[0]-1):
	delp=pres[ii]-pres[ii+1]
 	temp=temp-100*delp*GammaW(temp+273.15,(pres[ii]-delp/2)*100)
	t_out[ii+1]=temp

    return t_out


if __name__=='__main__':

    sounding=Sounding("../examples/94975.2013070900.txt")
    sounding.make_skewt_axes()
    sounding.add_profile(color='r',lw=2)
    sounding.lift_parcel(1033.,10.7,-0.9)

    sounding.fig.savefig("../examples/94975.2013070900.png")
    show()
    

