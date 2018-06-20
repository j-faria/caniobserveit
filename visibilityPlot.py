from __future__ import print_function
import numpy as np
from PyAstronomy import pyasl
from astropy.coordinates import SkyCoord
from astropy.coordinates import name_resolve
import ephem
import argparse
import matplotlib
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import MultipleLocator
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
rcParams['xtick.major.pad'] = 12
font0 = FontProperties()
font1 = font0.copy()
font0.set_family('sans-serif')
font0.set_weight('light')
font1.set_family('sans-serif')
font1.set_weight('medium')
from mpld3 import fig_to_html, plugins

jdbinsize = 1.0/24./20.


def _prepare_data(target):
    try:
        target = {'name': target, 'coord': SkyCoord.from_name(target)}
        return target
    except name_resolve.NameResolveError as e:
        return None


def _getJDarrays(date):
    jd = pyasl.jdcnv(date)
    jd_start = pyasl.jdcnv(date)-0.5
    jd_end = pyasl.jdcnv(date)+0.5
    jds = np.arange(jd_start, jd_end, jdbinsize)
    jdsub = jds - np.floor(jds[0])
    return jd, jds, jdsub


def _getAltAz(jds, ra, dec, obs):
    return pyasl.eq2hor(jds, np.ones(jds.size)*ra, np.ones(jds.size)*dec,
                        lon=obs['longitude'], lat=obs['latitude'], alt=obs['altitude'])


def visibilityPlot(date=None, target=None, observatory=None, plotLegend=False, showMoonDist=False, print2file=False, p=False):
  """
    Plot the visibility of target.
    Parameters
    ----------
    date: datetime
        The date for which to calculate the visibility.
    target: string
        Target to find visbility plot for.
    observatory: string
        Name of the observatory that pyasl.observatory can resolve.
        Basically, any of pyasl.listObservatories().keys()
    plotLegend: boolean, optional
        If True (default), show a legend.
    showMoonDist : boolean, optional
        If True (default), the Moon distance will be shown.
    print2file : boolean, optional
        If True, the plot will be dumped to a png-file.
  """
  target = _prepare_data(target)
  if target is None:
      return None

  if isinstance(observatory, dict):
    obs = observatory
  else:
    obs = pyasl.observatory(observatory)

  fig = plt.figure(figsize=(15, 10))
  fig.subplots_adjust(left=0.07, right=0.8, bottom=0.15, top=0.88)
  ax = host_subplot(111)

  target_coord = target['coord']
  target_ra = target_coord.ra.deg
  target_dec = target_coord.dec.deg

  # JD array
  jd, jds, jdsub = _getJDarrays(date)

  # Get alt/az of object
  altaz = _getAltAz(jds, target_ra, target_dec, obs)

  # Get alt/az of Sun
  sun_position = pyasl.sunpos(jd)
  sun_ra, sun_dec = sun_position[1], sun_position[2]
  sunpos_altaz = _getAltAz(jds, sun_ra, sun_dec, obs)

  # Define plot label
  plabel = "[{0:2d}]  {1!s}".format(1, target['name'])

  # Find periods of: day, twilight, and night
  day = np.where( sunpos_altaz[0] >= 0. )[0]
  twi = np.where( np.logical_and(sunpos_altaz[0] > -18., sunpos_altaz[0] < 0.) )[0]
  night = np.where( sunpos_altaz[0] <= -18. )[0]

  if (len(day) == 0) and (len(twi) == 0) and (len(night) == 0):
    print("\nVisibilityPlot - no points to draw\n")


  if showMoonDist:
    mpos = pyasl.moonpos(jds)
    mdist = pyasl.getAngDist(mpos[0], mpos[1], np.ones(jds.size)*target_ra,
                             np.ones(jds.size)*target_dec)
    bindist = int((2.0/24.)/jdbinsize)
    firstbin = np.random.randint(0,bindist)
    for mp in range(0, int(len(jds)/bindist)):
      bind = firstbin+mp*bindist
      if altaz[0][bind]-1. < 5.:
          continue
      ax.text(jdsub[bind], altaz[0][bind]-1., str(int(mdist[bind]))+r"$^\circ$", ha="center", va="top", \
              fontsize=8, stretch='ultra-condensed', fontproperties=font0, alpha=1.)


  if len(twi) > 1:
    # There are points in twilight
    linebreak = np.where((jdsub[twi][1:]-jdsub[twi][:-1]) > 2.0*jdbinsize)[0]
    if len(linebreak) > 0:
      plotrjd = np.insert(jdsub[twi], linebreak+1, np.nan)
      plotdat = np.insert(altaz[0][twi], linebreak+1, np.nan)
      ax.plot( plotrjd, plotdat, "-", color='#BEBEBE', linewidth=1.5)
    else:
      ax.plot( jdsub[twi], altaz[0][twi], "-", color='#BEBEBE', linewidth=1.5)

  ax.plot( jdsub[night], altaz[0][night], '.k', label=plabel)
  ax.plot( jdsub[day], altaz[0][day], '.', color='#FDB813')

  altmax = np.argmax(altaz[0])
  ax.text( jdsub[altmax], altaz[0][altmax], str(1), color="b", fontsize=14, \
           fontproperties=font1, va="bottom", ha="center")

  if 1 == 29:
    ax.text( 1.1, 1.0-float(1)*0.04, "too many targets", ha="left", va="top", transform=ax.transAxes, \
            fontsize=10, fontproperties=font0, color="r")
  else:
    ax.text( 1.1, 1.0-float(1)*0.04, plabel, ha="left", va="top", transform=ax.transAxes, \
            fontsize=12, fontproperties=font0, color="b")

  ax.text( 1.1, 1.03, "List of targets", ha="left", va="top", transform=ax.transAxes, \
          fontsize=12, fontproperties=font0, color="b")

  axrange = ax.get_xlim()
  ax.set_xlabel("UT [hours]")

  if axrange[1]-axrange[0] <= 1.0:
    jdhours = np.arange(0,3,1.0/24.)
    utchours = (np.arange(0,72,dtype=int)+12)%24
  else:
    jdhours = np.arange(0,3,1.0/12.)
    utchours = (np.arange(0,72, 2, dtype=int)+12)%24
  ax.set_xticks(jdhours)
  ax.set_xlim(axrange)
  ax.set_xticklabels(utchours, fontsize=18)

  # Make ax2 responsible for "top" axis and "right" axis
  ax2 = ax.twin()
  # Set upper x ticks
  ax2.set_xticks(jdhours)
  ax2.set_xticklabels(utchours, fontsize=18)
  ax2.set_xlabel("UT [hours]")

  # Horizon angle for airmass
  airmass_ang = np.arange(5.,90.,5.)
  geo_airmass = pyasl.airmass.airmassPP(90.-airmass_ang)
  ax2.set_yticks(airmass_ang)
  airmassformat = []
  for t in range(geo_airmass.size):
    airmassformat.append("{0:2.2f}".format(geo_airmass[t]))
  ax2.set_yticklabels(airmassformat, rotation=90)
  ax2.set_ylabel("Relative airmass", labelpad=32)
  ax2.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.015,-0.04, "Plane-parallel", transform=ax.transAxes, ha='left', \
           va='top', fontsize=10, rotation=90)

  ax22 = ax.twin()
  ax22.set_xticklabels([])
  ax22.set_frame_on(True)
  ax22.patch.set_visible(False)
  ax22.yaxis.set_ticks_position('right')
  ax22.yaxis.set_label_position('right')
  ax22.spines['right'].set_position(('outward', 25))
  ax22.spines['right'].set_color('k')
  ax22.spines['right'].set_visible(True)
  airmass2 = list(map(lambda ang: pyasl.airmass.airmassSpherical(90. - ang, obs['altitude']), airmass_ang))
  ax22.set_yticks(airmass_ang)
  airmassformat = []
  for t in airmass2:
      airmassformat.append("{0:2.2f}".format(t))
  ax22.set_yticklabels(airmassformat, rotation=90)
  ax22.tick_params(axis="y", pad=10, labelsize=10)
  plt.text(1.045,-0.04, "Spherical+Alt", transform=ax.transAxes, ha='left', va='top', \
           fontsize=10, rotation=90)

  ax3 = ax.twiny()
  ax3.set_frame_on(True)
  ax3.patch.set_visible(False)
  ax3.xaxis.set_ticks_position('bottom')
  ax3.xaxis.set_label_position('bottom')
  ax3.spines['bottom'].set_position(('outward', 50))
  ax3.spines['bottom'].set_color('k')
  ax3.spines['bottom'].set_visible(True)

  ltime, ldiff = pyasl.localtime.localTime(utchours, np.repeat(obs['longitude'], len(utchours)))
  jdltime = jdhours - ldiff/24.
  ax3.set_xticks(jdltime)
  ax3.set_xticklabels(utchours)
  ax3.set_xlim([axrange[0],axrange[1]])
  ax3.set_xlabel("Local time [hours]")

  ax.set_ylim([0, 91])
  ax.yaxis.set_major_locator(MultipleLocator(15))
  ax.yaxis.set_minor_locator(MultipleLocator(5))
  yticks = ax.get_yticks()
  ytickformat = []
  for t in range(yticks.size):
      ytickformat.append(str(int(yticks[t]))+r"$^\circ$")

  #TODO: This is where it breaks
  print('lalala')
  print(yticks)
  print(ytickformat)
  # ax.set_yticklabels(ytickformat, fontsize=16)
  ax.set_ylabel("Altitude", fontsize=18)
  yticksminor = ax.get_yticks(minor=True)
  ymind = np.where( yticksminor % 15. != 0. )[0]
  yticksminor = yticksminor[ymind]
  ax.set_yticks(yticksminor, minor=True)
  m_ytickformat = []
  for t in range(yticksminor.size):
      m_ytickformat.append(str(int(yticksminor[t]))+r"$^\circ$")
  ax.set_yticklabels(m_ytickformat, minor=True)
  ax.set_ylim([0, 91])

  ax.yaxis.grid(color='gray', linestyle='dashed')
  ax.yaxis.grid(color='gray', which="minor", linestyle='dotted')
  ax2.xaxis.grid(color='gray', linestyle='dotted')

  plt.text(0.5,0.95,"Visibility on {0!s}".format(date.date()), \
           transform=fig.transFigure, ha='center', va='bottom', fontsize=20)

  if plotLegend:
    line1 = matplotlib.lines.Line2D((0,0),(1,1), color='#FDB813', linestyle="-", linewidth=2)
    line2 = matplotlib.lines.Line2D((0,0),(1,1), color='#BEBEBE', linestyle="-", linewidth=2)
    line3 = matplotlib.lines.Line2D((0,0),(1,1), color='k', linestyle="-", linewidth=2)

    lgd2 = plt.legend((line1,line2,line3),("day","twilight","night",), \
                        bbox_to_anchor=(0.88, 0.13), loc='best', borderaxespad=0.,prop={'size':12}, fancybox=True)
    lgd2.get_frame().set_alpha(.5)

  obsco = r"Obs coord.: {0:8.4f}$^\circ$, {1:8.4f}$^\circ$, {2:4f} m".format(obs['longitude'], obs['latitude'], obs['altitude'])

  plt.text(0.01,0.97, obsco, transform=fig.transFigure, ha='left', va='center', fontsize=10)
  plt.text(0.01,0.95, obs['name'], transform=fig.transFigure, ha='left', va='center', fontsize=10)

  if p:
      plt.show()


  # fig = plt.figure(figsize=(10, 5))
  # ax = fig.add_subplot(111)
  # x = np.linspace(0, 2*np.pi, 100)
  # ax.plot(x, np.sin(x))
  # ax.plot(x, np.cos(x))


  return fig_to_html(fig)
