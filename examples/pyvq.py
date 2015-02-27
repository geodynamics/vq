#!/usr/bin/env python

from __future__ import print_function

import math
import sys
import argparse
import quakelib
import gc

scipy_available = True
try:
    import scipy.stats
except ImportError:
    scipy_available = False

matplotlib_available = True
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    import matplotlib.font_manager as mfont
    import matplotlib.colors as mcolor
    import matplotlib.colorbar as mcolorbar
    import matplotlib.lines as mlines
    from PIL import Image #TODO: Move this guy
    plt.switch_backend('agg') #Required for map plots
except ImportError:
    matplotlib_available = False

numpy_available = True
try:
    import numpy as np
except ImportError:
    numpy_available = False
    
MIN_LON_DIFF = 0.118   #corresponds to 10km at lat,lon = (40.35, -124.85)
MIN_LAT_DIFF = 0.090   #corresponds to 10km at lat,lon = (40.35, -124.85)

#-------------------------------------------------------------------------------
# Given a set of maxes and mins return a linear value betweem them.
# Used to compute cutoff for field value evaluation, cutoff scales with
# number of involved elements for FieldPlotter instances.
#-------------------------------------------------------------------------------
def linear_interp(x, x_min, x_max, y_min, y_max):
    return ((y_max - y_min)/(x_max - x_min) * (x - x_min)) + y_min

class SaveFile:
    def __init__(self, event_file, plot_type):
        self.filename = plot_type+"_"+event_file.split(".")[0]+".png"

class MagFilter:
    def __init__(self, min_mag=None, max_mag=None):
        self._min_mag = min_mag if min_mag is not None else -float("inf")
        self._max_mag = max_mag if max_mag is not None else float("inf")

    def test_event(self, event):
        return (event.getMagnitude() >= self._min_mag and event.getMagnitude() <= self._max_mag)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_mag != -float("inf"): label_str += str(self._min_mag)+"<"
        if self._max_mag != float("inf"): label_str += "M<"+str(self._max_mag)
        return label_str

class YearFilter:
    def __init__(self, min_year=None, max_year=None):
        self._min_year = min_year if min_year is not None else -float("inf")
        self._max_year = max_year if max_year is not None else float("inf")

    def test_event(self, event):
        return (event.getEventYear() >= self._min_year and event.getEventYear() <= self._max_year)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_year != -float("inf"): label_str += str(self._min_year)+"<"
        if self._max_year != float("inf"): label_str += "year<"+str(self._max_year)
        return label_str

class EventNumFilter:
    def __init__(self, min_event_num=None, max_event_num=None):
        self._min_event_num = min_event_num if min_event_num is not None else -sys.maxint
        self._max_event_num = max_event_num if max_event_num is not None else sys.maxint

    def test_event(self, event):
        return (event.getEventNumber() >= self._min_event_num and event.getEventNumber() <= self._max_event_num)

    def plot_str(self):
        label_str = ""
# TODO: change to <= character
        if self._min_event_num != -sys.maxint: label_str += str(self._min_event_num)+"<"
        if self._max_event_num != sys.maxint: label_str += "event num<"+str(self._max_event_num)
        return label_str

class SectionFilter:
    def __init__(self, model, section_list):
        self._section_list = section_list
        self._elem_to_section_map = {elem_num: model.element(elem_num).section_id() for elem_num in range(model.num_elements())}

    def test_event(self, event):
        event_elements = event.getInvolvedElements()
        for elem_num in event_elements:
            elem_section = self._elem_to_section_map[elem_num]
            if not elem_section in self._section_list: return True

        return False

    def plot_str(self):
        return "my_string"

class Events:
    def __init__(self, event_file, event_file_type, sweep_file = None):
        self._events = quakelib.ModelEventSet()
        if event_file_type == "hdf5":
            self._events.read_file_hdf5(event_file)
        elif event_file_type == "text" and sweep_file != None:
            self._events.read_file_ascii(event_file, sweep_file)
        else:
            raise "event_file_type must be hdf5 or text. If text, a sweep_file is required."
            
        self._filtered_events = range(len(self._events))
        self._plot_str = ""

    def plot_str(self):
        return self._plot_str

    def set_filters(self, filter_list):
        self._filtered_events = [evnum for evnum in range(len(self._events))]
        self._plot_str = ""
        for cur_filter in filter_list:
            new_filtered_events = [evnum for evnum in self._filtered_events if cur_filter.test_event(self._events[evnum])]
            self._filtered_events = new_filtered_events
            self._plot_str += cur_filter.plot_str()
        if len(self._filtered_events) == 0:
            raise "No events matching filters found!"

    def interevent_times(self):
        event_times = [self._events[evnum].getEventYear() for evnum in self._filtered_events]
        return [event_times[i+1]-event_times[i] for i in xrange(len(event_times)-1)]

    def event_years(self):
        return [self._events[evnum].getEventYear() for evnum in self._filtered_events]

    def event_rupture_areas(self):
        return [self._events[evnum].calcEventRuptureArea() for evnum in self._filtered_events]

    def event_magnitudes(self):
        return [self._events[evnum].getMagnitude() for evnum in self._filtered_events if self._events[evnum].getMagnitude() != float("-inf")]
# TODO: Handle -infinity magnitudes on the C++ side

    def event_mean_slip(self):
        return [self._events[evnum].calcMeanSlip() for evnum in self._filtered_events]
        
class FieldPlotter:
    def __init__(self, model, element_slips=None, event_id=None, events=None, cbar_max=None):
#TODO: Set fonts and figure size explicitly
        self.max_lat, self.max_lon = 0.0,0.0
        plot_height = 768.0
        max_map_width = 690.0
        max_map_height = 658.0
        map_res     = 'i'
        padding     = 0.08
        map_proj = 'cyl'
        # Define how the cutoff value scales if it is not explitly set.
        # Cutoff is the max distance away from elements to compute
        # the field given in units of element length.
        self.cutoff_min_size = 20.0
        self.cutoff_min = 20.0
        self.cutoff_p2_size = 65.0
        self.cutoff_p2 = 90.0
        # Read elements and slips into the SlippedElementList
        self.elements = quakelib.SlippedElementList()
        #if event_id is not None and events is not None and element_slips is None:
            # TODO: Implement event element accessor
            # get event_slips and element_ids
            # self.element_slips = events.get_event_element_slips(event_id)
        if event_id is None and events is None and element_slips is not None:
            self.element_ids = element_slips.keys()
            self.element_slips = element_slips
        else:
            raise "Must specify event_id for event fields or element_slips (dictionary of slip indexed by element_id) for custom field."
        if len(self.element_slips) != len(self.element_ids):
            raise "Must specify slip for all elements."
        for ele_id in self.element_ids:
            new_ele = model.create_slipped_element(ele_id)
            new_ele.set_slip(self.element_slips[ele_id])
            self.elements.append(new_ele)
        self.slip_map = quakelib.SlipMap()
        self.slip_map.add_elements(self.elements)
        # Grab base Lat/Lon from fault model, used for lat/lon <-> xyz conversion
        base = model.get_base()
        self.base_lat = self.min_lat = base[0]
        self.base_lon = self.min_lon = base[1]
        self.min_lat, self.max_lat, self.min_lon, self.max_lon = model.get_latlon_bounds()
        # Expand lat/lon range in the case of plotting a few elements
        if len(self.element_ids) < 10:
            self.min_lat = self.min_lat - MIN_LAT_DIFF*10
            self.max_lat = self.max_lat + MIN_LAT_DIFF*10
            self.min_lon = self.min_lon - MIN_LON_DIFF*10
            self.max_lon = self.max_lon + MIN_LON_DIFF*10
        # Adjust bounds for good framing on plot
        lon_range = self.max_lon - self.min_lon
        lat_range = self.max_lat - self.min_lat
        max_range = max((lon_range, lat_range))
        self.min_lon = self.min_lon - lon_range*padding
        self.min_lat = self.min_lat - lat_range*padding
        self.max_lon = self.max_lon + lon_range*padding
        self.max_lat = self.max_lat + lat_range*padding
        self.lat0, self.lon0  = (self.max_lat+self.min_lat)/2.0, (self.max_lon+self.min_lon)/2.0
        self.llcrnrlat = self.min_lat
        self.llcrnrlon = self.min_lon
        self.urcrnrlat = self.max_lat
        self.urcrnrlon = self.max_lon
        # We need a map instance to calculate the aspect ratio
        map = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        # Using the aspect ratio (h/w) to find the actual map width and height in pixels
        if map.aspect > max_map_height/max_map_width:
            map_height = max_map_height
            map_width = max_map_height/map.aspect
        else:
            map_width = max_map_width
            map_height = max_map_width*map.aspect
        # A conversion instance for doing the lat-lon to x-y conversions
        base_lld = quakelib.LatLonDepth(self.base_lat, self.base_lon, 0.0)
        self.convert = quakelib.Conversion(base_lld)
        self.lons_1d = np.linspace(self.min_lon, self.max_lon, num=int(map_width))
        self.lats_1d = np.linspace(self.min_lat, self.max_lat, num=int(map_height))
        _lons_1d = quakelib.FloatList()
        _lats_1d = quakelib.FloatList()
        for lon in self.lons_1d:
            _lons_1d.append(lon)
        for lat in self.lats_1d:
            _lats_1d.append(lat)
        # Set up the points for field evaluation, convert to xyz basis
        self.grid_1d = self.convert.convertArray2xyz(_lats_1d,_lons_1d)
        # Convert the fault traces to lat-lon
#TODO: find/create fault trace accessor
#        fault_traces = !!!!!!!!!!!.get_fault_traces()
#        fault_traces_latlon = {}
#        for secid in fault_traces.iterkeys():
#            fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])
        self._plot_str = ""
        #-----------------------------------------------------------------------
        # Gravity map configuration  #TODO: Put in switches for field_type
        #-----------------------------------------------------------------------
        if cbar_max is None:
            cbar_max = 20
        self.dmc = {
            'font':               mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='normal'),
            'font_bold':          mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='bold'),
            'cmap':               plt.get_cmap('seismic'),
        #water
            'water_color':          '#4eacf4',
        #map boundaries
            'boundary_color':       '#000000',
            'boundary_width':       1.0,
            'coastline_color':      '#000000',
            'coastline_width':      1.0,
            'country_color':        '#000000',
            'country_width':        1.0,
            'state_color':          '#000000',
            'state_width':          1.0,
        #rivers
            'river_width':          0.25,
        #faults
            'fault_color':          '#000000',
            'event_fault_color':    '#ff0000',
            'fault_width':          0.5,
        #lat lon grid
            'grid_color':           '#000000',
            'grid_width':           0.0,
            'num_grid_lines':       5,
        #map props
            'map_resolution':       map_res,
            'map_projection':       map_proj,
            'plot_resolution':      72.0,
            'map_tick_color':       '#000000',
            'map_frame_color':      '#000000',
            'map_frame_width':      1,
        #map_fontsize = 12
            'map_fontsize':         16.0,   # 12   THIS IS BROKEN
            'arrow_inset':          14.0, # 10
            'arrow_fontsize':       12.0, #  9
        #cb_fontsize = 12
            'cb_fontsize':          20.0, # 12
            'cb_fontcolor':         '#000000',
            'cb_height':            20.0,
            'cb_margin_t':          2.0, # 10
         #min/max gravity change labels for colorbar (in microgals)
            'cbar_min':             -cbar_max,
            'cbar_max':             cbar_max
        }
        self.norm = None
        #-----------------------------------------------------------------------
        # m1, fig1 is the oceans and the continents. This will lie behind the
        # masked data image.
        #-----------------------------------------------------------------------
        self.m1 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        #-----------------------------------------------------------------------
        # m2, fig2 is the plotted deformation data.
        #-----------------------------------------------------------------------
        self.m2 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        #-----------------------------------------------------------------------
        # m3, fig3 is the ocean land mask.
        #-----------------------------------------------------------------------
        self.m3 = Basemap(
            llcrnrlon=self.llcrnrlon,
            llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon,
            urcrnrlat=self.urcrnrlat,
            lat_0=self.lat0, 
            lon_0=self.lon0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        
    def compute_field(self, field_type, cutoff=None):
        self.lame_lambda = 3.2e10
        self.lame_mu     = 3.0e10
        #-----------------------------------------------------------------------
        # If the cutoff is none (ie not explicitly set) calculate the cutoff for
        # this event. 
        #-----------------------------------------------------------------------
        num_involved_elements = float(len(self.element_slips.keys()))
        if cutoff is None:
            if  num_involved_elements >= self.cutoff_min_size:
                cutoff = linear_interp(
                    num_involved_elements,
                    self.cutoff_min_size,
                    self.cutoff_p2_size,
                    self.cutoff_min,
                    self.cutoff_p2
                    )
            else:
                cutoff = self.cutoff_min
    
        sys.stdout.write('{:0.2f} cutoff : '.format(cutoff))
        sys.stdout.flush()
        if field_type.lower() == "gravity":
            print("Computing gravity field...")
            self.field_1d = self.slip_map.gravity_changes(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
        elif field_type.lower() == "displacement":
            self.field_1d = self.slip_map.displacements(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
            print("Computing displacement field...")
            # Grab only dz
        #elif field_type.lower() == "insar":
            #print("Computing InSAR field...")
#TODO: Have quakelib compute fringes
        elif field_type.lower() == "dilat_gravity":
            print("Computing dilatational gravity field...")
            self.field_1d = self.slip_map.dilat_gravity_changes(self.grid_1d, self.lame_lambda, self.lame_mu, cutoff)
        #elif field_type.lower() == "potential":
#TODO: Add potential calculation to QuakeLibElement.cpp, add here
        # Reshape field
        self.field = np.array(self.field_1d).reshape((self.lats_1d.size,self.lons_1d.size))
        
    def plot_str(self):
        return self._plot_str
        
    def create_field_image(self, fringes=True):
        #-----------------------------------------------------------------------
        # Set all of the plotting properties
        #-----------------------------------------------------------------------
        cmap            = self.dmc['cmap']
        water_color     = self.dmc['water_color']
        boundary_color  = self.dmc['boundary_color']
        land_color      = cmap(0.5)
        plot_resolution = self.dmc['plot_resolution']
        #-----------------------------------------------------------------------
        # Set the map dimensions
        #-----------------------------------------------------------------------
        mw = self.lons_1d.size
        mh = self.lats_1d.size
        mwi = mw/plot_resolution
        mhi = mh/plot_resolution
        #-----------------------------------------------------------------------
        # Fig1 is the background land and ocean.
        #-----------------------------------------------------------------------
        fig1 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m1.ax = fig1.add_axes((0,0,1,1))
        self.m1.drawmapboundary(
            color=boundary_color,
            linewidth=0,
            fill_color=water_color
        )
        self.m1.fillcontinents(
            color=land_color,
            lake_color=water_color
        )
        #-----------------------------------------------------------------------
        # Fig2 is the deformations.
        #-----------------------------------------------------------------------
        fig2 = plt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m2.ax = fig2.add_axes((0,0,1,1))
        
        # make sure the values are located at the correct location on the map
        dG_transformed = self.m2.transform_scalar(self.field, self.lons_1d, self.lats_1d, self.lons_1d.size, self.lats_1d.size)
        
        if self.norm is None:
            #self.norm = mcolor.Normalize(vmin=np.amin(dG_transformed), vmax=np.amax(dG_transformed))
            # Changed units to microgals (multiply MKS unit by 10^8)
            self.norm = mcolor.Normalize(vmin=self.dmc['cbar_min'], vmax=self.dmc['cbar_max'])
        
        #self.m2.imshow(dG_transformed, cmap=cmap, norm=self.norm)
        # Changed units to microgals (multiply MKS unit by 10^8)
        self.m2.imshow(dG_transformed*float(pow(10,8)), cmap=cmap, norm=self.norm)
        #-----------------------------------------------------------------------
        # Composite fig 1 - 2 together
        #-----------------------------------------------------------------------
        # FIGURE 1 draw the renderer
        fig1.canvas.draw()
        
        # FIGURE 1 Get the RGBA buffer from the figure
        w,h = fig1.canvas.get_width_height()
        buf = np.fromstring ( fig1.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 1 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im1 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        # FIGURE 2 draw the renderer
        fig2.canvas.draw()
        
        # FIGURE 2 Get the RGBA buffer from the figure
        w,h = fig2.canvas.get_width_height()
        buf = np.fromstring ( fig2.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 2 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im2 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        # Clear all three figures
        fig1.clf()
        fig2.clf()
        plt.close('all')
        gc.collect()
        return im2

    def plot_field(self, field_type, output_file=None, fringes=True, hi_res=False, cutoff=None):
        self.compute_field(field_type, cutoff=cutoff)
        map_image = self.create_field_image(fringes=fringes)
        
        sys.stdout.write('map overlay : ')
        sys.stdout.flush()
        # Convert the fault traces to lat-lon
        """
        fault_traces_latlon = {}
        for secid in fault_traces.iterkeys():
             fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(self.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])
        """
        #---------------------------------------------------------------------------
        # Plot all of the geographic info on top of the displacement map image.
        #---------------------------------------------------------------------------
        # Grab all of the plot properties that we will need.
        # properties that are fringes dependent
        if fringes and field_type == 'displacement':
            cmap            = self.dmc['cmap_f']
            coastline_color = self.dmc['coastline_color_f']
            country_color   = self.dmc['country_color_f']
            state_color     = self.dmc['state_color_f']
            fault_color     = self.dmc['fault_color_f']
            map_tick_color  = self.dmc['map_tick_color_f']
            map_frame_color = self.dmc['map_frame_color_f']
            grid_color      = self.dmc['grid_color_f']
            cb_fontcolor    = self.dmc['cb_fontcolor_f']
        else:
            cmap            = self.dmc['cmap']
            coastline_color = self.dmc['coastline_color']
            country_color   = self.dmc['country_color']
            state_color     = self.dmc['state_color']
            fault_color     = self.dmc['fault_color']
            map_tick_color  = self.dmc['map_tick_color']
            map_frame_color = self.dmc['map_frame_color']
            grid_color      = self.dmc['grid_color']
            cb_fontcolor    = self.dmc['cb_fontcolor']

        # properties that are not fringes dependent
        boundary_width  = self.dmc['boundary_width']
        coastline_width = self.dmc['coastline_width']
        country_width   = self.dmc['country_width']
        state_width     = self.dmc['state_width']
        river_width     = self.dmc['river_width']
        fault_width     = self.dmc['fault_width']
        map_frame_width = self.dmc['map_frame_width']
        map_fontsize    = self.dmc['map_fontsize']
        arrow_inset     = self.dmc['arrow_inset']
        arrow_fontsize  = self.dmc['arrow_fontsize']
        cb_fontsize     = self.dmc['cb_fontsize']
        cb_height       = self.dmc['cb_height']
        cb_margin_t     = self.dmc['cb_margin_t']
        grid_width      = self.dmc['grid_width']
        num_grid_lines  = self.dmc['num_grid_lines']
        font            = self.dmc['font']
        font_bold       = self.dmc['font_bold']

        map_resolution  = self.dmc['map_resolution']
        map_projection  = self.dmc['map_projection']
        plot_resolution = self.dmc['plot_resolution']

        # The sizing for the image is tricky. The aspect ratio of the plot is fixed,
        # so we cant set all of margins to whatever we want. We will set the anchor
        # to the top, left margin position. Then scale the image based on the
        # bottom/right margin, whichever is bigger.

        mw = self.lons_1d.size
        mh = self.lats_1d.size

        if mh > mw:
            ph = 768.0
            pw = mw + 70.0 + 40.0
        else:
            pw = 790.0
            ph = mh + 70.0 + 40.0

        width_frac = mw/pw
        height_frac = mh/ph
        left_frac = 70.0/pw
        bottom_frac = 70.0/ph

        pwi = pw/plot_resolution
        phi = ph/plot_resolution

        if hi_res:
            fig_res = plot_resolution*4.0
        else:
            fig_res = plot_resolution

        fig4 = plt.figure(figsize=(pwi, phi), dpi=fig_res)

        #---------------------------------------------------------------------------
        # m4, fig4 is all of the boundary data.
        #---------------------------------------------------------------------------
        m4 = Basemap(
            llcrnrlon=self.min_lon,
            llcrnrlat=self.min_lat,
            urcrnrlon=self.max_lon,
            urcrnrlat=self.max_lat,
            lat_0=(self.max_lat+self.min_lat)/2.0,
            lon_0=(self.max_lon+self.min_lon)/2.0,
            resolution=map_resolution,
            projection=map_projection,
            suppress_ticks=True
        )
        m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))

        # draw coastlines, edge of map.
        m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)

        # draw countries
        m4.drawcountries(linewidth=country_width, color=country_color)

        # draw states
        m4.drawstates(linewidth=state_width, color=state_color)

        # draw parallels.
        parallels = np.linspace(self.lats_1d.min(), self.lats_1d.max(), num_grid_lines+1)
        m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10], fontsize=map_fontsize)

        # draw meridians
        meridians = np.linspace(self.lons_1d.min(), self.lons_1d.max(), num_grid_lines+1)
        m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10], fontsize=map_fontsize)

        if field_type == 'displacement':
            box_size = 70.0
            # draw the azimuth look arrow
            az_width_frac    = box_size/pw
            az_height_frac   = box_size/ph
            az_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
            az_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac)/ph
            az_ax = fig4.add_axes((az_left_frac,az_bottom_frac,az_width_frac,az_height_frac))

            az_ax.set_xlim((0,1.0))
            az_ax.set_ylim((0,1.0))
            for item in az_ax.yaxis.get_ticklabels() + az_ax.xaxis.get_ticklabels() + az_ax.yaxis.get_ticklines() + az_ax.xaxis.get_ticklines():
                item.set_alpha(0)

            az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(self.look_azimuth)
            az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(self.look_azimuth)
            az_arrow_dx      = 0.8*math.sin(self.look_azimuth)
            az_arrow_dy      = 0.8*math.cos(self.look_azimuth)

            az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
            az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
            az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - self.convert.rad2deg(self.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
            az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(self.convert.rad2deg(self.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

            # draw the altitude look arrow
            al_width_frac    = box_size/pw
            al_height_frac   = box_size/ph
            al_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
            al_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac - ph*al_height_frac)/ph
            al_ax = fig4.add_axes((al_left_frac,al_bottom_frac,al_width_frac,al_height_frac))

            al_ax.set_xlim((0,1.0))
            al_ax.set_ylim((0,1.0))
            for item in al_ax.yaxis.get_ticklabels() + al_ax.xaxis.get_ticklabels() + al_ax.yaxis.get_ticklines() + al_ax.xaxis.get_ticklines():
                item.set_alpha(0)

            al_arrow_start_x    = 0.1 + 0.8*math.cos(self.look_elevation)
            al_arrow_start_y    = 0.1 + 0.8*math.sin(self.look_elevation)
            al_arrow_dx      = -0.8*math.cos(self.look_elevation)
            al_arrow_dy      = -0.8*math.sin(self.look_elevation)

            al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
            al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
            al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=self.convert.rad2deg(self.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
            al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(self.convert.rad2deg(self.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
            
            # draw the box with the magnitude
            mag_width_frac    = box_size/pw
            mag_height_frac   = 15.0/ph    # originally 10.0/ph
            mag_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
            mag_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac  - ph*az_height_frac - ph*mag_height_frac)/ph
            mag_ax = fig4.add_axes((mag_left_frac,mag_bottom_frac,mag_width_frac,mag_height_frac))

            mag_ax.set_xlim((0,1.0))
            mag_ax.set_ylim((0,1.0))
            for item in mag_ax.yaxis.get_ticklabels() + mag_ax.xaxis.get_ticklabels() + mag_ax.yaxis.get_ticklines() + mag_ax.xaxis.get_ticklines():
                item.set_alpha(0)
            
            mag_ax.text(0.5, 0.5, 'm = {:0.3f}'.format(float(event_data['event_magnitude'])), fontproperties=font_bold, size=arrow_fontsize, ha='center', va='center')

        # add the map image to the plot
        m4.imshow(map_image, origin='upper')

        # print faults on lon-lat plot
        """
        for sid, sec_trace in fault_traces_latlon.iteritems():
            trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
            
            if sid in event_sections:
                linewidth = fault_width + 2.5
            else:
                linewidth = fault_width

            m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')
        """

        #plot the cb
        left_frac = 70.0/pw
        bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
        width_frac = mw/pw
        height_frac = cb_height/ph

        cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
        norm = self.norm
        cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
               norm=norm,
               orientation='horizontal')
        if field_type == 'displacement':
            if fringes:
                cb_title = 'Displacement [m]'
            else:
                cb_title = 'Total displacement [m]'

        elif field_type == 'gravity' or field_type == 'dilat_gravity':
            cb_title        = r'Gravity changes [$\mu gal$]'
            # Make first and last ticks on colorbar be <MIN and >MAX
            cb_tick_labs    = [item.get_text() for item in cb_ax.get_xticklabels()]
            cb_tick_labs[0] = '<'+cb_tick_labs[0]
            cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
            cb_ax.set_xticklabels(cb_tick_labs)

        cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(font)
            label.set_fontsize(cb_fontsize)
            label.set_color(cb_fontcolor)
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)

        fig4.savefig(output_file, format='png', dpi=fig_res)
        sys.stdout.write('\nPlot saved: {}'.format(output_file))
        sys.stdout.write('\ndone\n')
        sys.stdout.flush()

class BasePlotter:
    def create_plot(self, plot_type, log_y, x_data, y_data, plot_title, x_label, y_label, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        if plot_type == "scatter":
            ax.scatter(x_data, y_data)
        elif plot_type == "line":
            ax.plot(x_data, y_data)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def multi_line_plot(self, x_data, y_data, colors, labels, linewidths, plot_title, x_label, y_label, legend_str, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if not (len(x_data) == len(y_data) and len(x_data) == len(colors) and len(colors) == len(labels) and len(linewidths) == len(colors)):
            raise "These lists must be the same length: x_data, y_data, colors, labels, linewidths."
        for i in range(len(x_data)):
            ax.plot(x_data[i], y_data[i], color=colors[i], label=labels[i], linewidth=linewidths[i])
        ax.legend(title=legend_str, loc='best')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def t0_vs_dt_plot(self, t0_dt_plot, wait_75, filename):
# TODO: Set fonts explicitly
        t0_dt_main_line_color   = '#000000'
        t0_dt_sub_line_color    = '#737373'
        t0_dt_main_line_width   = 2.0
        t0_dt_sub_line_width    = 1.0
        t0_dt_range_color       = plt.get_cmap('autumn')(0.99)
        years_since_line_color  = 'blue'
        legend_loc              = 'best'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(r't$_0$ [years]')
        ax.set_ylabel(r'$\Delta$t [years]')
        percents = t0_dt_plot.keys()
        ax.fill_between(t0_dt_plot[min(percents)]['x'], t0_dt_plot[min(percents)]['y'], y2=t0_dt_plot[max(percents)]['y'], linewidth=0, facecolor=t0_dt_range_color)
        for percent in t0_dt_plot.iterkeys():
            if percent == min(percents):
                linewidth = t0_dt_sub_line_width
                color = t0_dt_sub_line_color
                linestyle = '--'
            elif percent == max(percents):
                linewidth = t0_dt_sub_line_width
                color = t0_dt_sub_line_color
                linestyle = ':'
            else:
                linewidth = t0_dt_main_line_width
                color = t0_dt_main_line_color
                linestyle = '-'
            ax.plot(t0_dt_plot[percent]['x'], t0_dt_plot[percent]['y'], color=color, linewidth=linewidth, linestyle=linestyle, label='{}%'.format(percent))
        if wait_75 is not None:
            # Draw vertical dotted line where "today" is denoted by years_since
            ax.axvline(x=years_since,ymin=0,ymax=wait_75,color=years_since_line_color,linewidth=t0_dt_main_line_width,linestyle='--')
        ax.legend(title='event prob.', loc=legend_loc, handlelength=5)
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def scatter_and_errorbar(self, log_y, x_data, y_data, err_x, err_y, y_error, plot_title, x_label, y_label, filename, add_x = None, add_y = None, add_label = None):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data)
        ax.errorbar(err_x, err_y, yerr = y_error, label="UCERF2", ecolor='r')
        if add_x is not None:
            if log_y: ax.semilogy(add_x, add_y, label = add_label, c = 'k')
            if not log_y: ax.plot(add_x, add_y, label = add_label, c = 'k')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        ax.legend(loc = "best")
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

    def scatter_and_line(self, log_y, x_data, y_data, line_x, line_y, line_label, plot_title, x_label, y_label, filename):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(plot_title)
        if log_y:
            ax.set_yscale('log')
        ax.scatter(x_data, y_data)
        ax.plot(line_x, line_y, label = line_label, c = 'k')
        plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
        ax.legend(loc = "best")
        plt.savefig(filename,dpi=100)
        sys.stdout.write("Plot saved: {}\n".format(filename))

class MagnitudeRuptureAreaPlot(BasePlotter):
    def plot(self, events, filename):
        ra_list = events.event_rupture_areas()
        mag_list = events.event_magnitudes()
        ra_renorm_list = [quakelib.Conversion().sqm2sqkm(ra) for ra in ra_list]
        self.create_plot("scatter", True, mag_list, ra_renorm_list, events.plot_str(), "Magnitude", "Rupture Area (square km)", filename)

class MagnitudeMeanSlipPlot(BasePlotter):
    def plot(self, events, filename):
        slip_list = events.event_mean_slip()
        mag_list = events.event_magnitudes()
        self.create_plot("scatter", True, mag_list, slip_list, events.plot_str(), "Magnitude", "Mean Slip (meters)", filename)

class FrequencyMagnitudePlot(BasePlotter):
    def plot(self, events, filename, UCERF2 = False, b1 = False):
        # California observed seismicity rates and errorbars (UCERF2)
        x_UCERF = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
        y_UCERF = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
        y_error_UCERF = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
        add_x, add_y, add_label = None, None, None
        mag_list = events.event_magnitudes()
        cum_freq = {}
        freq_x, freq_y = [], []
        num_events = len(mag_list)
        years = events.event_years()
        if min(years) < max(years)*.01:
            # In most sims, it takes 30-100 years for first sim to occur, but the sim started at year 0. So if the first event occurs at a year that's before the first 1% of sim time, consider the filtered events to represent year=0 and onward. Needed for accurate # of events/yr
            year_range = max(years)
        else:
            year_range = max(years) - min(years)
        for num, mag in enumerate(sorted(mag_list)):
            cum_freq[mag] = num_events - (num + 1)
        for mag in sorted(cum_freq.iterkeys()):
            freq_x.append(mag)
            freq_y.append(float(cum_freq[mag])/year_range)
        if b1:
            add_x = np.linspace(min(freq_x),max(freq_x),10)
            add_y = 10**(math.log(freq_y[0],10)+freq_x[0]-add_x)
            add_label = "b==1"
        if UCERF2:
            self.scatter_and_errorbar(True, freq_x, freq_y, x_UCERF, y_UCERF, y_error_UCERF, events.plot_str(), "Magnitude (M)", "log(# events with mag > M /year)", filename, add_x=add_x, add_y=add_y, add_label=add_label)
        if b1 and not UCERF2:
            self.scatter_and_line(True, freq_x, freq_y, add_x, add_y, add_label, events.plot_str(), "Magnitude (M)", "log(# events with mag > M /year)", filename)
        if not UCERF2 and not b1:
            self.create_plot("scatter", True, freq_x, freq_y, events.plot_str(), "Magnitude (M)", "log(# events with mag > M /year)", filename)

class StressHistoryPlot(BasePlotter):
    def plot(self, stress_set, elements):
        stress_histories = {}
        for element in elements:
            stress_histories[element] = []
        for stress_state in stress_set:
            if stress_state.getSweepNum() == 0:
                stresses = stress_state.stresses()
                for stress in stresses:
                    if stress._element_id in elements:
                        stress_histories[element].append(stress._shear_stress)
        for element in elements:
            print(stress_histories[element])
        #self.create_plot("scatter", True, mag_vals, mag_norm, events.plot_str(), "Shear Stress", "Year")

class ProbabilityPlot(BasePlotter):
    def plot_p_of_t(self, events, filename):
        # Cumulative probability P(t) as a function of interevent time t
        intervals = np.array(events.interevent_times())
        prob = {}
        prob['x'] = np.sort(intervals)
        prob['y'] = np.arange(float(intervals.size))/float(intervals.size)
        self.create_plot("line", False, prob['x'], prob['y'], events.plot_str(),"t [years]", "P(t)", filename)

    def plot_conditional_fixed_dt(self, events, filename, fixed_dt=30.0):
        # P(t0 + dt, t0) vs. t0 for fixed dt
        intervals = np.array(events.interevent_times())
        prob_dt = {'x':[],'y':[]}
        t0_to_eval = np.arange(0.0,int(intervals.max())+.01,1.0)
        for t0 in t0_to_eval:
            int_t0_dt = intervals[np.where( intervals > t0+fixed_dt)]
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                prob_dt['x'].append(t0)
                prob_dt['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))
        self.create_plot("line", False, prob_dt['x'], prob_dt['y'], events.plot_str(),"t0 [years]", "P(t0 + dt, t0)", filename)

    def plot_p_of_t_multi(self, events, filename, beta=None, tau=None, num_t0=4, numPoints=200):
        # Cumulative conditional probability P(t,t0) as a function of
        # interevent time t, computed for multiple t0. Beta/Tau are Weibull parameters
        line_colormap = plt.get_cmap('autumn')
        intervals = np.array(events.interevent_times())
        conditional = {}
        weibull = {}
        max_t0 = int(intervals.max())
        t0_to_eval = np.linspace(0, max_t0, num=numPoints)
        for t0 in t0_to_eval:
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                conditional[t0] = {'x':[],'y':[]}
                weibull[t0]     = {'x':[],'y':[]}
                for dt in range(max_t0-int(t0)):
                    int_t0_dt  = intervals[np.where( intervals > t0+dt)]
                    prob_t0_dt = 1.0 - float(int_t0_dt.size)/float(int_t0.size)
                    conditional[t0]['x'].append(t0+dt)
                    conditional[t0]['y'].append(prob_t0_dt)
                    if beta is not None and tau is not None:
                        weibull[t0]['x'].append(t0+dt)
                        weibull_t0_dt = Distributions().cond_weibull(weibull[t0]['x'][-1],t0,beta,tau)
                        weibull[t0]['y'].append(weibull_t0_dt)
            else:
                conditional[t0] = None
                weibull[t0] = None
        t0_to_plot  = np.linspace(0, int(max_t0/2.0), num=num_t0)
        match_inds  = [np.abs(np.array(t0_to_eval-t0_to_plot[i])).argmin() for i in range(len(t0_to_plot))]
        t0_to_plot  = np.array([t0_to_eval[k] for k in match_inds])
        x_data_prob = [conditional[t0]['x'] for t0 in t0_to_plot]
        y_data_prob = [conditional[t0]['y'] for t0 in t0_to_plot]
        t0_colors   = [line_colormap(float(t0*.8)/t0_to_plot.max()) for t0 in t0_to_plot]
        prob_lw     = [2 for t0 in t0_to_plot]
        if beta is not None and tau is not None:
            x_data_weib = [weibull[t0]['x'] for t0 in t0_to_plot]
            y_data_weib = [weibull[t0]['y'] for t0 in t0_to_plot]
            weib_colors = ['k' for t0 in t0_to_plot]
            weib_labels = [None for t0 in t0_to_plot]
            weib_lw     = [1 for t0 in t0_to_plot]
            # List concatenation, not addition
            colors = t0_colors + weib_colors
            x_data = x_data_prob + x_data_weib
            y_data = y_data_prob + y_data_weib
            labels = [round(t0_to_plot[k],2) for k in range(len(t0_to_plot))] + weib_labels
            linewidths = prob_lw + weib_lw
        else:
            colors = t0_colors
            x_data = x_data_prob
            y_data = y_data_prob
            labels = [round(t0_to_plot[k],1) for k in range(len(t0_to_plot))]
            linewidths = prob_lw
        legend_string = r't$_0$='
        y_lab         = r'P(t, t$_0$)'
        x_lab         = r't = t$_0$ + $\Delta$t [years]'
        plot_title    = ""
        self.multi_line_plot(x_data, y_data, colors, labels, linewidths, plot_title, x_lab, y_lab, legend_string, filename)

    def plot_dt_vs_t0(self, events, filename, years_since=None):
        # Plot the waiting times corresponding to 25/50/75% conditional probabilities
        # as a function of t0 (time since last earthquake on the selected faults).
        # years_since is the number of years since the last observed (real) earthquake
        # on the selected faults.
        intervals = np.array(events.interevent_times())
        conditional = {}
        wait_75 = None
        max_t0 = int(intervals.max())
        # t0_to_eval used to evaluate waiting times with 25/50/75% probability given t0=years_since
        t0_to_eval = np.arange(0, max_t0, 1.0)
        # t0_to_plot is "smoothed" so that the plots aren't as jagged
        t0_to_plot = np.linspace(0, int(max_t0/2.0), num=5)
        t0_to_plot = [int(t0) for t0 in t0_to_plot]
        t0_dt      = {}
        t0_dt_plot = {}
        # First generate the conditional distributions P(t,t0) for each t0
        for t0 in t0_to_eval:
            int_t0 = intervals[np.where( intervals > t0)]
            if int_t0.size != 0:
                conditional[t0] = {'x':[],'y':[]}
                for dt in range(max_t0-int(t0)):
                    int_t0_dt = intervals[np.where( intervals > t0+dt)]
                    conditional[t0]['x'].append(t0+dt)
                    prob_t0_dt    = 1.0 - float(int_t0_dt.size)/float(int_t0.size)
                    conditional[t0]['y'].append(prob_t0_dt)
        # Loop over the probabilities whose waiting times we want to plot, invert P(t,t0)
        for percent in [0.25, 0.5, 0.75]:
            t0_dt[int(percent*100)]      = {'x':[],'y':[]}
            t0_dt_plot[int(percent*100)] = {'x':[],'y':[]}
            for t0 in t0_to_eval:
                if conditional[t0] is not None:
                    # Invert the conditional probabilities, find the recurrence time closest to
                    # the current percent
                    index   = (np.abs(np.array(conditional[t0]['y'])-percent)).argmin()
                    dt      = conditional[t0]['x'][index]-t0
                    t0_dt[int(percent*100)]['x'].append(t0)
                    t0_dt[int(percent*100)]['y'].append(dt)
                    if t0 in t0_to_plot:
                        t0_dt_plot[int(percent*100)]['x'].append(t0)
                        t0_dt_plot[int(percent*100)]['y'].append(dt)
        if years_since is not None:
            # Print out the "Forecast", the 25/50/75% probability given t0=years_since
            ind_25 = (np.abs(np.array(t0_dt[25]['x'])-years_since)).argmin()
            ind_50 = (np.abs(np.array(t0_dt[50]['x'])-years_since)).argmin()
            ind_75 = (np.abs(np.array(t0_dt[75]['x'])-years_since)).argmin()
            wait_25 = t0_dt[25]['y'][ind_25]
            wait_50 = t0_dt[50]['y'][ind_50]
            wait_75 = t0_dt[75]['y'][ind_75]
            sys.stdout.write('For t0 = {:.2f} years'.format(year_eval))
            sys.stdout.write('\n25% waiting time: {:.2f} years'.format(wait_25))
            sys.stdout.write('\n50% waiting time: {:.2f} years'.format(wait_50))
            sys.stdout.write('\n75% waiting time: {:.2f} years'.format(wait_75))
            sys.stdout.write('\n=======================================\n\n')
        self.t0_vs_dt_plot(t0_dt_plot, wait_75, filename)

class Distributions:
    def weibull(self, X, beta, tau):
        # Return the Weibull distribution at a point
        return 1-np.exp( -(X/float(tau))**beta)

    def cond_weibull(self, X, t0, beta, tau):
        # Return the conditional Weibull distribution at a single point
        return 1-np.exp( (t0/float(tau))**beta - (X/float(tau))**beta)

if __name__ == "__main__":
    # Specify arguments
    parser = argparse.ArgumentParser(description="PyVQ.")

    # Event/model file arguments
    parser.add_argument('--event_file', required=True,
            help="Name of event file to analyze.")
    parser.add_argument('--event_file_type', required=True,
            help="Event file type, either hdf5 or text.")
    parser.add_argument('--sweep_file', required=False,
            help="Name of sweep file to analyze.")
    parser.add_argument('--model_file', required=False,
            help="Name of model (geometry) file to use in analysis.")
    parser.add_argument('--stress_index_file', required=False,
            help="Name of stress index file to use in analysis.")
    parser.add_argument('--stress_file', required=False,
            help="Name of stress file to use in analysis.")

    # Event filtering arguments
    parser.add_argument('--min_magnitude', type=float, required=False,
            help="Minimum magnitude of events to process.")
    parser.add_argument('--max_magnitude', type=float, required=False,
            help="Maximum magnitude of events to process.")
    parser.add_argument('--min_year', type=float, required=False,
            help="Minimum year of events to process.")
    parser.add_argument('--max_year', type=float, required=False,
            help="Maximum year of events to process.")
    parser.add_argument('--min_event_num', type=float, required=False,
            help="Minimum event number of events to process.")
    parser.add_argument('--max_event_num', type=float, required=False,
            help="Maximum event number of events to process.")
    parser.add_argument('--use_sections', type=int, nargs='+', required=False,
            help="List of model sections to use (all sections used if unspecified).")

    # Event plotting arguments
    parser.add_argument('--plot_freq_mag', required=False, type=int,
            help="Generate frequency magnitude plot. 1: Only event data, 2: Plot b=1 Gutenberg-Richter relation, 3: Plot UCERF2 observed seismicity rates, 4: Plot UCERF2 and the b=1 line.")
    parser.add_argument('--plot_mag_rupt_area', required=False, action='store_true',
            help="Generate magnitude vs rupture area plot.")
    parser.add_argument('--plot_mag_mean_slip', required=False, action='store_true',
            help="Generate magnitude vs mean slip plot.")

    # Probability plotting arguments
    parser.add_argument('--plot_prob_vs_t', required=False, action='store_true',
            help="Generate earthquake recurrence probability at time t plot.")
    parser.add_argument('--plot_prob_vs_t_fixed_dt', required=False, action='store_true',
            help="Generate earthquake recurrence probability at time t + dt vs t plot.")
    parser.add_argument('--plot_cond_prob_vs_t', required=False, action='store_true',
            help="Generate earthquake recurrence conditional probabilities at time t = t0 + dt for multiple t0.")
    parser.add_argument('--plot_waiting_times', required=False, action='store_true',
            help="Generate waiting times until the next earthquake as function of time since last earthquake.")
    parser.add_argument('--beta', required=False, type=float,
            help="Beta parameter for the Weibull distribution, must also specify Tau")
    parser.add_argument('--tau', required=False, type=float,
            help="Tau parameter for the Weibull distribution, must also specify Beta")
            
    # Field plotting arguments
    parser.add_argument('--field_plot', required=False, action='store_true',
            help="Plot surface field for a specified event, e.g. gravity changes or displacements.")
    parser.add_argument('--field_type', required=False, help="Field type: gravity, dilat_gravity, displacement, insar")
    parser.add_argument('--field_savefile', required=False, help="File name to save ")
    parser.add_argument('--colorbar_max', required=False, type=float, help="Max unit for colorbar")

    # Stress plotting arguments
    parser.add_argument('--stress_elements', type=int, nargs='+', required=False,
            help="List of elements to plot stress history for.")

    # Validation/testing arguments
    parser.add_argument('--validate_slip_sum', required=False,
            help="Ensure the sum of mean slip for all events is within 1 percent of the specified value.")
    parser.add_argument('--validate_mean_interevent', required=False,
            help="Ensure the mean interevent time for all events is within 2 percent of the specified value.")

    args = parser.parse_args()

    # Read the event and sweeps files
    events = Events(args.event_file, args.event_file_type, args.sweep_file)

    # Read the geometry model if specified
    if args.model_file:
        model = quakelib.ModelWorld()
        model.read_file_ascii(args.model_file)
# TODO: add HDF5 compatibility
    else:
        if args.use_sections:
            raise "Model file required if specifying fault sections."
        else:
            model = None

    # Check that if either beta or tau is given then the other is also given
    if (args.beta and not args.tau) or (args.tau and not args.beta):
        raise "Must specify both beta and tau."
        
    # Check that field_type is one of the supported types
    if args.field_type:
        type = args.field_type.lower()
        if type != "gravity" and type != "dilat_gravity" and type != "displacement" and type != "insar":
            raise "Field type is one of gravity, dilat_gravity, displacement, insar"

    # Read the stress files if specified
    if args.stress_index_file and args.stress_file:
        stress_set = quakelib.ModelStressSet()
        stress_set.read_file_ascii(args.stress_index_file, args.stress_file)
    else:
        stress_set = None

    # Set up filters
    event_filters = []
    if args.min_magnitude or args.max_magnitude:
        event_filters.append(MagFilter(min_mag=args.min_magnitude, max_mag=args.max_magnitude))

    if args.min_year or args.max_year:
        event_filters.append(YearFilter(min_mag=args.min_year, max_mag=args.max_year))

    if args.min_event_num or args.max_event_num:
        event_filters.append(EventNumFilter(min_mag=args.min_event_num, max_mag=args.max_event_num))

    if args.use_sections:
        event_filters.append(SectionFilter(model, args.use_sections))

    events.set_filters(event_filters)

    # Generate plots
    if args.plot_freq_mag:
        filename = SaveFile(args.event_file, "freq_mag").filename
        if args.plot_freq_mag == 1: UCERF2,b1 = False, False
        if args.plot_freq_mag == 2: UCERF2,b1 = False, True
        if args.plot_freq_mag == 3: UCERF2,b1 = True, False
        if args.plot_freq_mag == 4: UCERF2,b1 = True, True
        FrequencyMagnitudePlot().plot(events, filename, UCERF2=UCERF2, b1=b1)
    if args.plot_mag_rupt_area:
        filename = SaveFile(args.event_file, "mag_rupt_area").filename
        MagnitudeRuptureAreaPlot().plot(events, filename)
    if args.plot_mag_mean_slip:
        filename = SaveFile(args.event_file, "mag_mean_slip").filename
        MagnitudeMeanSlipPlot().plot(events, filename)
    if args.plot_prob_vs_t:
        filename = SaveFile(args.event_file, "prob_vs_time").filename
        ProbabilityPlot().plot_p_of_t(events, filename)
    if args.plot_prob_vs_t_fixed_dt:
        filename = SaveFile(args.event_file, "p_vs_t_fixed_dt").filename
        ProbabilityPlot().plot_conditional_fixed_dt(events, filename)
    if args.plot_cond_prob_vs_t:
        filename = SaveFile(args.event_file, "cond_prob_vs_t").filename
        if args.beta:
            ProbabilityPlot().plot_p_of_t_multi(events, filename, beta=args.beta, tau=args.tau)
        else:
            ProbabilityPlot().plot_p_of_t_multi(events, filename)
    if args.plot_waiting_times:
        filename = SaveFile(args.event_file, "waiting_times").filename
        ProbabilityPlot().plot_dt_vs_t0(events, filename)
    if args.field_plot:
        element_ids = model.getElementIDs()
        ele_slips = {}
        for ele_id in element_ids:
            ele_slips[ele_id] = 5.0      
            sys.stdout.write("Loading {} elements...".format(len(element_ids))) 
        if args.colorbar_max: 
            cbar_max = args.colorbar_max
        else:
            cbar_max = None
        FP = FieldPlotter(model, element_slips=ele_slips, cbar_max=cbar_max)
        FP.plot_field(args.field_type, output_file=args.field_savefile, cutoff=1000)

    # Generate stress plots
    if args.stress_elements:
# TODO: check that stress_set is valid
        StressHistoryPlot().plot(stress_set, args.stress_elements)

    # Validate data if requested
    err = False
    if args.validate_slip_sum:
        mean_slip = sum(events.event_mean_slip())
        if abs(mean_slip-args.validate_slip_sum)/args.validate_slip_sum > 0.01: err = True
        print("Calculated mean slip:", mean_slip, "vs. expected:", args.validate_slip_sum)

    if args.validate_mean_interevent:
        ie_times = events.interevent_times()
        mean_ie = sum(ie_times)/len(ie_times)
        if abs(mean_ie-args.mean_interevent)/args.mean_interevent > 0.02: err = True
        print("Calculated mean interevent:", mean_interevent, "vs. expected:", args.mean_interevent)

    if err: exit(1)

