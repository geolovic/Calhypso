# -*- coding: iso-8859-15 -*-
#
#  CalHypso.py
#
#  Copyright (C) 2017  J. Vicente Perez, Universidad de Granada
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

#  For additional information, contact to:
#  Jose Vicente Perez Pena
#  Dpto. Geodinamica-Universidad de Granada
#  18071 Granada, Spain
#  vperez@ugr.es // geolovic@gmail.com

#  Version: 1.0
#  November 10th, 2017

#  Last modified November 10th, 2017

import gdal
import ogr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import matplotlib.cm as cmap


COLOR_RAMPS = {"BrBG":cmap.BrBG, "bwr":cmap.bwr,"coolwarm":cmap.coolwarm, "PuOr":cmap.PuOr,
               "PRGn":cmap.PRGn, "PiYG":cmap.PiYG, "RdYlBu":cmap.RdYlBu, "RdYlGn":cmap.RdYlGn,
               "RdBu":cmap.RdBu, "RdGy":cmap.RdGy, "seismic":cmap.seismic, "Spectral":cmap.Spectral}

DIRECTIONS = {"LEFT": -1, "RIGHT": 1}
FILETYPES = {1: ".eps", 2: ".pdf", 3: ".png"}


# DEBUG PARAMETERS
# ================
basin_shapefile = "../SampleData/gisdata/basins.shp"
dem = "../SampleData/gisdata/dem40.tif"
attr_field = "id"
name_field = "name"


def main(dem, basin, attr_field, name_field):
    out_dir = os.path.dirname(basin_shapefile) + "/calhypso_output_files/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    curves = get_hypsometric_curves(dem, basin, attr_field, name_field)
    fig = plt.figure()
    return CalHypsoApp(curves, fig, out_dir)
    

def get_hypsometric_curves(dem, basin_shapefile, attr_field, name_field=""):
    """
    This function extracts hipsometric curves from a basin shapefile and a Digital Elevation Model
    :param dem: Digital Elevation Model
    :param basin_shapefile: Polygon shapefile with the basins
    :param attr_field: Integer field that defines each shapefile
    :param name_field: Field with names to label curves
    :return: List of CurvaHipsometrica objects
    """
    # Open the DEM and get attributes
    dem_raster = gdal.Open(dem)
    dem_band = dem_raster.GetRasterBand(1)
    dem_arr = dem_band.ReadAsArray()
    geot = dem_raster.GetGeoTransform()
    proj = dem_raster.GetProjection()
    xsize = dem_raster.RasterXSize
    ysize = dem_raster.RasterYSize

    # Clean-up
    dem_band = None
    dem_raster = None
    
    # Create a new blank "in-memory" raster (output_ras)
    driver = gdal.GetDriverByName("MEM")
    output_ras = driver.Create("", xsize, ysize, 1, gdal.GDT_Int16)
    output_ras.SetGeoTransform(geot)
    output_ras.SetProjection(proj)
    output_ras.GetRasterBand(1).SetNoDataValue(0)
    
    # Open the basin_shapefile and get layer
    dataset = ogr.Open(basin_shapefile)
    layer = dataset.GetLayer(0)
    
    # Rasterize vector layer by attribute (ensure that field exists and is type integer)
    layerdef = layer.GetLayerDefn()
    fields = [layerdef.GetFieldDefn(n).GetName() for n in range(layerdef.GetFieldCount())]
    field_type = layerdef.GetFieldDefn(layerdef.GetFieldIndex(attr_field)).GetType()
    if attr_field in fields and (field_type == 0 or field_type == 12):
        gdal.RasterizeLayer(output_ras, [1], layer, options=["ATTRIBUTE=%s" % attr_field])
        banda = output_ras.GetRasterBand(1)
        basin_arr = banda.ReadAsArray()
    else:
        return []
    
    # Get basins names
    names = {}
    if name_field in fields:
        for feat in layer:
            names[feat.GetField(attr_field)] = str(feat.GetField(name_field))
    else:
        for feat in layer:
            names[feat.GetField(attr_field)] = str(feat.GetField(attr_field))

    # Clean-up Memory
    banda = None
    output_ras = None
    layer = None
    dataset = None
    
    # Get unique values of basin_arr except for "0"
    basin_values = np.unique(basin_arr)
    basin_values.sort()
    if 0 in basin_values:
        basin_values = basin_values[1:]

    curvas = []
    for value in basin_values:
        elev_values = dem_arr[basin_arr == value]
        # Avoid small basins (probably errors) and flat ones min == max
        if len(elev_values) > 10 and elev_values.max() > elev_values.min():
            curvas.append(CurvaHypsometrica(elev_values, names[value]))
    
    return curvas

    
class CurvaHypsometrica:
    """
    Class that define a simple hypsometric curve.
    """
    def __init__(self, elev_values, name=""):
        elev_values = np.sort(elev_values)[::-1]
        hh = (elev_values - elev_values[-1])/float(elev_values[0]-elev_values[-1])
        aa = (np.arange(elev_values.size) + 1)/float(elev_values.size)
        max_h = np.max(elev_values)
        min_h = np.min(elev_values)
        mean_h = np.mean(elev_values)
        
        # simplify elevation and area data to 10 bits (1024 values)
        x = np.linspace(0, 1, 1024)
        y = np.interp(x, aa, hh)
        self._data = np.array((x, y)).T

        # get hi values and name
        self._HI = (mean_h - min_h) / float(max_h - min_h)
        self._HI2 = self._calculate_hi2()
        self._name = name
        
        # get statistical moments
        self.moments = self._get_moments()
        
    def _get_moments(self):
        """
        Get statistical moments given the 4 coefficients of the curve polynomial regression

        This code has been translated from a FORTRAN code developed by Harlin (1979)
        Harlin, J. M. (1978). Statistical Moments of the Hypsometric Curve and Its Density
        Function. Mathematical Geology, Vol. 10 (1), 59-72.
        """
        coeficients = np.polyfit(self._data[:, 0], self._data[:, 1], 3)[::-1]
        moments = []
        ia = 0
        m10 = 0
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        for idx, coef in enumerate(coeficients):
            ia += coef / (idx + 1)
            m10 += coef / (idx + 2)
            tmp1 += coef / (idx + 3)
            tmp2 += coef / (idx + 4)
            tmp3 += coef / (idx + 5)
       
        tmp1 = tmp1/ia
        tmp2 = tmp2/ia
        tmp3 = tmp3/ia
        moments.append(ia)
        
        m10 = m10 / ia
        m20 = tmp1 - (m10**2)
        m30 = tmp2 - (3 * m10 * tmp1) + (2 * m10**3)
        m40 = tmp3 - (4 * m10 * tmp2) + (6 * (m10**2) * tmp1) - (3 * m10**4)
        moments.append(m30 / (m20**1.5))
        moments.append(m40 / (m20**2))
        
        demon = 0
        for coef in coeficients[1:]:
            demon += coef
        
        tmp1 = 0
        tmp2 = 0
        ex = 0
        ex4 = 0
        for n in range(3):
            ex += (n+1) * coeficients[n+1] / (n+2)
            ex4 += (n+1) * coeficients[n+1] / (n+5)
            tmp1 = tmp1 + (n + 1) * coeficients[n + 1] / (n + 3)
            tmp2 = tmp2 + (n + 1) * coeficients[n + 1] / (n + 4)
        
        ex = ex / demon
        tmp1 = tmp1 / demon
        tmp2 = tmp2 / demon
        ex2 = tmp1 - ex**2
        ex3 = tmp2 - 3.0 * ex * tmp1 + 2.0 * ex**3
        tk = ex2**0.5
        ts = tk**3
        sdk = ex3 / ts
        ex4 = ex4 / demon - 4.0 * ex * tmp2 + 6.0 * ex**2 * tmp1 - 3.0 * ex**4
        krd = ex4 / ex2**2 
        
        moments.append(sdk)
        moments.append(krd)
        
        return moments

    def _calculate_hi2(self):
        """
        Calculates HI by integrating the hypsometric curve
        """
        area_accum = 0
        for n in range(len(self._data[:, 0])-2):
            a1 = self._data[n, 0]
            a2 = self._data[n+1, 0]
            h1 = self._data[n, 1]
            h2 = self._data[n+1, 1]

            if h1 < h2:
                aux = h1
                h1 = h2
                h2 = aux

            area_accum += (a2-a1) * h2 + ((a2-a1) * (h1-h2))/2
        return area_accum

    def get_aa(self):
        """
        Return relative area values of the hypsometric curve
        """
        return self._data[:, 0]

    def get_hh(self):
        """
        Return relative elevation values of the hypsometric curve
        """
        return self._data[:, 1]

    def get_name(self):
        """
        Return the name of the hypsometric curve
        """
        return self._name

    def get_hi(self):
        """
        Return the hipsometric integral calculated by aproximation (hmed-hmin) / (hmax-hmin)
        """
        return self._HI

    def get_hi2(self):
        """
        Return the hipsometric integral calculated by integrating the curve
        """
        return self._HI2
    
    def get_ku(self):
        return self.moments[1]
    
    def get_sw(self):
        return self.moments[2]
    
    def get_dku(self):
        return self.moments[3]
    
    def get_dsw(self):
        return self.moments[4]


class CalHypsoApp:
    
    def __init__(self, curves, figure, out_dir = "", cm = "BrBG", HI_colors = True):
        
        # Drawing variables
        self.ax = figure.add_axes((0.05, 0.1, 0.8, 0.8))
        self.curves = curves
        self.n_curves = len(curves)
        self.active_id = 0
        self.gr_type = 1
        self.show_moments = False
        
        # Draws the legend if there are less than 20 curves
        if len(curves) < 20:
            self.draw_legend = True
        else:
            self.draw_legend = False
        
        # Color variables
        self._cmap = cm
        self._HI_colors = HI_colors
        # Get HIMax and HIMin to scale colors with HI values
        self.HIMax = 0
        self.HIMin = 1
        for curva in curves:
            if curva.get_hi() > self.HIMax:
                self.HIMax = curva.get_hi()
            if curva.get_hi() < self.HIMin:
                self.HIMin = curva.get_hi()
        
        self.out_dir = out_dir
        # Variables to save and log 
        self.saving = False
        self.info_text = figure.text(0.76, 0.9, "", size = 9, verticalalignment="top")
        
        # Connect the figure with a key_press_event
        self.b_cid = figure.canvas.mpl_connect("key_press_event", self.key_press)
        self.draw_all()

    def format_axe(self, ax, labels=True, grids=True):
        """
        This function formats the self.ax instance to represent hypsometric curves

        :param labels: (bool) If True, labels will be shown for X and Y axes
        :param grids: (bool) If True, grid lines will be shown
        :return: None
        """
        ax.clear()
        ax.set_aspect("equal")
        ax.set_xlim((0,1))
        ax.set_ylim((0,1))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        if grids:
            ax.grid(True, which='minor', axis="both", linestyle ="--", c="0.4", lw=0.25)
            ax.grid(True, which='major', axis="both", linestyle ="-", c="0.8", lw=0.75)
        if labels:
            ax.set_xlabel('Relative area (a/A)')
            ax.set_ylabel('Relative elevation (h/H)')
    
    def key_press(self, event):
        """
        Interactive key events for hypsometric curves:
            A >> Draws all the curves together. Color mode can be selected with the UP key
            LEFT / RIGHT >> Moves to the previous/next curve and represents it in a single graphic
            UP >> Changes the color sytle. Colors in the color ramp will be scaled by HI / Ordering
            C >> Changes the ramp color. There are 12 different diverging color ramps. 
                ['BrBG', 'bwr', 'coolwarm', 'PuOr', 'PRGn', 'PiYG', 
                'RdYlBu', 'RdYlGn', 'RdBu', 'RdGy', 'seismic', 'Spectral']
            D >> Displays/hides the legend when all the curves are represented together
            
        """
        self.info_text.set_text("")
        if event.key.upper() == "LEFT" or event.key.upper() == "RIGHT":
            self.active_id += DIRECTIONS[event.key.upper()]
            self.active_id = self.active_id % self.n_curves
            self.draw_active()
            
        elif event.key.upper() == "A":
            self.draw_all()
            
        elif event.key.upper() == "D" and self.gr_type == 0:
            if self.draw_legend == False:
                self.draw_legend = True
            else:
                self.draw_legend = False
            self.draw_all()
            
        elif event.key.upper() == "C" and self.gr_type == 0:
            self.next_color()
            
        elif event.key.upper() == "UP" and self.gr_type == 0:
            self.info_text.set_text("HI color mode")
            self._HI_colors = True
            self.draw_all()
            
        elif event.key.upper() == "DOWN" and self.gr_type == 0:
            self.info_text.set_text("Order color mode")
            self._HI_colors = False
            self.draw_all()
            
        elif event.key.upper() == "X":
            mensaje = "Prepared to save\nselect filetype:\n"
            mensaje += "1 > *.eps\n2 > *.pdf\n3 > *.png"
            self.info_text.set_text(mensaje)
            self.saving = True
            self.ax.figure.canvas.draw()
            
        elif event.key in ["1", "2", "3"] and self.saving:
            self.saving = False
            ftype = FILETYPES[int(event.key)]
            self.save_all(ftype)
            self.info_text.set_text("Successfully saved")
            self.ax.figure.canvas.draw()
            
        elif event.key.upper() == "E":
            self.export_all()
            self.info_text.set_text("Successfully exported")
            self.ax.figure.canvas.draw()
            
        elif event.key.upper() == "M" and self.gr_type == 1:
            if self.show_moments:
                self.show_moments = False
            else:
                self.show_moments = True
                
            self.draw_active()

    def draw_all(self):
        """
        This function combines and draws all the hypsometric curves in a single graphics.
        Curves will be represented with a color ramp according their HI value or their order. 
        """
        ax = self.ax
        self.format_axe(ax)
        cm = COLOR_RAMPS[self._cmap]
        
        for idx, curva in enumerate(self.curves):
            if self._HI_colors:
                curve_color = (curva.get_hi() - self.HIMin) / (self.HIMax - self.HIMin)
            else:
                curve_color = idx / len(self.curves)
                
            ax.plot(curva.get_aa(), curva.get_hh(), label=curva.get_name(), c=cm(curve_color))
        
        if self.draw_legend:
            ax.legend(loc="lower left", bbox_to_anchor=(1, 0))
        ax.set_title("Hypsometric curves")
        self.gr_type = 0
        ax.figure.canvas.draw()   
    
    def draw_active(self):
        """
        This function draws the active curve in the self.ax Axe object
        """
        ax = self.ax
        self.format_axe(ax)
        curva = self.curves[self.active_id]
        
        ax.set_title(curva.get_name())
        
        cadena = "HI: {0:.3f}\n".format(curva.get_hi())
        
        if self.show_moments:
            cadena += "KU: {0:.3f}\n".format(curva.get_ku())
            cadena += "SW: {0:.3f}\n".format(curva.get_sw())
            cadena += "DK: {0:.3f}\n".format(curva.get_dku())
            cadena += "DS: {0:.3f}\n".format(curva.get_dsw())
        
        ax.plot(curva.get_aa(), curva.get_hh(), c="navy")    
        ax.text(0.77, 0.97, cadena, size=11, verticalalignment="top")
        self.gr_type = 1
        ax.figure.canvas.draw()
        
            
    def next_color(self):
        """
        Changes to the next color ramp in COLOR_RAMPS
        """
        cc = list(COLOR_RAMPS.keys())
        idx = cc.index(self._cmap)
        idx += 1
        idx = idx % len(cc)
        self._cmap = cc[idx]
        self.draw_all()
             
    def save_all(self, filetype = ".eps"):
        """
        Save all the curves in the graphic (one figure per curve)
        Figures will be saved in the basedir directory
        """
        aux_fig = plt.figure()
        ax = aux_fig.add_axes((0.05, 0.1, 0.8, 0.8))

        for curva in self.curves:
            self.format_axe(ax)
            ax.plot(curva.get_aa(), curva.get_hh(), c="navy")
            ax.set_title(curva.get_name())
            cadena = "HI: {0:.3f}\n".format(curva.get_hi())
        
            if self.show_moments:
                cadena += "KU: {0:.3f}\n".format(curva.get_ku())
                cadena += "SW: {0:.3f}\n".format(curva.get_sw())
                cadena += "DK: {0:.3f}\n".format(curva.get_dku())
                cadena += "DS: {0:.3f}\n".format(curva.get_dsw())
          
            ax.text(0.77, 0.97, cadena, size=11, verticalalignment="top")
            ax.figure.canvas.draw()
            aux_fig.savefig(self.out_dir + "/" + curva.get_name() + filetype)
            
        plt.close(aux_fig)
        
    def export_all(self):
        """
        Export all the curves in the graphic
        """
        ff = open(self.out_dir + "/curve_moments.txt", "w")
        ff.write("Basin;HI;Kurtosis;Skewness;DKurtosis;DSkewness\n")
        for cc in self.curves:
            out_path = self.out_dir + "/" + cc.get_name() + ".txt"
            np.savetxt(out_path, cc._data)
            hi = cc.get_hi()
            ku = cc.get_ku()
            sw = cc.get_sw()
            dku = cc.get_dku()
            dsw = cc.get_dsw()
            momentos = "{0:.3f};{1:.3f};{2:.3f};{3:.3f};{4:.3f}".format(hi, ku, sw, dku, dsw)
            ff.write(cc.get_name() + ";" + momentos + "\n")
            
        ff.close()     
    
if __name__ == "__main__":
    # execute only if run as a script
    my_graph = main(dem, basin_shapefile, attr_field, name_field)
