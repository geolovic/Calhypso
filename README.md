# CalHypso
QGIS toolbox to create hypsometric curves from a 
Digital Elevation Model and a polygon shapefile with basins. 

## Instalation

To install the toolbox, copy the **Calhypso.py** and **Calhypso.py.help** files in the script folder:
- **Users/user_name/.qgis2/processing/scripts**  [MAC] (Hidden folder)
- **C:\Users\user_name\.qgis2/processing/scripts**   [Windows]

An alternative to install the toolbox is using the **"Add script from file"** tool located in the Processing Toolbox > Scripts > Tools. (This way does not include the help file)

## Usage
Launch the QGIS Toolbox and fill in the required parameters:

![Calhypso Toolbox](https://geolovic.github.io/Calhypso/images/CalypsoToolbox.jpg)

1. Digital Elevation Model (DEM) to take elevations. 
2. Polygon shapefile with basins
3. Id field in the shapefile. It MUST be an integer field identifying the different basins. Avoid use the zero as an id number.
4. Name Field in the shapefile to take basin labels

Once the tool have processed all the basins, it will show a graphic 
window with the hipsometric curves. The graphic window can show all the curves
together (G1) or individual basins (G2).

![Calhypso Toolbox](https://geolovic.github.io/Calhypso/images/Graphic_types.jpg)

In this graphic window we can use the following keys:

* **A** : Shows all the hipsometric curves in a single graphic (G1)
* **LEFT / RIGHT** : Navigates throught single hypsometric curves (G2)
* **M** : Shows/hides alternative hypsometric metrics in G2 (Kurtois, 
Skewness, Density Kurtosis and Density Skewness, see 
**Harlin, 1989** and **Pérez-Peña et al., 2009**)
* **D** : Shows/hides the legend in G1
* **UP** : Selects the HI color mode. Curve colors for G1 are selected by their HI 
values according the used ramp color
* **DOWN**: Selects the order color mode. Curve colores for G1 are selected by 
their order (by useing the Id field) 
according the used ramp color
* **C** : Change the color ramp for G1. There are 8 different diverging color 
ramps available.
* **X** : Save all the profiles in a specific format. It creates a folder called *Calhypso_output_files* in the same location that
the basin shapefile to put all figures there. Format has to be selected with keys 1-3. 
single graphic per curve. If alternative metrics are on, they would be printed. 
* **1 - 3** - Keys to select formats (only active if **X** was pressed first) 
* **E** - Export all curves in text files. It creates a folder called *Calhypso_output_files* in the same location that
the basin shapefile to put all figures there. It also generated a file named **"curve_moments.txt"** with all the hypsometric metrics values.
