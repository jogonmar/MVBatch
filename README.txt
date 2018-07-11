------------------------------------------------------------------------------------------------------
MVBATCH toolbox for bilinear batch process modeling
------------------------------------------------------------------------------------------------------

Contents

1. Key features
2. Requirements
3. Installation
4. Getting started
5. License
6. Citing
7. Feedback 

------------------------------------------------------------------------------------------------------
1. Key features
------------------------------------------------------------------------------------------------------

The MVBatch toolbox is a graphical user-friendly interface for process understanding, troubleshooting
and monitoring. The analysis is performed by following the main steps of the bilinear modeling cycle: 
from the missing data imputation and synchronization of batch trajectories, to the calibration of latent
variable models based on Principal Component Analysis (PCA) and monitoring of historical batches. 
In the market, there exist commercial software packages addressing these modeling steps from different 
perspectives, such as ProMV Batch by ProSensus, Unscrambler X Batch Modeling by CAMO and SIMCA by MSK. 
The intention of this toolbox is not to replace these standard tools, but proving the scientific 
community with the most recent analytical solutions that overcome limitations of methods already 
available. Besides, the MVBatch Toolbox has been designed and implemented as a free software package 
to enable the generation of new ideas in the scientific community. Contributions, additions and 
improvements from the chemoetric comminuty is encouraged. 

The MVBatch toolbox integrates the most recent developments with the well-known techniques for batch 
missing data imputation, synchronization, calibration and monitoring. 


The package includes the following functionality:

(1.1) Missing data imputation
- TSR model building: Chemometrics and Intelligent Laboratory Systems, 146:77-88, 2015.
- TSR model exploitation: Journal of Chemometrics, 19: 439-447, 2005.

(1.2) Batch data synchronization
- Indicator Variable: AIChE Journal, 40:1361{1375, 1994.
- Dynamic Time Warping (DTW-Kassidas): AIChE Journal, 44:864{875, 1998.
- Dynamic Time Warping (DTW-Ramaker): Analytica Chimica Acta, 498:133-153, 2003.
- Dynamic Time Warping (variable constraints):  PhD thesis,  2015, 10.4995/Thesis/10251/55684.
- Relaxed Greedy Time Warping: Chemometrics and Intelligent Laboratory Systems, 105:195-206, 2011 .
- Multi-synchro: Journal of Chemometrics, 28:462-475, 2014.

(1.3) Pre-processing
- Trajectory centering
- Trajectory centering and scaling  
- Trajectory centering and variable scaling
- variable centering
- variable centering and scaling

(1.4) Cross-validation.
- CKF: Journal of Chemometrics, 29(8):467-478, 2015
- EKF: Chemometrics and Intelligent Laboratory Systems 131:37-50, 2014.

(1.5) Transformation of the three-way array into a two-way array (calibration based on Principal 
Component Analysis)
- Variable-wise unfolding: Chemometrics and Intelligent Laboratory Systems, 44:331-340, 1998.
- Batch Dynamic unfolding: Chemical Engineering Science, 57:63-75, 2002.
- Observation-wise unfolding: Comprehensive Chemometrics - Chemical and Biochemical Data Analysis, 
  Elsevier:Oxford, 1:163{195, 2009.
- Multi-phase: Journal of Chemometrics, 22:632-643, 2008.

(1.6) Exploratory data analysis
- MEDA toolbox: Chemometrics and Intelligent Laboratory Systems, 143:49–57, 2015.


The MVBatch toolbox for MATLAB is freely available for download at
https://github.com/jogonmar/MVBatch/releases/tag/v1.0. 

------------------------------------------------------------------------------------------------------
2. Requirements
------------------------------------------------------------------------------------------------------

A functional installation of MATLAB 7.0 or later is required. Compatibility has been checked until
Matlab 9.1. To get full functionality, the MATLAB statistics toolbox and the Multivariate Exploratory 
Data Analysis (MEDA) Toolbox. The latter is an open-source software used for visualization and exploratory
data analysis through the latent structures. A stable version of the MEDA toolbox can be download at 
https://github.com/josecamachop/MEDA-Toolbox/releases.

Note that the MVBatch toolbox has been tested in Matlab running under Window Vista and Windows 7 Operative 
Systems. Due to potential incompatibilities with Windows 10 API, the visualization of the MVBatch toolbox 
might be affected. In case of API issues, the authors suggest running the code via the console of Matlab.  

------------------------------------------------------------------------------------------------------
3. Installation
------------------------------------------------------------------------------------------------------

a) Unzip the MATLAB source files (*.m and *.mat) of the MVBatch toolbox into any directory of choice 
<directory_path>

. 

b) Unzip the Matlab files of the MEDA toolbox into the directory <directory_path>.
 
c) Add to the MATLAB path the following directories (use command addpath, e.g. addpath '<path>'):
		
	- <directory_path>
		
	- <directory_path>/examples	
	- <directory_path>/GUI
	- <directory_path>/MEDA-Toolbox
	- <directory_path>/modeling
	- <directory_path>/synchronization

   Another alternative is to add folders and subfolders via the Matlab user interface: Home-> Set Path -> 
   -> Add with Subfolders.

------------------------------------------------------------------------------------------------------
4. Getting started
------------------------------------------------------------------------------------------------------

There are two ways to work with the MVbatch Toolbox: using the graphical user interface (GUI) for starting
users and using the command line for expert users.The main GUI is self-explanatory and follows the modeling
cycle of batch processes: alignment, preprocessing, transformation of the three-way array into a two-way 
array and calibration, and monitoring. It incorporates help messages to guide its use. To launch the main 
GUI, type "MVBatch" in the command line of Matlab after the installation, as explained in Section 3. 

Once the MVBatch user interface is launched, the user will have to load a .mat file which contains:
 
i)   data: (Ix1) cell array containing the 2-way arrays (KixJ) cotaining the batch trajectories
ii)  batchNames: (Ix1) string cell array containing the Batch ID. 
iii) varNames: (Jx1) string cell array containing the tagnames.

Examples of this data structure can be found in folder "Examples". Note that a different naming of these structures
will lead to an error prompted in the command window of Matlab.

The Matlab function provide includes helping information, that can be seen by typing "help <command>" in 
the command line of Matlab. This helping information includes examples on how to execute the corresponfing 
function. 

The package includes a demonstration of the functionality available in the toolbox, which is intended to be 
a basic guideline for future use. The demonstration is based on a simulated data set of the Saccharomyces 
Cerevisiae cultivation, which contains batches affected by four different types of asynchronisms. The description
of the bilinear modeling of this data set can be found in Section 4 (Case study: Saccharomyces cerevisiae) of
the scientific article. 

------------------------------------------------------------------------------------------------------
5. License
------------------------------------------------------------------------------------------------------

The MVBatch toolbox is free software: you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public 
License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, 
see <http://www.gnu.org/licenses/>.


------------------------------------------------------------------------------------------------------
6. Citing 
------------------------------------------------------------------------------------------------------

Please, acknowledge the use of this software by referring it: 

"Gonzalez-Martinez, J.M., Camacho, J., Ferrer, A. MVBatch: a Matlab toolbox for batch process modeling and monitoring. Submitted to Chemometrics and Intelligent Laboratory Systems, 
available at https://github.com/jogonmar/MVBatch/releases/tag/v1.0" 

------------------------------------------------------------------------------------------------------
7. Feedback
------------------------------------------------------------------------------------------------------

Please send comments to José María González Martínez <jogonmar@gmail.com> or
José Camacho Páez <josecamacho@go.ugr.es>
