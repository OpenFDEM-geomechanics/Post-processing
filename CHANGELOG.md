# Current Version 13 - Dated April 14, 2021

# THINGS TO IMPROVE
1. Integrate Qi Matlab Code for "seismic tool" <+> Maybe work on returning the file is CSV for Matlab (V#07)? # Needs to be verified

# OUTPUTS
1. Stress/Strain
1. Stress/Volumetric Strain
1. Prints information about the model, mesh, and materials in the model
1. Takes screenshots at predefined frequencies and calculates the Rosette diagram for it
1. Makes the Rosette diagram for the Mesh and Initial DFN if it exists
1. No user input required, the script identifies everything
1. Produces animation of the screenshots with the rosette diagrams
1. Has a histogram showing
   1.1 ax1 = Histogram chart of broken Phase Interfaces. 
   1.1 ax2 = Stacked Bar Chart different broken Phase Interfaces based on Time Step. 
   1.1 ax3 = Stacked Bar Chart of Crack Type (IntER/IntRA granular cracks)
   1.1 ax4 = Stacked Bar Chart for Fracture Mode Type (Tensile Dominant / Shear Dominant / Mixed Mode)
1. CSV => Clusters/consolidates the seismic data is case you want to do some post-processing data
1. CSV => Clusters/consolidates the broken joints in case you wish to do more post-processing

# UPDATES
# V01
1.Seismic Clustering Added.
1.Graphical representation of crack type.

# V02
1.Graphical representation of crack type - Enhanced.
1.Input routine enhanced.
1.Graphical display detaches as a Child Window.

# V03
1.Added Graphical Display (Failure Mode Type).
1.Code Enhancements & General Cleanup.
1.Percentage Materials displayed (takes into account the platens).

#V04
1.Locate Broken Joint along DFN.
1.Added Failure_Mode dictionary & Changed algorithm.
1.Added crack type dictionary.
1.Time stamp of each output streamlined.
1.Color scale for ax3 and ax4 amended.
1.Range of Color Scheme rectified.
1.dfn_line_list has set to allow only unique values.
1.pvpython cycler error highlighted.

#VO5
1.seismicclustering function updated:
        # 1) To read the #11 (12) ID in the reader.
        # 2) Return the time step as an integer so it can do the max/min operation.
1.Allow error free execution if seismic data is not available.

#V06
1. Visualization code runs separately from the rest of the code to overcome pvpython/python limitation on Glass computers.

#V07
1. Optimized Timer Display.
1. Optimized the pvpython/python Visualization problem.
1. Better Validations and Terminates if no broken joints.
1. Checks for missing modules prior to execution to save on processing time.

#V08
1. Table name inserted instead of Array Number.
1. Identifies 2D or 3D Simulation.
1. Based on 2D or 3D, skip variable added to iterate over points of the broken joints.
1. Batch processing added.
1. Return 'Area' instead of length.

#V09
1. Returns the angle of the crack (orientation).
1. WindRose (Radar) Diagram of the crack orientation binned based on failuremode. Bins can be changed to parameter of choice.

#V11
1. os.walk() added to collect all entire directory tree.
1. gather information on the file being processed.
1. create a "post_processing" folder that saves all the files to it.
1. Added outputs for the animation of the fractures.
1. Updated the history.csv files.
1. Attempt to engage all processors.
1. Fixed memory leak with pvpython.
1. Changing Algorithm.

#V12
1. Identifies platen ID without user input.
1. BD and UCS are now identified.
1. Points on boundary of BD and Area calculated correctly.
1. Debugged the memory loss for the animation.
1. Updated Algorithm.
1. Added consistency to terminology usage:
    1. Mode 1, Mode 2, Mixed Mode = Tensile Dominant, Shear Dominant, Mixed Mode (I - II).
    1. Mode 1 - Failure = [1.0 to 1.5].
    1. Mode 2 - Failure = (1.5 to 2.0].
    1. Mixed Mode = 3.0.
1. IntER and IntRA granular cracks properly identified.
1. Uses broken joints to trace cracks and the DFN is used for verification.
1. Overall code improvements and comments.
1. Added Stress/Strain (UCS / BD) Stress/Volumetric Strain (UCS).
1. Updated BD Analysis.
1. Additional checks on the Platen ID.

#V13
1. Separate pvpython and python modules.
   
#V14
1. Revamped the strain gauge option.
1. Added more illustrations to include the seismic data on the stress-strain.
1. Understands the failure modes from 1 to 8 based on the new Irazu Manual.
1. Runs on venv 
1. Obseleted the specimen center by using modelbounds()