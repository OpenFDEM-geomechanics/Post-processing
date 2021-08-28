# OpenFDEM Objectives and Future plans

## DONE for UCS

### Platen tracking (Should track any body of intact elements)
1. Force (x y z directions)
2. Displacement (x y z directions)
3. Velocity xyz
4. Contact area (should be separate function?)

### Strain gauge analysis
1. input: Location, direction, length (DIRECTION And Length <=> PENDING: Center Location)
2. output: Strain gauge simulator?

### Model processing utilities
1. Stress/strain diagram (COMPLETED)
2. Volumetric strain (Lateral Strain (COMPLETED))
3. Elastic mod. (COMPLETED)
4. Plotters (COMPLETED)
5. Extract individual cells data into a panda.DataFrame

## TODO
1. Ability to define the center point of the SG.
2. Function to extract stress/strain for the element
 

### Need full OpenFDEM outputs for UCS/BD

### BD processing

### Fracture/fragment characterization
1. Clustering
2. Fracture timeline
3. Fragment element groups

### Seismic processing

### Shear test
1. Contact locations
1. Monitoring points for seismicity asperity breakage

## Documentation
1. website
2. combine with preprocessing and openfdem processing?
3. creating a manual PDF from Sphinx?

## Testing scripts?

## Pre-processing UI
1. Mesh from Geo file?