# sembdip3d
Computes dip and coherency at each sample in a 3D volume.


```
********************************************************************************
Description: 
  Computes the dip and coherancy at each sample via semblance. 
  Program expects the following command line:  
  sembDip <in.sgy><dipI.sgy><dipX.sgy><sem.sgy><Num><maxDip><win><thresh>
  <reso><iOrigin><xOrigin><iMax><xMax>
 
Inputs: 
  in.sgy: Input segy filename. 
  dipI.sgy: Output Inline dip sgy filename.
  dipX.sgy: Output Xline dip sgy filename.
  sem.sgy: Output semblance  sgy filename.
  Num: Number of traces for semblance given as  
       distance from trace being analysed. 
       ex. numTr = 1 uses 9 traces. 3 by 3 block.
		o--o--o
	        |  |  |
		o--x--o
		|  |  |
		o--o--o
  maxDip: Maximum dip in ms/trace.
  win: Height of window in samples (along time dimension) for semblance analysis.
  thresh: minimum value for semblance - dips zerored if below this value.
  reso: resolution in points/samp. Ex. reso = 1 calculates dip to within one 
	point per sample per trace.  reso = 2 calcs to within 2 points/samp.
  iOrigin: Starting inline number.
  xOrigin: Starting xnline number.
  iMax: Max iline. 
  xMax: Max xline.
  
********************************************************************************
```
