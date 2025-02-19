# Beachball Focal Mechanism Visualization

This project computes and visualizes earthquake focal mechanisms (commonly called "beachballs") using strike, dip, and rake parameters or a moment tensor. It provides several visualization methods including a 3D plot, a 2D (Wulff) projection, and an ASCII art output.

## Features

- **Focal Mechanism Calculation:** Converts fault plane parameters to a moment tensor and derives auxiliary vectors (tension, pressure, and null directions).
- **3D Visualization:** Displays the focal mechanism on a sphere with annotated P, B, and T axes.
- **2D (Wulff) Projection:** Projects the spherical data onto a 2D plane.
- **ASCII Plot:** Generates a simple ASCII art representation for terminal viewing.

## Requirements

- Python >= 3.11
- [NumPy]
- [Matplotlib]

## Installation

You can install the required dependencies using pip:

```bash
pip install numpy matplotlib
```

## Usage
Example 1
```python
bb = Beachball(strike=30, dip=40, rake=50)
bb.plot_beachball_3d()
bb.plot_beachball_2d()
bb.plot_beachball_ascii()
```

Example 2
```python
bb = Beachball(mt=[1,2,3,4,5,6]) # mt in ned axis
bb.plot_beachball_3d()
bb.plot_beachball_2d()
bb.plot_beachball_ascii()
```

Example 3
```bash
python beachball.py 30 40 50

# output
            ................            
       ......              ......       
     ...                      ++++.     
   ...                         ++++++   
 ...                            ++ ++.. 
..               + ++++++++++++++++ +++.
.        + ++++++++++++++++++++        .
.     ++++++++++++++++++++++++         .
.  +++++++++++++++++++++++++++         .
.+++++++++++++++++++++++++++          ..
 .++++++++++++++++++++++++          ... 
   +.++ +++++++++++++++           ...   
     .+.++ ++++++++             ...     
       ..+.+++             ......       
            ................
```
