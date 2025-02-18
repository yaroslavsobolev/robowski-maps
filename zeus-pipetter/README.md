## Reagent pipetting platform ##

## Introduction ##
These scripts control an automated platform that utilizes a pipetting module 
and an XY gantry to precisely transfer chemical reagents between designated containers.

#### Zeus.py ####
This script manages communication between the server/PC and the pipetting module via a CAN bus, 
sending control commands and configuring liquid class parameters.

#### gantry.py ####  
This script controls the XY gantry by an Arduino-running firmware: [GRBL](https://github.com/grbl/grbl).

#### breadboard.py ####
This script configures the objects on the breadboard, including
* **wellplates**
* **containers:** vials (2 mL), wells (200 microliter), tubes (1.5 mL), bottles (20 mL) and jars (100 mL).
* **racks** for tips: 50, 300, 1000 microliter tips
* **balance**
* **trash bin**

#### pipetter.py #### 
This script orchestrates the interaction between Zeus.py and gantry.py, 
ensuring seamless coordination between the pipetting module and the XY gantry 
for precise reagent transfer and system control."

#### planner.py #### 
The planner handles the execution of liquid transfers by 
processing a spreadsheet input, converting it into structured dataframes 
for each transfer. These transfer events are then interpreted as objects, 
which are subsequently sent to pipetter.py for precise execution.

#### main.py ####
The central script that oversees and manages the entire system, 
coordinating all modules to ensure seamless operation and  
precise execution of liquid transfers.
