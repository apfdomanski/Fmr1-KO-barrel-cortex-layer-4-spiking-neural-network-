### Distorted Ex/In dynamics and pattern classification in a spiking model of Fmr1-KO cortical layer 4
---
![Cartoon Network](Screenshots/Picture7.png)

This MATLAB code accompanies the following BiorXiv Preprint: 

Domanski APF, Booker S, Wyllie DJA, Isaac JTR, Kind PC (2018)
"_Cellular and Synaptic Compensations Limit Circuit Disruption in Fmr1-KO Mouse but Fail to Prevent Deficits in Information Processing_"

BiorXiv preprint: https://doi.org/10.1101/403725

Author: Aleksander PF Domanski 2015-2019 University of Bristol, UK aleks.domanski@bristol.ac.uk

Copyright: (C) Aleksander PF Domanski 2019 University of Bristol, UK

## License: 
GNU General Public License version 2
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
http://www.gnu.org/copyleft/gpl.html

## Usage:
This code builds a spiking network model of early postnatal _Fmr1-KO_ cortical layer 4, provides stimulation and analyses single cell and population output.
Some features of this model:
- Leaky Integrate-and-Fire model neurons with Conductance-based synapses and realistic intrinsic properties
- Sparse, random connectivity between Ex and In neurons
- Realistic short-term depression at all synapses
- _All_ parameters tuned by experimental results from the above paper


`L4sim_DesignNetwork.m`

This function builds intrinsic parameters and a synaptic connectivity  matrix for a recurrent spiking neural network with external stimulation and synapse-specific short-term plasticity.

---


![Cartoon Network](Screenshots/Picture11.png) 




