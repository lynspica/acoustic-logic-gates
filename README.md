# acoustic-logic-gates
## **Modelling the acoustic logic gates by solving the wave equation over the prescribed geometries by Wang et al. (10.1038/s41598-019-44769-0)**

This project presents the numerical modelling of acoustic logic gates in accordance with Wang et al. (10.1038/s41598-019-44769-0). Wang et al. introduced OR, XOR, NOT, AND, XNOR, XOR and XAND gates based on manipulating the incoming acoustic waves with Helmholtz resonators and used COMSOL to validate their results. This project is only limited with OR and XOR gates, and requires further work on the XOR gate as the introduced geometry does not satisfy the output requirements within the prescribed frequency range.

The figure below shows the OR and XOR logic gates geometries constructed in MATLAB. The upper geometry belongs to the OR gate, and the lower gate belongs to the XOR gate. Helmholtz resonators and input and output regions are depicted in the figure.

![](https://github.com/lynspica/acoustic-logic-gates/blob/main/figs/geometries.png)

The animations below show the propogation of an acoustic wave (sinusoidal signal with arbitrary frequency for demonstration) over the XOR gates, where the upper simulation corresponds to Input = {1,0} and the below corresponds to Input = {0,1}.  

![](https://github.com/lynspica/acoustic-logic-gates/blob/main/figs/XOR_input10prefactor_19.gif)

![](https://github.com/lynspica/acoustic-logic-gates/blob/main/figs/XOR_input01prefactor_19.gif)

The figure below shows the average pressure obtained at the outlet, over a range of input frequencies. For the OR gate, the numerical model is in accordance with the expected behavior of OR gates, while the XOR gate with Input = {0,1} does not show the expected behavior of the XOR gate.

![](https://github.com/lynspica/acoustic-logic-gates/blob/main/figs/results.png)
