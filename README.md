# acoustic-logic-gates
## **Modelling the acoustic logic gates by solving the wave equation over the prescribed geometries by Wang et al. (10.1038/s41598-019-44769-0)**

This project presents the numerical modelling of acoustic logic gates in accordance with Wang et al. (10.1038/s41598-019-44769-0). Wang et al. introduced OR, XOR, NOT, AND, XNOR, XOR and XAND gates based on manipulating the incoming acoustic waves with Helmholtz resonators. This project is only limited with OR and XOR gates, and requires further work on the XOR gate as the introduced geometry does not satisfy the output requirements within the prescribed frequency range.

The model is tested for an exemplary microfluidic channel, with rigid and stationary walls (bounce-back boundary condition) and periodicity in x-direction.

The animation below shows the transport of the droplets inside the microfluidic channel to perform the first coalescence and then splitting of droplets;
and the color bar denotes velocity of the base fluid, and therefore maximum value is present at the interface. For the simulation details, two droplets are placed at the beginning of the channels. For t < 1300, a flow $\Delta p > 0$ is introduced to the fluid at the entrances of two channels to allow the droplets to meet at the junction. As the coalescence is observed, the flow is resetted $\Delta p = 0$ to let the larger droplet undergo relaxation, for 1300 < t < 1800. After t > 1800, a reversed-flow $\Delta p < 0$ is introduced to demonstrate splitting the larger droplet into two parts.

![](https://github.com/lynspica/droplet-microfluidics-lbm/blob/main/figs/channel.gif)
