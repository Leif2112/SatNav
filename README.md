# SatNav

This report is written in support of the MATLAB code developed for satellite orbital mechanics. In essence, this report should lead the reader through the implications of the two-body problem, demonstrate how they affect a satellite’s motion in space, and provide clear derivation of the underlying physics. It should also discuss the various methods used to determine and predict a satellite’s orientation in space, considering a torque-free context, all the while studying their relative merits. Ultimately, the objective of this paper is to explore the fundamentals of space dynamics, specifically, the principles of orbital and attitude mechanics. These hold great significance to spacecraft deployment and navigation in space.

## Two Body Problem
In order to derive simple analytical results for the motion of a spacecraft about the Earth, it is assumed that any external forces from additional celestial bodies are disregarded - no appreciable force is exerted on the spacecraft from a third body.
Likewise, any additional effects of orbit perturbations (i.e., solar wind pressure, atmospheric drag, Earth’s oblateness etc)are neglected, further simplifying the model. Under these assumptions, it can be stated that the motion of a satellite is purely due to the gravitational interaction with Earth.

The following orbital parameters were considered during the development of the software:
𝑎 = 7151.6 km
𝑒 = 0.0008
𝑖 = 98.39 deg
Ω = 10.0 deg
𝜔 = 233.0 deg
𝑀0 = 127.0 deg
