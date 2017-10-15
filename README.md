# Fanuc Robot Point Translator
This program converts robot points between 6-axis Fanuc robots with different configurations/environments.

Copying robot programs from one robot to another is easy if they have the same configuration and the workpiece is in the exact same location with respect to each robot. However, if the workpiece location changes from one robot to another (for example, if robot A has to reach one inch farther than robot B to touch the workpiece), it is not possible to directly copy programs between the robots. In addition, some robots may have different configurations or end-of-arm tooling which prevents a direct program copy.

This program is my solution to the above problem. Using a static center point (e.g. the center of a turntable) as a reference point, each robot point can be re-defined based on its location and orientation with respect to this reference point. As long as the workpiece location remains the same with respect to the center point, this method can be used to copy programs from robots with very different environments.

## Files
The "program_translator.py" file is the main file, which calls the other files and runs the whole program. The only files not used during normal program execution are "utool_converter.py" and "test_point_converter.py". These files have many duplicate functions from "point_converter.py" and are used for testing.

## Notes
This program is designed for Fanuc robots. A similar approach could be used for other types of robots, but some differences must be accounted for that are not included in this code. Additionally, this program was developed for a specific application (thermal spray) and I am using this repository mostly as a backup, so I do not intend modify the code for more general applications. This is also the reason why some variable names may seem unusual to other Fanuc robot users; for example, "booth" and "gun" are application-specific terms.

Finally, for this program to work, information must be known about the each robot's configuration and environment. I use an Excel document that contains this information, which includes: robot utools and uframes, reference point coordinates, end-of-arm tooling orientations, and some minor additional environment information.