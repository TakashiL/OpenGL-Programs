# Robot Arm Animation

## How to Run the Game:
* Use command "make" to create the executable file.
* Use command "./myrobot" or "./myrobot old_x old_y old_z new_x new_y new_z <-tv|-sv>" to run the program.
* User can press 'q', 'Q' or 'esc' key to terminate the program.

## Viewing mode:
* Use command "./myrobot" to run the viewing mode.
* User can press key '1' to choose the base, press key '2' to choose the lower arm, or press key '3' to choose the upper arm.
* User can press "right" and "left" key or using mouse clicks to rotate the part of the robot that has been chosen.
* User can press "t" key to change to the top view, or press "s" key to change back to the side view.

## Fetching mode:
* User command "./myrobot old_x old_y old_z new_x new_y new_z <-tv | -sv>" to run the fetching mode.
* (old_x, old_y, old_z) is the initial position for the sphere.
* (new_x, new_y, new_z) is the desination for the sphere.
* "-tv" for top view, "-sv" for side view.
* User can press "t" key to change to the top view, or press "s" key to change back to the side view while program is running.
