# ParticleImaging
This software generates an image of a spherical particle given a set of parameters.
Solves scattering of light in a sphere with geometrical optical approximation.
## User guide
To generate an example image run `example.py`. This creates three images side by side from three cameras capturing the particle.

To change the output of the images, the following parameters can be changed:

`n_m` - Refractive index of the medium around the particle (ie. air).\
`n_p` - Refractive index of particle medium (ie. water).\
`r` - Particle radius [m].\
`rho` - Ray density per square meter.\
`particle_center_x` - Particle center x-coordinates.\
`particle_center_y` - Particle center y-coordinates.\
`particle_center_z` - Particle center z-coordinates.\
`pol_X` - Polarization x-position.\
`pol_Y` - Polarization y-position.\
`pol_Z` - Polarization z-position.\
`pol_Vx` - Polarization x-direction.\
`pol_Vy` - Polarization y-direction.\
`pol_Vz` - Polarization z-direction.\
`scattering_max_number_of_iterations` - Maximum number of iterations the rays are scattered within the particle.\
`noise_intensity` - Perlin noise with an intensity is generated.\
`lens_size` - Camera lens size [m]. This is equal to the image width. Also assuming quadratic lens/image.\
`resolution` - Resolution of the generated image. Assuming quadratic image. Number of pixels in image becomes resolution*resolution.\
`n_of_cameras` - Number of cameras to capture images of droplet.\
`x_rotations` - Array with the rotation (radians around the x-axis) required to flip the coordinates from one camera to the other.\
`y_rotations` - Array with the rotation (radians around the y-axis) required to flip the coordinates from one camera to the other.\
