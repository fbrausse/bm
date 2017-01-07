#version 120

uniform mat3x4 transform;

attribute vec2  in_position;
attribute float in_diameter;
attribute vec3  in_color;

varying vec3 ex_color;
varying vec2 ex_position;
varying float ex_radius;

void main(void)
{
	gl_Position  = transform * vec3(in_position, 1);
	gl_PointSize = in_diameter;

	ex_color = in_color;
	ex_position = in_position;
	ex_radius = in_diameter / 2;
}
