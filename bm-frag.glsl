#version 120

uniform mat3x4 transform;

varying vec3 ex_color;
varying vec2 ex_position;
varying float ex_radius;

void main(void)
{
	vec3 color = ex_color;
	float dist = 2*distance(gl_PointCoord,vec2(.5,.5));
	if (dist > 1)
		discard;
	float delta = 1/ex_radius;
	float alpha = 1-smoothstep(1-delta,1,dist);
	color *= alpha;
	gl_FragColor = vec4(color, 1);
}
